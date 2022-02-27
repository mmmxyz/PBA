#include <iostream>
#include <cmath>
#include <cstdint>
#include <vector>
#include <deque>
#include <algorithm>
#include <random>
#include <omp.h>

#include "opengl/visualizer.hpp"
#include "opengl/vertarray.hpp"
#include "opengl/drawobject.hpp"
#include "opengl/renderer2d.hpp"

#include "utils/mathfunc/mathfunc.hpp"
#include "utils/mathfunc/polardecompose.hpp"

#include "utils/meshgenerator/meshgenerator.hpp"

#include "utils/fem/fem.hpp"

#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>

constexpr uint32_t FPS = 60;
constexpr float dt     = 1.0 / FPS;
//glfwSwapInterval(1); なので60FPS以下になる。
//これはモニタやGPUに依るかもしれない。

float mu     = 40;
float lambda = 80;
float rho    = 0.001;

bool rightwall = false;
float trans;

class ClothMesh {
    public:
	std::vector<fvec2> PositionList;
	std::vector<fvec2> RestPositionList;
	std::vector<fvec2> VelocityList;

	uint32_t vertsize;
	fvec2* vertdata;
	uint32_t* tridata;
	uint32_t trisize;
	uint32_t* edgedata;
	uint32_t edgesize;
	uint32_t* boundarydata;
	uint32_t boundarysize;

	linevertarray lva;
	trianglevertarray tva;

	std::vector<float> ElasticLamdalist;

	std::vector<fmat2> AList;
	std::vector<float> VList;
	std::vector<float> oList;

	const uint32_t N, M;

	int32_t MaterialInd = 0;

	float mass;

	ClothMesh(uint32_t N, uint32_t M, float lengthx, float lengthy, const fvec2& bias)
	    : N(N)
	    , M(M)
	    , lva()
	    , tva()
	{

		RectTriangle(N, M, lengthx, lengthy, &vertdata, vertsize, &tridata, trisize, &edgedata, edgesize, &boundarydata, boundarysize, bias);

		PositionList.resize(vertsize);
		RestPositionList.resize(vertsize);
		VelocityList.resize(vertsize);

		AList.resize(trisize / 3);
		VList.resize(trisize / 3);
		ElasticLamdalist.resize(trisize / 3);
		oList.resize(trisize / 3);

		mass = 0.001;

		for (uint32_t i = 0; i < vertsize; i++) {
			PositionList[i]	    = vertdata[i];
			RestPositionList[i] = vertdata[i];
			VelocityList[i]	    = fvec2(0.0);
		}

		for (uint32_t i = 0; i < trisize / 3; i++) {
			fvec2 X0 = RestPositionList[tridata[3 * i + 0]];
			fvec2 X1 = RestPositionList[tridata[3 * i + 1]];
			fvec2 X2 = RestPositionList[tridata[3 * i + 2]];

			AList[i] = fmat2(X1 - X0, X2 - X0).inverse();
			VList[i] = (X1 - X0).cross(X2 - X0) / 2.0;
			oList[i] = 0.0;
		}

		tva.resetvertarray(vertsize, trisize, tridata);
		tva.setcolor(0.8, 0.8, 0.8, 1.0);
		tva.settype(1);

		lva.resetvertarray(vertsize, edgesize, edgedata);
		lva.setcolor(0.0, 0.0, 0.0, 1.0);

		this->UpdataVertarray();
		this->Setdata();
	}

	void
	UpdataVertarray()
	{

#pragma omp parallel
		for (uint32_t i = 0; i < vertsize; i++) {
			const fvec2& v	   = PositionList[i];
			tva[i].position[0] = v.x;
			tva[i].position[1] = v.y;
			tva[i].position[2] = 0.0;
			lva[i].position[0] = v.x;
			lva[i].position[1] = v.y;
			lva[i].position[2] = -0.1;
		}
	}

	void
	Setdata()
	{
		tva.vboupdate();
		lva.vboupdate();
	}

	void ClearLamda()
	{
#pragma omp parallel
		for (auto& x : ElasticLamdalist) {
			x = 0.0;
		}
	}

	void FemElasticProjectGS(std::vector<fvec2>& TempPosition)
	{
		for (uint32_t i = 0; i < trisize / 3; i++) {
			const fvec2& X0 = RestPositionList[tridata[3 * i + 0]];
			const fvec2& X1 = RestPositionList[tridata[3 * i + 1]];
			const fvec2& X2 = RestPositionList[tridata[3 * i + 2]];

			const fvec2& x0 = TempPosition[tridata[3 * i + 0]];
			const fvec2& x1 = TempPosition[tridata[3 * i + 1]];
			const fvec2& x2 = TempPosition[tridata[3 * i + 2]];

			const fmat2& A = AList[i];
			const fmat2 F  = mat2(x1 - x0, x2 - x0) * A;

			const fmat2 E = 0.5 * (F.transpose() * F - fmat2::identity());

			const float& V = VList[i];

			float W;
			fvec2 dx0, dx1, dx2;

			if (MaterialInd == 0) {
				FemElasticDxStVenant(F, E, A, V, lambda, mu, W, dx0, dx1, dx2);
			} else if (MaterialInd == 1) {
				FemElasticDxNeoHookean(F, E, A, V, lambda, mu, W, dx0, dx1, dx2);
			} else if (MaterialInd == 2) {
				FemElasticDxCoRotational(F, E, A, oList[i], V, lambda, mu, W, dx0, dx1, dx2);
			}

			if (W > 0.0001) {

				float C = std::sqrt(2.0 * W);

				//fmat2 BAt = B * A.transpose();

				fvec2 dC1 = (1.0 / C) * dx1;
				fvec2 dC2 = (1.0 / C) * dx2;
				fvec2 dC0 = (1.0 / C) * dx0;

				//std::cout << dC1 << std::endl;

				float dtdtdlambda = (-C - ElasticLamdalist[i]) / ((dC0.sqlength() + dC1.sqlength() + dC2.sqlength()) / mass + 1.0 / (dt * dt));

				TempPosition[tridata[3 * i + 0]] = TempPosition[tridata[3 * i + 0]] + dtdtdlambda * (1.0 / mass) * dC0;
				TempPosition[tridata[3 * i + 1]] = TempPosition[tridata[3 * i + 1]] + dtdtdlambda * (1.0 / mass) * dC1;
				TempPosition[tridata[3 * i + 2]] = TempPosition[tridata[3 * i + 2]] + dtdtdlambda * (1.0 / mass) * dC2;

				ElasticLamdalist[i] += dtdtdlambda / (dt * dt);
			}

			if (F.det() < 0.00001 && lambda > 0.00001 && mu > 0.00001) {
				float omega = ExtractRotation(F, 1, oList[i]);
				oList[i]    = omega;
				fmat2 R(omega);

				//fmat2 S	  = R.transpose() * F;
				//float trS = S.trace();

				float W = 0.5 * (F - R).sqlength();
				fmat2 B = (F - R);

				float C = std::sqrt(2.0 * W);

				fmat2 BAt = B * A.transpose();

				fvec2 dC1 = fvec2(BAt.m[0], BAt.m[2]);
				fvec2 dC2 = fvec2(BAt.m[1], BAt.m[3]);
				fvec2 dC0 = -(dC1 + dC2);

				float dtdtdlambda = -C / ((dC0.sqlength() + dC1.sqlength() + dC2.sqlength()));
				dtdtdlambda *= 0.35;

				TempPosition[tridata[3 * i + 0]] = TempPosition[tridata[3 * i + 0]] + dtdtdlambda * dC0;
				TempPosition[tridata[3 * i + 1]] = TempPosition[tridata[3 * i + 1]] + dtdtdlambda * dC1;
				TempPosition[tridata[3 * i + 2]] = TempPosition[tridata[3 * i + 2]] + dtdtdlambda * dC2;
			}
		}
	}

	void FemElasticProjectJC(std::vector<fvec2>& TempPosition)
	{

		std::vector<fvec2> dx(trisize);

#pragma omp parallel for
		for (uint32_t i = 0; i < trisize / 3; i++) {
			const fvec2& X0 = RestPositionList[tridata[3 * i + 0]];
			const fvec2& X1 = RestPositionList[tridata[3 * i + 1]];
			const fvec2& X2 = RestPositionList[tridata[3 * i + 2]];

			const fvec2& x0 = TempPosition[tridata[3 * i + 0]];
			const fvec2& x1 = TempPosition[tridata[3 * i + 1]];
			const fvec2& x2 = TempPosition[tridata[3 * i + 2]];

			const fmat2& A = AList[i];
			const fmat2 F  = mat2(x1 - x0, x2 - x0) * A;

			const fmat2 E = 0.5 * (F.transpose() * F - fmat2::identity());

			const float& V = VList[i];

			float W;
			fvec2 dx0, dx1, dx2;

			if (MaterialInd == 0) {
				FemElasticDxStVenant(F, E, A, V, lambda, mu, W, dx0, dx1, dx2);
			} else if (MaterialInd == 1) {
				FemElasticDxNeoHookean(F, E, A, V, lambda, mu, W, dx0, dx1, dx2);
			} else if (MaterialInd == 2) {
				FemElasticDxCoRotational(F, E, A, oList[i], V, lambda, mu, W, dx0, dx1, dx2);
			}

			if (W > 0.0001) {
				float C = std::sqrt(2.0 * W);

				fvec2 dC1 = (1.0 / C) * dx1;
				fvec2 dC2 = (1.0 / C) * dx2;
				fvec2 dC0 = (1.0 / C) * dx0;

				float dtdtdlambda = (-C - ElasticLamdalist[i]) / ((dC0.sqlength() + dC1.sqlength() + dC2.sqlength()) / mass + 1.0 / (dt * dt));
				dtdtdlambda *= 0.2;

				dx[3 * i + 0] = dtdtdlambda * (1.0 / mass) * dC0;
				dx[3 * i + 1] = dtdtdlambda * (1.0 / mass) * dC1;
				dx[3 * i + 2] = dtdtdlambda * (1.0 / mass) * dC2;

				ElasticLamdalist[i] += dtdtdlambda / (dt * dt);
			} else {
				dx[3 * i + 0] = fvec2(0.0);
				dx[3 * i + 1] = fvec2(0.0);
				dx[3 * i + 2] = fvec2(0.0);
			}

			if (F.det() < 0.00001 && lambda > 0.00001 && mu > 0.00001) {
				float omega = ExtractRotation(F, 1, oList[i]);
				oList[i]    = omega;
				fmat2 R(omega);

				//fmat2 S	  = R.transpose() * F;
				//float trS = S.trace();

				float W = 0.5 * (F - R).sqlength();
				fmat2 B = (F - R);

				float C = std::sqrt(2.0 * W);

				fmat2 BAt = B * A.transpose();

				fvec2 dC1 = fvec2(BAt.m[0], BAt.m[2]);
				fvec2 dC2 = fvec2(BAt.m[1], BAt.m[3]);
				fvec2 dC0 = -(dC1 + dC2);

				float dtdtdlambda = -C / ((dC0.sqlength() + dC1.sqlength() + dC2.sqlength()));
				dtdtdlambda *= 0.35;

				dx[3 * i + 0] = dtdtdlambda * dC0;
				dx[3 * i + 1] = dtdtdlambda * dC1;
				dx[3 * i + 2] = dtdtdlambda * dC2;
			}
		}

#pragma omp parallel for
		for (uint32_t i = 0; i < trisize; i++) {
			TempPosition[tridata[i]] = TempPosition[tridata[i]] + dx[i];
		}
	}

	void FixedProjection(std::vector<fvec2>& TempPosition)
	{

		//RectMeshに依存

		for (uint32_t i = 0; i <= N * (M - 1); i += N)
			TempPosition[i] = RestPositionList[i];

		if (rightwall) {
			for (uint32_t y = 0; y < M; y++) {
				TempPosition[N * y + N - 1].x = RestPositionList[N * y + N - 1].x + trans;
				TempPosition[N * y + N - 1].y = RestPositionList[N * y + N - 1].y;
			}
		}
	}
};

namespace Physics {

int32_t solver = 0;

void timestep(ClothMesh& CM)
{

	CM.ClearLamda();

	uint32_t NodeSize = CM.vertsize;
	std::vector<fvec2> tempp(NodeSize);

#pragma omp parallel for
	for (uint32_t i = 0; i < NodeSize; i++) {
		fvec2 velocity = CM.VelocityList[i] + dt * fvec2(0.0, -9.8);
		tempp[i]       = CM.PositionList[i] + dt * velocity;
	}

	if (solver == 0)
		for (uint32_t x = 0; x < 100; x++) {
			CM.FemElasticProjectGS(tempp);
			CM.FixedProjection(tempp);
		}
	else if (solver == 1)
		for (uint32_t x = 0; x < 200; x++) {
			CM.FemElasticProjectJC(tempp);
			CM.FixedProjection(tempp);
		}

#pragma omp parallel for
	for (uint32_t i = 0; i < NodeSize; i++) {
		CM.VelocityList[i] = (tempp[i] - CM.PositionList[i]) / dt;
		CM.VelocityList[i] = 0.99 * CM.VelocityList[i];
		CM.PositionList[i] = tempp[i];
	}
}
};

int main(int argc, char const* argv[])
{
	using namespace std;

	Visualizer::Init(1024, 1024);

	IMGUI_CHECKVERSION();
	ImGui::CreateContext();
	ImGuiIO& io    = ImGui::GetIO();
	io.IniFilename = NULL;
	(void)io;

	ImGui::StyleColorsDark();

	ImGui_ImplOpenGL3_Init("#version 460 core");
	ImGui_ImplGlfw_InitForOpenGL(Visualizer::GetWindowPtr(), true);

	Renderer2D::Init();

	omp_set_num_threads(10);

	//init

	ClothMesh CM0(40, 8, 50.0, 10.0, fvec2(0.0));

	std::vector<Renderer2D::drawobject> renderlist;

	renderlist.emplace_back(CM0.tva);
	renderlist.emplace_back(CM0.lva);

	//rendering loop

	double ctime = 0.0;
	double vtime = 0.0;

	float width  = 5.0;
	float height = 5.0;
	fvec2 center;

	while (!Visualizer::Is_Closed()) {
		Renderer2D::Clear();
		//imgui reset
		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();

		static bool is_stop = true;

		static float scale = 5.0f;

		static bool nextframe = false;

		static int32_t rendermethod = 0;

		//physics
		ctime = Visualizer::GetTime();
		if (!is_stop || nextframe) {
			Physics::timestep(CM0);
			vtime += dt;
			nextframe = false;

			CM0.UpdataVertarray();
			CM0.Setdata();
		}

		//imgui
		{

			ImGui::Begin("2D FEM Elasticity");

			ImGui::Text("FPS : %.1f", ImGui::GetIO().Framerate);

			ImGui::Checkbox("stop", &is_stop);
			if (ImGui::Button("nextframe")) {
				nextframe = true;
			}

			ImVec2 mousemove = ImGui::GetMouseDragDelta(1);

			center.x += scale * mousemove.x / (1024 * 5);
			center.y += scale * mousemove.y / (1024 * 5);
			ImGui::Text(" centerx = %.1f", center.x);
			ImGui::Text(" centery = %.1f", center.y);

			float dlength = io.MouseWheel;
			scale += dlength * 1.0f;
			if (scale < 1.0)
				scale = 1.0;
			if (scale > 100.0)
				scale = 100.0;

			width  = scale;
			height = scale;

			ImGui::Text("scale = %.1f", scale);

			if (ImGui::Combo("RenderMethod", &(rendermethod), "Triangle\0Vertex\0Surface\0\0")) {
				if (rendermethod == 0) {
					renderlist[0].renderswitch = true;
					renderlist[1].renderswitch = true;
				} else if (rendermethod == 1) {
					renderlist[0].renderswitch = false;
					renderlist[1].renderswitch = false;

				} else if (rendermethod == 2) {
					renderlist[0].renderswitch = false;
					renderlist[1].renderswitch = false;
				}
			}

			ImGui::Combo("Solver", &(Physics::solver), "Gauss-Seidel\0Jacobi\0\0");
			ImGui::Combo("Material", &(CM0.MaterialInd), "Saint Venant-Kirchhoff\0Neo-Hookean\0Co-Rotational\0\0");
			ImGui::SliderFloat("mu", &mu, 0.0f, 500.0f, "%.lf", ImGuiSliderFlags_Logarithmic);
			ImGui::SliderFloat("lambda", &lambda, 0.0f, 500.0f, "%.lf", ImGuiSliderFlags_Logarithmic);

			ImGui::Checkbox("rightwall", &rightwall);
			ImGui::SliderFloat("trans", &trans, -50.0f, 50.0f);

			ImGui::Text("realtime = %.1f", ctime);
			ImGui::Text("virtualtime = %.1f", vtime);

			if (ImGui::Button("reset")) {

				vtime = 0.0;

				for (uint32_t i = 0; i < CM0.vertsize; i++) {
					CM0.PositionList[i].x = CM0.RestPositionList[i].x;
					CM0.PositionList[i].y = CM0.RestPositionList[i].y;
					CM0.VelocityList[i]   = fvec2(0.0);
					for (auto& x : CM0.oList)
						x = 0.0;
				}
			}

			ImGui::End();
		}

		if (rendermethod == 1) {
			for (uint32_t i = 0; i < CM0.trisize / 3; i++) {
				fvec2 x0 = CM0.PositionList[CM0.tridata[3 * i + 0]];
				fvec2 x1 = CM0.PositionList[CM0.tridata[3 * i + 1]];
				fvec2 x2 = CM0.PositionList[CM0.tridata[3 * i + 2]];

				fvec2 cm = (x0 + x1 + x2) / 3.0;

				x0 = (x0 - cm) * 0.8 + cm;
				x1 = (x1 - cm) * 0.8 + cm;
				x2 = (x2 - cm) * 0.8 + cm;

				Renderer2D::DrawLine(x0, x1);
				Renderer2D::DrawLine(x1, x2);
				Renderer2D::DrawLine(x2, x0);
			}
		}
		if (rendermethod == 2) {
			for (uint32_t i = 0; i < CM0.boundarysize / 2; i++) {
				fvec2 x0 = CM0.PositionList[CM0.boundarydata[2 * i + 0]];
				fvec2 x1 = CM0.PositionList[CM0.boundarydata[2 * i + 1]];

				fvec2 cm = (x0 + x1) / 2.0;

				x0 = (x0 - cm) * 0.8 + cm;
				x1 = (x1 - cm) * 0.8 + cm;

				Renderer2D::DrawLine(x0, x1);
			}
		}

		//renderer set

		Renderer2D::setcenter(center);
		Renderer2D::setWH(width, height);
		Renderer2D::updateUniformobj();

		//rendering

		Renderer2D::Draw(renderlist);

		ImGui::Render();
		ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

		//swapbuff
		Visualizer::SwapBuffer();
		//wait event
		Visualizer::PollEvent();
	}

	return 0;
}
