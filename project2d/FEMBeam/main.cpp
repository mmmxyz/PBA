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
	std::vector<uint32_t> TriIndList;
	std::vector<uint32_t> InnerEdgeIndList;

	linevertarray lva;
	uint32_t* lilist;
	uint32_t lisize;

	trianglevertarray tva;
	uint32_t* tilist;
	uint32_t tisize;

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

		PositionList.reserve(N * M);
		RestPositionList.reserve(N * M);
		VelocityList.reserve(N * M);
		TriIndList.reserve(3 * 2 * (N - 1) * (M - 1));

		AList.reserve(2 * (N - 1) * (M - 1));
		VList.reserve(2 * (N - 1) * (M - 1));

		ElasticLamdalist.resize(3 * 2 * (N - 1) * (M - 1));

		for (uint32_t y = 0; y < M; y++) {
			float vy = -0.5 * lengthy + lengthy * (y / float(M - 1));
			for (uint32_t x = 0; x < N; x++) {
				float vx = -0.5 * lengthx + lengthx * (x / float(N - 1));
				PositionList.emplace_back(fvec2(vx, vy) + bias);
				RestPositionList.emplace_back(fvec2(vx, vy) + bias);
				VelocityList.emplace_back(fvec2(0.0));
			}
		}

		mass = rho * lengthx * lengthy / (M * N);

		for (uint32_t y = 0; y < M - 1; y++) {
			for (uint32_t x = 0; x < N - 1; x++) {
				uint32_t Ind0 = N * y + x;
				uint32_t Ind1 = Ind0 + 1;
				uint32_t Ind2 = Ind0 + N;
				uint32_t Ind3 = Ind0 + N + 1;

				TriIndList.emplace_back(Ind0);
				TriIndList.emplace_back(Ind1);
				TriIndList.emplace_back(Ind2);

				TriIndList.emplace_back(Ind2);
				TriIndList.emplace_back(Ind1);
				TriIndList.emplace_back(Ind3);
			}
		}

		for (uint32_t i = 0; i < 2 * (M - 1) * (N - 1); i++) {
			fvec2 X0 = RestPositionList[TriIndList[3 * i + 0]];
			fvec2 X1 = RestPositionList[TriIndList[3 * i + 1]];
			fvec2 X2 = RestPositionList[TriIndList[3 * i + 2]];

			AList.emplace_back(fmat2(X1 - X0, X2 - X0).inverse());
			VList.emplace_back((X1 - X0).cross(X2 - X0) / 2.0);
			oList.emplace_back(0.0);
		}

		uint32_t tvsize = N * M;

		tisize = 3 * 2 * (N - 1) * (M - 1);
		tilist = new uint32_t[tisize];

		for (uint32_t y = 0; y < M - 1; y++) {
			for (uint32_t x = 0; x < N - 1; x++) {
				uint32_t Ind0			  = N * y + x;
				tilist[6 * ((N - 1) * y + x) + 0] = Ind0;
				tilist[6 * ((N - 1) * y + x) + 1] = Ind0 + 1;
				tilist[6 * ((N - 1) * y + x) + 2] = Ind0 + N;

				tilist[6 * ((N - 1) * y + x) + 3] = Ind0 + N;
				tilist[6 * ((N - 1) * y + x) + 4] = Ind0 + 1;
				tilist[6 * ((N - 1) * y + x) + 5] = Ind0 + N + 1;
			}
		}

		tva.resetvertarray(tvsize, tisize, tilist);

		for (uint32_t i = 0; i < M * N; i++) {
			tva[i].color[0] = 0.8;
			tva[i].color[1] = 0.8;
			tva[i].color[2] = 0.8;
			tva[i].color[3] = 1.0;
			tva[i].type	= 1;
		}

		uint32_t lvsize = N * M;

		lisize = (6 * (N - 1) + 2) * (M - 1) + 2 * (N - 1);
		lilist = new uint32_t[lisize];

		for (uint32_t y = 0; y < M - 1; y++) {
			for (uint32_t x = 0; x < N - 1; x++) {
				uint32_t Ind0				  = N * y + x;
				lilist[(6 * (N - 1) + 2) * y + 6 * x + 0] = Ind0;
				lilist[(6 * (N - 1) + 2) * y + 6 * x + 1] = Ind0 + 1;

				lilist[(6 * (N - 1) + 2) * y + 6 * x + 2] = Ind0 + 1;
				lilist[(6 * (N - 1) + 2) * y + 6 * x + 3] = Ind0 + N;

				lilist[(6 * (N - 1) + 2) * y + 6 * x + 4] = Ind0 + N;
				lilist[(6 * (N - 1) + 2) * y + 6 * x + 5] = Ind0;
			}
			lilist[(6 * (N - 1) + 2) * y + 6 * (N - 1) + 0] = N * y + N - 1;
			lilist[(6 * (N - 1) + 2) * y + 6 * (N - 1) + 1] = N * y + 2 * N - 1;
		}
		for (uint32_t x = 0; x < N - 1; x++) {
			lilist[(6 * (N - 1) + 2) * (M - 1) + 2 * x + 0] = (M - 1) * N + x;
			lilist[(6 * (N - 1) + 2) * (M - 1) + 2 * x + 1] = (M - 1) * N + x + 1;
		}

		lva.resetvertarray(lvsize, lisize, lilist);

		for (uint32_t i = 0; i < M * N; i++) {
			lva[i].color[0] = 0.0;
			lva[i].color[1] = 0.0;
			lva[i].color[2] = 0.0;
			lva[i].color[3] = 1.0;
			lva[i].type	= 0;
		}

		this->UpdataVertarray();

		this->Setdata();
	}

	void
	UpdataVertarray()
	{

#pragma omp parallel
		for (uint32_t y = 0; y < M; y++) {
			for (uint32_t x = 0; x < N; x++) {
				fvec2 v			   = PositionList[N * y + x];
				tva[N * y + x].position[0] = v.x;
				tva[N * y + x].position[1] = v.y;
				tva[N * y + x].position[2] = 0.0;
				lva[N * y + x].position[0] = v.x;
				lva[N * y + x].position[1] = v.y;
				lva[N * y + x].position[2] = -0.1;
			}
		}
	}

	void Setdata()
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
		for (uint32_t i = 0; i < 2 * (M - 1) * (N - 1); i++) {
			fvec2 X0 = RestPositionList[TriIndList[3 * i + 0]];
			fvec2 X1 = RestPositionList[TriIndList[3 * i + 1]];
			fvec2 X2 = RestPositionList[TriIndList[3 * i + 2]];

			fvec2 x0 = TempPosition[TriIndList[3 * i + 0]];
			fvec2 x1 = TempPosition[TriIndList[3 * i + 1]];
			fvec2 x2 = TempPosition[TriIndList[3 * i + 2]];

			fmat2 A = AList[i];
			fmat2 F = mat2(x1 - x0, x2 - x0) * A;

			fmat2 E = 0.5 * (F.transpose() * F - fmat2::identity());

			float V = VList[i];
			V *= 1.0f;

			float W;
			fmat2 B;
			if (MaterialInd == 0) {
				W = V * (mu * E.sqlength() + 0.5 * lambda * E.trace() * E.trace());
				B = V * (2 * mu * F * E + lambda * E.trace() * F);
			} else if (MaterialInd == 1) {
				float J	   = F.det();
				float logJ = std::log(J);
				if (J < 0.0)
					logJ = 0.0;
				W = V * (0.5 * mu * (F.sqlength() - 2) - mu * logJ + 0.5 * lambda * logJ * logJ);
				B = V * (mu * F - mu * (F.inverse()).transpose() + lambda * logJ * (F.inverse()).transpose());
			} else if (MaterialInd == 2) {
				float omega = ExtractRotation(F, 3, oList[i]);
				oList[i]    = omega;
				fmat2 R(omega);
				fmat2 S	  = R.transpose() * F;
				float trS = S.trace();
				W	  = V * (0.5 * mu * ((F - R).sqlength()) + 0.5 * lambda * (trS * trS - 4 * trS + 4));
				B	  = V * (mu * (F - R) + lambda * (trS - 2) * R);
			}

			//std::cout << W << std::endl;
			//std::cout << V << std::endl;

			if (W > 0.0) {
				float C = std::sqrt(2.0 * W);

				fmat2 BAt = B * A.transpose();

				fvec2 dC1 = (1.0 / C) * fvec2(BAt.m[0], BAt.m[2]);
				fvec2 dC2 = (1.0 / C) * fvec2(BAt.m[1], BAt.m[3]);
				fvec2 dC0 = -(dC1 + dC2);

				//std::cout << dC1 << std::endl;

				float dtdtdlambda = (-C - ElasticLamdalist[i]) / ((dC0.sqlength() + dC1.sqlength() + dC2.sqlength()) / mass + 1.0 / (dt * dt));

				TempPosition[TriIndList[3 * i + 0]] = TempPosition[TriIndList[3 * i + 0]] + dtdtdlambda * (1.0 / mass) * dC0;
				TempPosition[TriIndList[3 * i + 1]] = TempPosition[TriIndList[3 * i + 1]] + dtdtdlambda * (1.0 / mass) * dC1;
				TempPosition[TriIndList[3 * i + 2]] = TempPosition[TriIndList[3 * i + 2]] + dtdtdlambda * (1.0 / mass) * dC2;

				ElasticLamdalist[i] += dtdtdlambda / (dt * dt);
			}
		}
	}

	void FemElasticProjectJC(std::vector<fvec2>& TempPosition)
	{

		std::vector<fvec2> dx(3 * 2 * (M - 1) * (N - 1));

#pragma omp parallel for
		for (uint32_t i = 0; i < 2 * (M - 1) * (N - 1); i++) {
			fvec2 X0 = RestPositionList[TriIndList[3 * i + 0]];
			fvec2 X1 = RestPositionList[TriIndList[3 * i + 1]];
			fvec2 X2 = RestPositionList[TriIndList[3 * i + 2]];

			fvec2 x0 = TempPosition[TriIndList[3 * i + 0]];
			fvec2 x1 = TempPosition[TriIndList[3 * i + 1]];
			fvec2 x2 = TempPosition[TriIndList[3 * i + 2]];

			fmat2 A = AList[i];
			fmat2 F = mat2(x1 - x0, x2 - x0) * A;

			fmat2 E = 0.5 * (F.transpose() * F - fmat2::identity());

			float V = VList[i];
			V *= 1.0f;

			float W;
			fmat2 B;
			if (MaterialInd == 0) {
				W = V * (mu * E.sqlength() + 0.5 * lambda * E.trace() * E.trace());
				B = V * (2 * mu * F * E + lambda * E.trace() * F);
			} else if (MaterialInd == 1) {
				float J	   = F.det();
				float logJ = std::log(J);
				if (J < 0.0)
					logJ = 0.0;
				W = V * (0.5 * mu * (F.sqlength() - 2) - mu * logJ + 0.5 * lambda * logJ * logJ);
				B = V * (mu * F - mu * (F.inverse()).transpose() + lambda * logJ * (F.inverse()).transpose());
			} else if (MaterialInd == 2) {
				float omega = ExtractRotation(F, 3, oList[i]);
				oList[i]    = omega;
				fmat2 R(omega);
				fmat2 S	  = R.transpose() * F;
				float trS = S.trace();
				W	  = V * (0.5 * mu * ((F - R).sqlength()) + 0.5 * lambda * (trS * trS - 4 * trS + 4));
				B	  = V * (mu * (F - R) + lambda * (trS - 2) * R);
			}

			//std::cout << W << std::endl;
			//std::cout << V << std::endl;

			if (W > 0.0) {
				float C = std::sqrt(2.0 * W);

				fmat2 BAt = B * A.transpose();

				fvec2 dC1 = (1.0 / C) * fvec2(BAt.m[0], BAt.m[2]);
				fvec2 dC2 = (1.0 / C) * fvec2(BAt.m[1], BAt.m[3]);
				fvec2 dC0 = -(dC1 + dC2);

				//std::cout << dC1 << std::endl;

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
		}

#pragma omp parallel for
		for (uint32_t i = 0; i < 3 * 2 * (M - 1) * (N - 1); i++) {
			TempPosition[TriIndList[i]] = TempPosition[TriIndList[i]] + dx[i];
		}
	}

	void FixedProjection(std::vector<fvec2>& TempPosition)
	{

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

	uint32_t NodeSize = CM.N * CM.M;
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
		for (uint32_t x = 0; x < 100; x++) {
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

	ClothMesh CM0(50, 10, 50.0, 10.0, fvec2(0.0));

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

				for (uint32_t i = 0; i < (CM0.M * CM0.N); i++) {
					CM0.PositionList[i].x = CM0.RestPositionList[i].x;
					CM0.PositionList[i].y = CM0.RestPositionList[i].y;
					CM0.VelocityList[i]   = fvec2(0.0);
					for (auto& x : CM0.oList)
						x = 0.0;
				}
			}

			ImGui::End();
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
