#include <iostream>
#include <cmath>
#include <cstdint>
#include <vector>
#include <omp.h>

#include "opengl/visualizer.hpp"
#include "opengl/vertarray.hpp"
#include "opengl/drawobject.hpp"
#include "opengl/renderer3d.hpp"

#include "utils/mathfunc/mathfunc.hpp"

#include "utils/fileloader/OBJLoader.hpp"

#include "utils/meshprocessing/MeshConv.hpp"

#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>

constexpr uint32_t FPS = 60;
constexpr float dt     = 1.0 / FPS;
//glfwSwapInterval(1); なので60FPS以下になる。
//これはモニタやGPUに依るかもしれない。

float muD = 0.1f;
float muC = 0.5f;

int32_t optimstatus = 0;

int32_t numitr = 5;

int32_t renderstatus = 1;

fvec3 heatmap(float x, float max, float min)
{
	float midpoint = 0.5 * (max + min);
	if (max < x)
		return fvec3(1.0, 0.0, 0.0);
	else if (x < min)
		return fvec3(0.0, 0.0, 1.0);
	else if (x < midpoint) {
		float t = (x - min) / (midpoint - min);
		return fvec3(0.0, t, 1.0 - t);
	} else {
		float t = (x - midpoint) / (max - midpoint);
		return fvec3(t, 1.0 - t, 0.0);
	}
}

class Mesh {
    public:
	std::vector<fvec3> PositionList;
	std::vector<fvec3> RestPositionList;
	std::vector<fvec3> dx;

	std::vector<float> Hmap;

	uint32_t vertsize;
	uint32_t* tilist;
	uint32_t tisize; //三角形数x3

	uint32_t* edgedata; //辺x2
	uint32_t edgesize;

	//triangle
	uint32_t* VtoTlist;
	uint32_t* VtoTind;

	trianglevertarray tva;
	linevertarray lva;

	Mesh()
	    : tva()
	    , lva()
	{

		fvec3* vertdata;

		LoadOBJtoPhysicsTriangleMesh("../../../resource/Bunny.obj", &vertdata, vertsize, &tilist, tisize, fvec3(0.0, -10.0, 0.0), 8.0);
		ConvertPTMtoPEM(vertsize, tilist, tisize, &edgedata, edgesize);

		RestPositionList.resize(vertsize);
		PositionList.resize(vertsize);
		Hmap.resize(vertsize);
		dx.resize(vertsize);

		for (uint32_t i = 0; i < vertsize; i++) {
			RestPositionList[i] = vertdata[i];
			PositionList[i]	    = vertdata[i];
			Hmap[i]		    = 0.0;
		}

		ConvertEVtoVE(vertsize, tilist, tisize, &VtoTind, &VtoTlist);

		tva.resetvertarray(vertsize, tisize, tilist);
		tva.settype(1);
		tva.setcolor(1.0, 1.0, 1.0, 1.0);

		lva.resetvertarray(vertsize, edgesize, edgedata);
		lva.settype(0);
		lva.setcolor(0.2, 0.2, 0.2, 0.5);

		this->setHmap();
		this->update();

		delete[] vertdata;
	}

	void update()
	{

		std::vector<fvec3> NormalSet(vertsize);

#pragma omp parallel for
		for (auto& x : NormalSet)
			x = fvec3(0.0);

#pragma omp parallel for
		for (uint32_t i = 0; i < vertsize; i++) {
			for (uint32_t j = VtoTind[i]; j < VtoTind[i + 1]; j++) {
				const uint32_t Tind = VtoTlist[j] / 3;

				const fvec3& v0 = PositionList[tilist[3 * Tind + 0]];
				const fvec3& v1 = PositionList[tilist[3 * Tind + 1]];
				const fvec3& v2 = PositionList[tilist[3 * Tind + 2]];

				fvec3 normal = (v1 - v0).cross(v2 - v0);

				NormalSet[i] = NormalSet[i] + normal;
			}
		}

#pragma omp parallel for
		for (uint32_t i = 0; i < vertsize; i++) {
			NormalSet[i] = NormalSet[i].normalize();
		}

#pragma omp parallel for
		for (uint32_t i = 0; i < vertsize; i++) {
			tva[i].normal[0] = NormalSet[i].x;
			tva[i].normal[1] = NormalSet[i].y;
			tva[i].normal[2] = NormalSet[i].z;

			tva[i].position[0] = PositionList[i].x;
			tva[i].position[1] = PositionList[i].y;
			tva[i].position[2] = PositionList[i].z;

			if (renderstatus == 2) {
				fvec3 hm	= heatmap(Hmap[i], 1.5, -1.5);
				tva[i].color[0] = hm.x;
				tva[i].color[1] = hm.y;
				tva[i].color[2] = hm.z;
			} else {
				tva[i].color[0] = 1.0;
				tva[i].color[1] = 1.0;
				tva[i].color[2] = 1.0;
			}

			lva[i].position[0] = PositionList[i].x + 0.001 * NormalSet[i].x;
			lva[i].position[1] = PositionList[i].y + 0.001 * NormalSet[i].y;
			lva[i].position[2] = PositionList[i].z + 0.001 * NormalSet[i].z;
		}

		tva.vboupdate();
		lva.vboupdate();
	}

	void setHmap()
	{

		std::vector<fvec3> HNSet(vertsize);
		std::vector<fvec3> NormalSet(vertsize);
		std::vector<float> AreaSet(vertsize);

#pragma omp parallel for
		for (auto& x : NormalSet)
			x = fvec3(0.0);

#pragma omp parallel for
		for (auto& x : HNSet)
			x = fvec3(0.0);

#pragma omp parallel for
		for (auto& x : AreaSet)
			x = 0.0;

#pragma omp parallel for
		for (auto& x : Hmap)
			x = 0.0;

#pragma omp parallel for
		for (uint32_t i = 0; i < vertsize; i++) {
			for (uint32_t j = VtoTind[i]; j < VtoTind[i + 1]; j++) {
				uint32_t Tind = VtoTlist[j] / 3;

				const fvec3& Rvi = RestPositionList[tilist[3 * Tind + (VtoTlist[j] % 3 + 0) % 3]];
				const fvec3& Rvj = RestPositionList[tilist[3 * Tind + (VtoTlist[j] % 3 + 1) % 3]];
				const fvec3& Rvk = RestPositionList[tilist[3 * Tind + (VtoTlist[j] % 3 + 2) % 3]];
				const fvec3& vi	 = PositionList[tilist[3 * Tind + (VtoTlist[j] % 3 + 0) % 3]];
				const fvec3& vj	 = PositionList[tilist[3 * Tind + (VtoTlist[j] % 3 + 1) % 3]];
				const fvec3& vk	 = PositionList[tilist[3 * Tind + (VtoTlist[j] % 3 + 2) % 3]];

				float cotj = (Rvk - Rvj).dot(Rvi - Rvj) / ((Rvk - Rvj).cross(Rvi - Rvj)).length();
				float cotk = (Rvi - Rvk).dot(Rvj - Rvk) / ((Rvi - Rvk).cross(Rvj - Rvk)).length();

				fvec3 Hn = 0.5 * (-0.5 * (cotj + cotk) * PositionList[i] + 0.5 * cotj * vk + 0.5 * cotk * vj);
				HNSet[i] = HNSet[i] + Hn;
				AreaSet[i] += (0.5 * ((Rvj - Rvi).cross(Rvk - Rvi)).length()) / 3.0;

				fvec3 normal = (vj - vi).cross(vk - vi);
				NormalSet[i] = NormalSet[i] + normal;
			}
		}

#pragma omp parallel for
		for (uint32_t i = 0; i < vertsize; i++) {
			NormalSet[i] = NormalSet[i].normalize();
		}

#pragma omp parallel for
		for (uint32_t i = 0; i < vertsize; i++) {
			Hmap[i] += HNSet[i].dot(NormalSet[i]) / AreaSet[i];
		}
	}

	void projection()
	{

		for (auto& x : dx)
			x = fvec3(0.0);

#pragma omp parallel for
		for (uint32_t i = 0; i < vertsize; i++) {

			fvec3 HN = fvec3(0.0);

			for (uint32_t j = VtoTind[i]; j < VtoTind[i + 1]; j++) {
				uint32_t Tind = VtoTlist[j] / 3;

				const fvec3& Rvi = RestPositionList[tilist[3 * Tind + (VtoTlist[j] % 3 + 0) % 3]];
				const fvec3& Rvj = RestPositionList[tilist[3 * Tind + (VtoTlist[j] % 3 + 1) % 3]];
				const fvec3& Rvk = RestPositionList[tilist[3 * Tind + (VtoTlist[j] % 3 + 2) % 3]];
				const fvec3& vi	 = PositionList[tilist[3 * Tind + (VtoTlist[j] % 3 + 0) % 3]];
				const fvec3& vj	 = PositionList[tilist[3 * Tind + (VtoTlist[j] % 3 + 1) % 3]];
				const fvec3& vk	 = PositionList[tilist[3 * Tind + (VtoTlist[j] % 3 + 2) % 3]];

				if (optimstatus == 0) {
					if (((Rvj - Rvi).cross(Rvk - Rvi)).sqlength() / ((Rvj - Rvi).sqlength() * (Rvk - Rvi).sqlength()) > 0.01) {
						const float cotj = (Rvk - Rvj).dot(Rvi - Rvj) / ((Rvk - Rvj).cross(Rvi - Rvj)).length();
						const float cotk = (Rvi - Rvk).dot(Rvj - Rvk) / ((Rvi - Rvk).cross(Rvj - Rvk)).length();
						fvec3 dC	 = 0.5 * (cotj + cotk) * PositionList[i] - 0.5 * cotj * vk - 0.5 * cotk * vj;
						HN		 = HN + muC * dC;
					}
				} else if (optimstatus == 1) {
					fvec3 dC = 2.0 * PositionList[i] - vk - vj;
					HN	 = HN + muD * dC;
				}
			}

			dx[i] = -HN;
		}

#pragma omp parallel for
		for (uint32_t i = 0; i < vertsize; i++) {
			PositionList[i] = PositionList[i] + (1.0 / numitr) * dx[i];
		}
	}

	void smooth()
	{
		for (uint32_t i = 0; i < numitr; i++) {
			projection();
		}

		this->setHmap();
		this->update();
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

	Renderer3D::Init();

	omp_set_num_threads(20);

	uint32_t cagelist[24] = { 0, 1, 1, 2, 2, 3, 3, 0, 4, 5, 5, 6,
		6, 7, 7, 4, 0, 4, 3, 7, 1, 5, 2, 6 };
	linevertarray cage(8, 24, cagelist);
	cage.setposition(0, -15.0f, -15.0f, -15.0f);
	cage.setposition(1, 15.0f, -15.0f, -15.0f);
	cage.setposition(2, 15.0f, 15.0f, -15.0f);
	cage.setposition(3, -15.0f, 15.0f, -15.0f);
	cage.setposition(4, -15.0f, -15.0f, 15.0f);
	cage.setposition(5, 15.0f, -15.0f, 15.0f);
	cage.setposition(6, 15.0f, 15.0f, 15.0f);
	cage.setposition(7, -15.0f, 15.0f, 15.0f);
	cage.settype(0);
	cage.vboupdate();

	uint32_t floorlist[6] = { 0, 1, 2, 0, 2, 3 };
	trianglevertarray floor(6, 6, floorlist);
	floor.setposition(0, -12.0f, -12.0f, -12.0f);
	floor.setposition(1, -12.0f, -12.0f, 12.0f);
	floor.setposition(2, 12.0f, -12.0f, 12.0f);
	floor.setposition(3, 12.0f, -12.0f, -12.0f);
	for (uint32_t i = 0; i < 4; i++) {
		floor.setnormal(i, 0.0f, 1.0f, 0.0f);
	}
	floor.setuv(0, 0.0f, 1.0f);
	floor.setuv(1, 0.0f, 0.0f);
	floor.setuv(2, 1.0f, 0.0f);
	floor.setuv(3, 1.0f, 1.0f);
	floor.settype(2);
	floor.vboupdate();

	uint8_t* images0 = new uint8_t[1024 * 4 * 1024];
	for (uint32_t i = 0; i < 1024; i++) {
		for (uint32_t j = 0; j < 1024; j++) {
			if ((((i / 32) % 2) + ((j / 32) % 2)) == 1) {
				images0[1024 * 4 * i + 4 * j + 0] = 255;
				images0[1024 * 4 * i + 4 * j + 1] = 255;
				images0[1024 * 4 * i + 4 * j + 2] = 255;
				images0[1024 * 4 * i + 4 * j + 3] = 255;
			} else {
				images0[1024 * 4 * i + 4 * j + 0] = 50;
				images0[1024 * 4 * i + 4 * j + 1] = 50;
				images0[1024 * 4 * i + 4 * j + 2] = 50;
				images0[1024 * 4 * i + 4 * j + 3] = 255;
			}
		}
	}
	texture simasima0(1024, 1024, images0);

	fvec3 camerap(0.0, 15.0, 20.0);

	Mesh M;

	//init

	std::vector<Renderer3D::drawobject> shadowlist;
	std::vector<Renderer3D::drawobject> edgelist;
	std::vector<Renderer3D::drawobject> renderlist;

	shadowlist.emplace_back(floor);
	shadowlist.emplace_back(M.tva);

	renderlist.emplace_back(floor, simasima0);
	renderlist.emplace_back(cage);
	renderlist.emplace_back(M.tva);
	renderlist.emplace_back(M.lva);

	//random

	Renderer3D::setLightint(600.0);
	Renderer3D::setclookat(fvec3(0.0, -7.0, 0.0));

	//rendering loop

	double ctime = 0.0;
	double vtime = 0.0;

	double prevtime = 0.0;

	while (!Visualizer::Is_Closed()) {
		Renderer3D::Clear();
		//imgui reset
		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();
		ImGui::SetNextWindowSize(ImVec2(300, 360), ImGuiCond_FirstUseEver);

		static float mousemovex = 0.5 * 3.14;
		static float mousemovey = 0.5 * 3.14;
		static float cameraL	= 30.0f;

		float camerasinp;
		float cameracosp;
		float camerasint;
		float cameracost;

		static float Lposix = 5.0;
		static float Lposiy = 15.0;
		static float Lposiz = 18.0;

		//imgui
		{

			ImGui::Begin("Mesh Smoothing");

			ImGui::Text("Rendering FPS : %.1f", ImGui::GetIO().Framerate);

			if (ImGui::Combo("render status", &renderstatus, "Surface\0Surface Edge\0Mean Curvature\0\0")) {
				if (renderstatus == 0) {
					M.tva.settype(1);
					renderlist[3].renderswitch = false;
				} else if (renderstatus == 1) {
					renderlist[3].renderswitch = true;
				} else if (renderstatus == 2) {
					M.tva.settype(0);
					renderlist[3].renderswitch = true;
				}
				M.update();
			}

			ImVec2 mousemove = ImGui::GetMouseDragDelta(1);

			mousemovex += mousemove.x / (1024 * 5);
			mousemovey += mousemove.y / (1024 * 5);
			ImGui::Text(" theta = %.1f", mousemovey);
			ImGui::Text(" phi = %.1f", mousemovex);

			camerasinp = std::sin(mousemovex);
			cameracosp = std::cos(mousemovex);
			camerasint = std::sin(mousemovey);
			cameracost = std::cos(mousemovey);

			float dlength = io.MouseWheel;
			cameraL += dlength * 1.0f;
			if (cameraL < 10.0)
				cameraL = 10.0;
			if (cameraL > 80.0)
				cameraL = 80.0;

			camerap.x = cameraL * camerasint * cameracosp;
			camerap.z = cameraL * camerasint * camerasinp;
			camerap.y = cameraL * cameracost;

			ImGui::Text(" x = %.1f", camerap.x);
			ImGui::Text(" y = %.1f", camerap.y);
			ImGui::Text(" z = %.1f", camerap.z);

			ImGui::Text("camera length = %.1f", cameraL);

			ImGui::SliderFloat("Lightx", &Lposix, -20.0f, 20.0f);
			ImGui::SliderFloat("Lighty", &Lposiy, 0.0f, 20.0f);
			ImGui::SliderFloat("Lightz", &Lposiz, -20.0f, 20.0f);

			ImGui::Combo("laplacian", &optimstatus, "Cotangent\0Discrete\0\0");
			if (optimstatus == 0)
				ImGui::SliderFloat("Cotangent mu", &muC, 0.0f, 2.0f);
			else if (optimstatus == 1)
				ImGui::SliderFloat("Discrete mu", &muD, 0.0f, 2.0f);

			ImGui::SliderInt("num iteration", &numitr, 1, 200);

			if (ImGui::Button("smoothing")) {
				M.smooth();
			}
			if (ImGui::Button("reset")) {
				for (uint32_t i = 0; i < M.vertsize; i++) {
					M.PositionList[i] = M.RestPositionList[i];
				}
				M.setHmap();
				M.update();
			}

			ImGui::End();
		}

		//renderer parameters
		Renderer3D::setLposi(fvec3(Lposix, Lposiy, Lposiz));
		Renderer3D::DrawPoint(fvec3(Lposix, Lposiy, Lposiz), 1.0, 1.0, 1.0);
		Renderer3D::setcposi(camerap);
		Renderer3D::updateUniformobj();

		//rendering
		Renderer3D::Draw(shadowlist, edgelist, renderlist);

		ImGui::Render();
		ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

		//swapbuff
		Visualizer::SwapBuffer();
		//wait event
		Visualizer::PollEvent();
	}

	return 0;
}
