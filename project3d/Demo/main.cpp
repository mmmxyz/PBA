#include <iostream>
#include <cmath>
#include <cstdint>
#include <vector>
#include <deque>
#include <algorithm>
#include <random>

#include "opengl/visualizer.hpp"
#include "opengl/vertarray.hpp"
#include "opengl/drawobject.hpp"
#include "opengl/renderer3d.hpp"

#include "utils/mathfunc/mathfunc.hpp"
#include "utils/fileloader/OBJLoader.hpp"

#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>

constexpr uint32_t FPS = 60;
constexpr float dt     = 1.0 / FPS;

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

	uint32_t cagelist[24] = { 0, 1, 1, 2, 2, 3, 3, 0, 4, 5, 5, 6,
		6, 7, 7, 4, 0, 4, 3, 7, 1, 5, 2, 6 };
	linevertarray cage(8, nullptr, 24, cagelist);
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
	trianglevertarray floor(6, nullptr, 6, floorlist);
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

	trianglevertarray Bunnytva;

	LoadOBJtoRenderTriangleMesh("../../../resource/Cube.obj", Bunnytva, fvec3(0.0, -4.0, 0.0), 5.00);

	Bunnytva.settype(2);
	Bunnytva.vboupdate();

	fvec3 cm = fvec3(0.0);
	fquaternion rotq(fvec3(0.0, 2.0 * 3.1415, 0.0));

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

	uint8_t* images1 = new uint8_t[1024 * 4 * 1024];
	for (uint32_t i = 0; i < 1024; i++) {
		for (uint32_t j = 0; j < 1024; j++) {
			if ((((i / 64) % 2) + ((j / 64) % 2)) == 1) {
				images1[1024 * 4 * i + 4 * j + 0] = 255;
				images1[1024 * 4 * i + 4 * j + 1] = 255;
				images1[1024 * 4 * i + 4 * j + 2] = 255;
				images1[1024 * 4 * i + 4 * j + 3] = 255;
			} else {
				images1[1024 * 4 * i + 4 * j + 0] = 50;
				images1[1024 * 4 * i + 4 * j + 1] = 50;
				images1[1024 * 4 * i + 4 * j + 2] = 50;
				images1[1024 * 4 * i + 4 * j + 3] = 255;
			}
		}
	}
	texture simasima1(1024, 1024, images1);

	fvec3 camerap(0.0, 15.0, 20.0);

	//init

	std::vector<Renderer3D::drawobject> shadowlist;
	std::vector<Renderer3D::drawobject> edgelist;
	std::vector<Renderer3D::drawobject> renderlist;

	shadowlist.emplace_back(floor);
	shadowlist.emplace_back(Bunnytva, rotq, cm);

	renderlist.emplace_back(floor, simasima0);
	renderlist.emplace_back(Bunnytva, simasima1, rotq, cm);
	renderlist.emplace_back(cage);

	//rendering loop

	double ctime = 0.0;
	double vtime = 0.0;

	float cmx = 0.0;

	while (!Visualizer::Is_Closed()) {
		Renderer3D::Clear();
		//imgui reset
		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();

		static bool is_stop = true;

		static float mousemovex = 0.5 * 3.14;
		static float mousemovey = 0.5 * 3.14;
		static float cameraL	= 30.0f;

		float camerasinp;
		float cameracosp;
		float camerasint;
		float cameracost;

		static bool nextframe = false;

		//physics
		ctime = Visualizer::GetTime();

		if (!is_stop || nextframe) {
			vtime += dt;

			cmx += 0.2;
			if (cmx > 5.0)
				cmx = -5.0;
			cm.x = cmx;

			rotq = fquaternion(fvec3(0.1, 0.01 * cmx, 0.1)) * rotq;

			Renderer3D::DrawLine(fvec3(0.0), fvec3(0.1 * cmx, 10.0, 0.1 * cmx));
			Renderer3D::DrawPoint(fvec3(0.1 * cmx, 10.0, 0.1 * cmx));

			constexpr uint32_t hogesize = 256;
			fvec3 hogef[hogesize];
			for (uint32_t i = 0; i < hogesize; i++)
				hogef[i] = fvec3(i * 0.001 - 9.0, 10 * sin(i * 0.005 * cmx), 0.0);
			hogef[hogesize - 1] = fvec3(10.0, 10.0, 5.0);

			Renderer3D::DrawPolyLine(hogef, hogesize);

			nextframe = false;
		}

		//imgui
		{

			ImGui::Begin("Demo");

			ImGui::Text("FPS : %.1f", ImGui::GetIO().Framerate);

			ImGui::Checkbox("stop", &is_stop);
			if (ImGui::Button("nextframe")) {
				nextframe = true;
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

			ImGui::Text("realtime = %.1f", ctime);
			ImGui::Text("virtualtime = %.1f", vtime);

			if (ImGui::Button("reset")) {

				vtime = 0.0;
			}

			ImGui::End();
		}

		//renderer set

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
