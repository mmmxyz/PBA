#include <cstdint>
#include <iostream>
#include <cmath>
#include <omp.h>

#include "opengl/visualizer.hpp"
#include "opengl/vertarray.hpp"
#include "opengl/drawobject.hpp"
#include "opengl/renderer3d.hpp"

#include "utils/mathfunc/mathfunc.hpp"
#include "utils/geometry/curve.hpp"

#include "utils/mathfunc/mathfunc.hpp"

#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>

constexpr uint32_t FPS = 60;
constexpr float dt     = 1.0 / FPS;
//glfwSwapInterval(1); なので60FPS以下になる。
//これはモニタやGPUに依るかもしれない。

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

	fvec3 camerap(0.0, 0.0, 20.0);
	fvec3 cameralookat(0.0, 0.0, 0.0);

	//init

	std::vector<Renderer3D::drawobject> shadowlist;
	std::vector<Renderer3D::drawobject> edgelist;
	std::vector<Renderer3D::drawobject> renderlist;

	//init curve

	const uint32_t N = 20;
	fvec3 vecs[N];
	for (uint32_t i = 0; i < N; i++) {
		vecs[i] = fvec3(-8.0 + 20 * i / float(N), 10 * exp(-0.001 * float(i)) * sin(float(i)), 10.0 * sin(4.0 * float(i)));
	}

	bezier_curve bc(N, vecs);
	bezier_curve diffbc1 = bc.diff();
	bezier_curve diffbc2 = diffbc1.diff();

	b_spline bs(10, N - 1, vecs, nullptr);
	b_spline diffbs1 = bs.diff();
	b_spline diffbs2 = diffbs1.diff();

	//piecewise_cubic_hermite pch(N - 1, vecs, derib);
	catmull_rom catmull(N - 1, vecs);

	//rendering loop

	while (!Visualizer::Is_Closed()) {
		Renderer3D::Clear();
		//imgui reset
		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();
		ImGui::SetNextWindowSize(ImVec2(300, 250), ImGuiCond_FirstUseEver);

		static float mousemovex = 0.5 * 3.14;
		static float mousemovey = 0.5 * 3.14;
		static float cameraL	= 30.0f;

		static float vx = 0.0;

		static int32_t curvetype = 0;
		//0 bezier
		//1 b-spline
		//2 catmull

		static int32_t bsplinedegree = 10;

		float camerasinp;
		float cameracosp;
		float camerasint;
		float cameracost;

		//imgui
		{

			ImGui::Begin("Curves");

			ImGui::Text("Rendering FPS : %.1f", ImGui::GetIO().Framerate);

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

			ImGui::Combo("Curve", &(curvetype), "bezier\0b-spline\0catmull\0\0");

			ImGui::SliderFloat("t", &vx, 0.0f, 1.0f);

			if (curvetype == 0) {
			} else if (curvetype == 1) {
				int32_t hoge = bsplinedegree;
				ImGui::InputInt("degree", &hoge);
				if (hoge < 3)
					hoge = 3;
				if (hoge >= N)
					hoge = N - 1;
				bsplinedegree = hoge;
				bs.setdegree(bsplinedegree, nullptr);
				diffbs1 = bs.diff();
				diffbs2 = diffbs1.diff();
			} else if (curvetype == 2) {
			}

			ImGui::End();
		}

		uint32_t psize = 1000;
		fvec3 plotline[psize];

		if (curvetype == 0) {
			for (uint32_t i = 0; i < psize; i++)
				plotline[i] = bc(i / float(psize - 1));
			Renderer3D::DrawPoint(bc(0.0));
			Renderer3D::DrawPoint(bc(vx));
			Renderer3D::DrawPoint(bc(1.0));

			fvec3 e1 = diffbc1(vx).normalize();
			fvec3 e2 = (diffbc1(vx).cross(diffbc2(vx))).cross(diffbc1(vx)).normalize();
			fvec3 e3 = (diffbc1(vx).cross(diffbc2(vx))).normalize();

			Renderer3D::DrawLine(bc(vx), 3.0 * e1 + bc(vx), 1.0, 0.2, 0.2);
			Renderer3D::DrawLine(bc(vx), 3.0 * e2 + bc(vx), 0.2, 1.0, 0.2);
			Renderer3D::DrawLine(bc(vx), 3.0 * e3 + bc(vx), 0.2, 0.2, 1.0);

		} else if (curvetype == 1) {
			for (uint32_t i = 0; i < psize; i++)
				plotline[i] = bs(i / float(psize - 1));
			Renderer3D::DrawPoint(bs(0.0));
			Renderer3D::DrawPoint(bs(vx));
			Renderer3D::DrawPoint(bs(1.0));

			fvec3 e1 = diffbs1(vx).normalize();
			fvec3 e2 = (diffbs1(vx).cross(diffbs2(vx))).cross(diffbs1(vx)).normalize();
			fvec3 e3 = (diffbs1(vx).cross(diffbs2(vx))).normalize();

			Renderer3D::DrawLine(bs(vx), 3.0 * e1 + bs(vx), 1.0, 0.2, 0.2);
			Renderer3D::DrawLine(bs(vx), 3.0 * e2 + bs(vx), 0.2, 1.0, 0.2);
			Renderer3D::DrawLine(bs(vx), 3.0 * e3 + bs(vx), 0.2, 0.2, 1.0);
		} else if (curvetype == 2) {
			for (uint32_t i = 0; i < psize; i++)
				plotline[i] = catmull(i / float(psize - 1));
			Renderer3D::DrawPoint(catmull(0.0));
			Renderer3D::DrawPoint(catmull(vx));
			Renderer3D::DrawPoint(catmull(1.0));
		}

		Renderer3D::DrawPolyLine(plotline, psize);

		for (uint32_t i = 0; i < N; i++)
			Renderer3D::DrawPoint(vecs[i], 1.0, 0.0, 0.0);
		Renderer3D::DrawPolyLine(vecs, N, 0.0, 1.0, 0.0);

		//renderer set

		Renderer3D::setcposi(camerap);
		Renderer3D::setclookat(cameralookat);
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
