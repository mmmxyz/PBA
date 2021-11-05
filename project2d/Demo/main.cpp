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
#include "opengl/renderer2d.hpp"

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

	Renderer2D::Init();

	uint32_t ilist[3] = { 0, 1, 2 };
	trianglevertarray tritva(3, nullptr, 3, ilist);
	tritva.setposition(0, 0.0, 0.0, 0.0);
	tritva.setposition(1, 1.0, 0.0, 0.0);
	tritva.setposition(2, 0.5, 0.5, 0.0);
	tritva.setcolor(0, 1.0, 0.0, 0.0, 1.0);
	tritva.setcolor(1, 0.0, 1.0, 0.0, 1.0);
	tritva.setcolor(2, 0.0, 0.0, 1.0, 1.0);

	//init

	std::vector<Renderer2D::drawobject> renderlist;

	renderlist.emplace_back(Renderer2D::drawobject { tritva });

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

		static float scale = 1.0f;

		static bool nextframe = false;

		//imgui
		{

			ImGui::Begin("Demo");

			ImGui::Text("FPS : %.1f", ImGui::GetIO().Framerate);

			ImGui::Checkbox("stop", &is_stop);
			if (ImGui::Button("nextframe")) {
				nextframe = true;
			}

			ImVec2 mousemove = ImGui::GetMouseDragDelta(1);

			center.x += scale * mousemove.x / (1024 * 10);
			center.y += scale * mousemove.y / (1024 * 10);
			ImGui::Text(" width = %.1f", center.x);
			ImGui::Text(" heighti = %.1f", center.y);

			float dlength = io.MouseWheel;
			scale += dlength * 1.0f;
			if (scale < 1.0)
				scale = 1.0;
			if (scale > 10.0)
				scale = 10.0;

			width  = scale;
			height = scale;

			ImGui::Text("scale = %.1f", scale);

			ImGui::Text("realtime = %.1f", ctime);
			ImGui::Text("virtualtime = %.1f", vtime);

			if (ImGui::Button("reset")) {

				vtime = 0.0;
			}

			ImGui::End();
		}

		ctime = Visualizer::GetTime();
		if (!is_stop || nextframe) {
			vtime += dt;

			nextframe = false;
		}

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
