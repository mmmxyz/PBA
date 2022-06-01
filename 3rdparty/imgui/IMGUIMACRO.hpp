#pragma once

#define InitializeIMGUI                              \
	IMGUI_CHECKVERSION();                        \
	ImGui::CreateContext();                      \
	ImGuiIO& io    = ImGui::GetIO();             \
	io.IniFilename = NULL;                       \
	(void)io;                                    \
	ImGui::StyleColorsDark();                    \
	ImGui_ImplOpenGL3_Init("#version 460 core"); \
	ImGui_ImplGlfw_InitForOpenGL(Visualizer::GetWindowPtr(), true);

#define NewframeIMGUI(x, y)           \
	ImGui_ImplOpenGL3_NewFrame(); \
	ImGui_ImplGlfw_NewFrame();    \
	ImGui::NewFrame();            \
	ImGui::SetNextWindowSize(ImVec2((x), (y)), ImGuiCond_FirstUseEver);
