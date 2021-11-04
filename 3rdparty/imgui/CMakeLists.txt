cmake_minimum_required(VERSION 3.10)

project(imguiLib)

add_library(imguiLib 
	imgui.cpp
	imgui_draw.cpp
	imgui_tables.cpp
	imgui_widgets.cpp
	backend/imgui_impl_glfw.cpp
	backend/imgui_impl_opengl3.cpp
)

target_include_directories(imguiLib PUBLIC
	.
	./backend
)


