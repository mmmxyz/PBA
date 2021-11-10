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

find_package(glfw3 REQUIRED)

set(EX_INC "")
if(APPLE)
	set(EX_INC "/usr/local/include")
endif()

target_include_directories(imguiLib PUBLIC
	.
	./backend
)

target_include_directories(imguiLib PRIVATE
	${glfw3_INCLUDE_DIRS}
	${EX_INC}
)
