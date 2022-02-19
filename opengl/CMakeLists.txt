cmake_minimum_required(VERSION 3.10)

project(glLib CXX)

set(CMAKE_VERBOSE_MAKEFILE 1)
set(CMAKE_CXX_STANDARD 17)

if(POLICY CMP0072)
    cmake_policy(SET CMP0072 NEW)
endif()
find_package(OpenGL REQUIRED)
find_package(glfw3 REQUIRED)

set(EX_INC "")
if(APPLE)
	set(EX_INC "/usr/local/include")
endif()

#message(STATUS "CMAKE_CURRENT_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}")

add_library(glLib STATIC
	drawobject.cpp
	vertarray.cpp
	visualizer.cpp
	renderer3d.cpp
	renderer2d.cpp
	texture.cpp
	shader.cpp
)

target_compile_definitions(glLib PRIVATE ABSOLUTEPATH=\"${CMAKE_CURRENT_SOURCE_DIR}\")

target_include_directories(glLib PRIVATE
	../3rdparty/glad/include
	..
	${OPENGL_INCLUDE_DIRS}
	${glfw3_INCLUDE_DIRS}
	${EX_INC}
)
