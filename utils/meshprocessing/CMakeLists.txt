cmake_minimum_required(VERSION 3.10)

project(meshprosLib CXX)

set(CMAKE_VERBOSE_MAKEFILE 1)
set(CMAKE_CXX_STANDARD 17)

add_library(meshprosLib STATIC
	IntOnMesh.cpp
	MeshConv.cpp
)

target_include_directories(meshprosLib PRIVATE
	../..
)
