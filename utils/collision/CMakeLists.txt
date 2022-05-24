cmake_minimum_required(VERSION 3.10)

project(collisionLib CXX)

set(CMAKE_VERBOSE_MAKEFILE 1)
set(CMAKE_CXX_STANDARD 17)

add_library(collisionLib STATIC
	RigidObbBvh.cpp
	#primitive.cpp
	primitive2d.cpp
	#DeformableBvh.cpp
	DeformableBvh2D.cpp
	DeformableBvh3D_Tetra.cpp
	DeformableBvh3D_Tri.cpp
	basic3d.cpp
)

target_include_directories(collisionLib PRIVATE
	../..
)
