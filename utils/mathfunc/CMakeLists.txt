cmake_minimum_required(VERSION 3.10)

project(mathLib CXX)

set(CMAKE_VERBOSE_MAKEFILE 1)
set(CMAKE_CXX_STANDARD 17)

add_library(mathLib STATIC
	mathfunc.cpp
	polardecompose.cpp
)

target_include_directories(mathLib PRIVATE
	../..
)
