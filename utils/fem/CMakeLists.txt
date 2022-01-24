cmake_minimum_required(VERSION 3.10)

project(femLib CXX)

set(CMAKE_VERBOSE_MAKEFILE 1)
set(CMAKE_CXX_STANDARD 17)

add_library(femLib STATIC
	fem.cpp
)

target_include_directories(femLib PRIVATE
	../..
)
