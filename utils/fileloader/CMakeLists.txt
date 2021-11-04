cmake_minimum_required(VERSION 3.10)

project(fileloaderLib CXX)

set(CMAKE_VERBOSE_MAKEFILE 1)
set(CMAKE_CXX_STANDARD 17)

include_directories(
	../..
)

add_library(fileloaderLib STATIC
	OBJLoader.cpp
)
