cmake_minimum_required(VERSION 3.10)

project(meshgeneratorLib CXX)

set(CMAKE_VERBOSE_MAKEFILE 1)
set(CMAKE_CXX_STANDARD 17)

add_library(meshgeneratorLib STATIC
	meshgenerator.cpp
)

include_directories(meshgeneratorLib PRIVATE
	../..
)
