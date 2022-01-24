cmake_minimum_required(VERSION 3.10)

project(geoLib CXX)

set(CMAKE_VERBOSE_MAKEFILE 1)
set(CMAKE_CXX_STANDARD 17)

add_library(geoLib STATIC
	curve.cpp
)

target_include_directories(geoLib PRIVATE
	../..
)
