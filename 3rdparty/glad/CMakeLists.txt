cmake_minimum_required(VERSION 3.10)

project(gladLib)

add_library(gladLib 
	src/glad.c
)

target_include_directories(gladLib PRIVATE
	include/
)

