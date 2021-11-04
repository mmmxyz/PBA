cmake_minimum_required(VERSION 3.10)

project(RigidCollision CXX)

set(CMAKE_VERBOSE_MAKEFILE 1)
set(CMAKE_CXX_STANDARD 17)


if(POLICY CMP0072)
    cmake_policy(SET CMP0072 NEW)
endif()
find_package(OpenGL REQUIRED)
find_package(glfw3 REQUIRED)


add_subdirectory(../../opengl ../../../opengl/build)
add_subdirectory(../../utils/mathfunc ../../../utils/mathfunc/build)
add_subdirectory(../../utils/fileloader ../../../utils/fileloader/build)
add_subdirectory(../../utils/collision ../../../utils/collision/build)
add_subdirectory(../../utils/meshprocessing ../../../utils/meshprocessing/build)
add_subdirectory(../../3rdparty/glad ../../../3rdparty/glad/build)
add_subdirectory(../../3rdparty/imgui ../../../3rdparty/imgui/build)

add_executable(Program main.cpp)

target_include_directories(Program PRIVATE
	../..
)

if (APPLE)
	link_directories(Program
		"/usr/local/lib"
		)
endif()

set(EX_GLFWLIB "")
if (APPLE)
	set(EX_GLFWLIB "-framework Cocoa -framework OpenGL -framework IOKit")
endif()

target_link_libraries(Program PRIVATE
	glLib
	imguiLib
	collisionLib
	fileloaderLib
	meshprosLib
	mathLib
	${EX_GLFWLIB}
	${OPENGL_LIBRARY}
	glfw
	gladLib
	${CMAKE_DL_LIBS}
)
