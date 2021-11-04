#include <glad/glad.h>
#include <KHR/khrplatform.h>
#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>

#include <iostream>
#include <cstdint>

#include "opengl/visualizer.hpp"

namespace Visualizer {

static GLFWwindow* window;
static uint32_t width, height;

bool Init(const uint32_t& w, const uint32_t& h)
{

	width  = w;
	height = h;

	if (glfwInit() == GL_FALSE) {
		std::cout << "fail to initialize glfw" << std::endl;
		return false;
	}

	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 6);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

	atexit(glfwTerminate);

	window = glfwCreateWindow(width, height, "Physically Based Animation", NULL, NULL);

	if (window == NULL) {
		std::cerr << "fail to create window" << std::endl;
		return false;
	}

	glfwMakeContextCurrent(window);

	if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
		std::cout << "fail to create OpenGL context" << std::endl;
		return false;
	}

	//Version
	std::cout << "OpenGL Version : " << glGetString(GL_VERSION) << std::endl;
	std::cout << "GLSL Version : " << glGetString(GL_SHADING_LANGUAGE_VERSION) << std::endl;
	std::cout << "Vender : " << glGetString(GL_VENDOR) << std::endl;
	std::cout << "Renderer : " << glGetString(GL_RENDERER) << std::endl;
	std::cout << std::endl;

	return true;
}

bool Is_Closed()
{
	return glfwWindowShouldClose(window);
}

void SwapBuffer()
{
	glfwSwapBuffers(window);
}

void PollEvent()
{
	glfwPollEvents();
}

void WaitEvent(const double& sec)
{
	glfwWaitEventsTimeout(sec);
}

double GetTime()
{
	return glfwGetTime();
}

void setViewport(const uint32_t w, const uint32_t h)
{
	glViewport(0, 0, w, h);
}

void setViewport()
{
	glViewport(0, 0, width, height);
}

GLFWwindow* GetWindowPtr()
{
	return window;
}
}
