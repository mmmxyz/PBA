#pragma once

struct GLFWwindow;

namespace Visualizer {

bool Init(const uint32_t& w = 1024, const uint32_t& h = 1024);

bool Is_Closed();

void SwapBuffer();

void PollEvent();

void WaitEvent(const double& sec);

double GetTime();

void setViewport(const uint32_t w, const uint32_t h);
void setViewport();

GLFWwindow* GetWindowPtr();

}
