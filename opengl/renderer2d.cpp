#include <glad/glad.h>
#include <KHR/khrplatform.h>
#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>

#include "utils/mathfunc/mathfunc.hpp"
#include "opengl/visualizer.hpp"
#include "opengl/shader.hpp"
#include "opengl/texture.hpp"
#include "opengl/renderer2d.hpp"
#include "opengl/drawobject.hpp"

namespace Renderer2D {

uint32_t persLoc;

static fvec2 center;
static float width, height;

static shader rendershader;

static linevertarray Line;
static pointvertarray Point;
static linestripvertarray PolyLine;

static trianglevertarray Triangle;

fmat4 makeperspective()
{
	fmat4 pmat;
	pmat.m[0]  = 1.0 / width;
	pmat.m[5]  = 1.0 / height;
	pmat.m[10] = 1.0;
	pmat.m[15] = 1.0;

	pmat.m[3] = -center.x / width;
	pmat.m[7] = -center.y / height;

	return pmat;
}

bool Init()
{
	rendershader.setprogram("vert2d.c", "frag2d.c");

	//mics OpenGL setting for 3D rendering
	glPointSize(10.0f);

	//clear color set
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	//clear depth set (same as default value)
	glClearDepth(1.0);

	//alpha blend on
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_BLEND);

	//depth test
	glDisable(GL_DEPTH_TEST);

	//vsync on
	glfwSwapInterval(1);

	rendershader.useprogram();
	persLoc = glGetUniformLocation(rendershader.program, "pers");

	width  = 5.0;
	height = 5.0;
	center = fvec2(0.0, 0.0);

	updateUniformobj();

	Line.resetvertarray(2);
	Point.resetvertarray(1);
	PolyLine.resetvertarray(1024);
	Triangle.resetvertarray(3);

	return true;
}

void setcenter(const fvec2& c)
{
	center = c;
}

void setWH(const float W, const float H)
{
	width  = W;
	height = H;
}

void updateUniformobj()
{
	glUniformMatrix4fv(persLoc, 1, false, makeperspective().transpose().m);
}

void Draw(const std::vector<drawobject>& renderlist)
{
	updateUniformobj();

	for (const auto& d : renderlist) {
		if (d.renderswitch)
			d.Va.draw();
	}
}

void Clear()
{
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}

void DrawLine(const fvec2& x0, const fvec2& x1, const float& r, const float& g, const float& b)
{
	Line.setposition(0, x0);
	Line.setposition(1, x1);
	Line.setcolor(r, g, b, 1.0);
	Line.vboupdate();
	Line.draw();
}

void DrawLine(const fvec2& x0, const fvec2& x1, const fvec3& color)
{
	DrawLine(x0, x1, color.x, color.y, color.z);
}

void DrawLine(const fvec2& x0, const fvec2& x1)
{
	DrawLine(x0, x1, 0.0, 0.0, 0.0);
}

void DrawPoint(const fvec2& x, const float& r, const float& g, const float& b)
{
	Point.setposition(0, x);
	Point.setcolor(r, g, b, 1.0);
	Point.vboupdate();
	Point.draw();
}

void DrawPoint(const fvec2& x, const fvec3& color)
{
	DrawPoint(x, color.x, color.y, color.z);
}

void DrawPoint(const fvec2& x)
{
	DrawPoint(x, 0.0, 0.0, 0.0);
}

void DrawPolyLine(const fvec2* const X, const uint32_t size, const float& r, const float& g, const float& b)
{

	uint32_t cnt = 0;
	for (; size - cnt > 1024; cnt += 1024) {
		PolyLine.setallposition(X + cnt);
		PolyLine.setcolor(r, g, b, 1.0);
		PolyLine.draw();
	}
	PolyLine.setallposition(X + cnt, 0, size - cnt);
	PolyLine.setallposition(X[size - 1], size - cnt, 1024);
	PolyLine.setcolor(r, g, b, 1.0);
	PolyLine.vboupdate();
	PolyLine.draw();
}

void DrawPolyLine(const fvec2* const X, const uint32_t size, const fvec3& color)
{
	DrawPolyLine(X, size, color.x, color.y, color.z);
}

void DrawPolyLine(const fvec2* const X, const uint32_t size)
{
	DrawPolyLine(X, size, 0.0, 0.0, 0.0);
}

void DrawTriangle(const fvec2& x0, const fvec2& x1, const fvec2& x2, const float& r, const float& g, const float& b)
{

	Triangle.setposition(0, x0);
	Triangle.setposition(1, x1);
	Triangle.setposition(2, x2);

	Triangle.setcolor(r, g, b, 1.0);
	Triangle.vboupdate();
	Triangle.draw();
}

void DrawTriangle(const fvec2& x0, const fvec2& x1, const fvec2& x2, const fvec3& color)
{

	DrawTriangle(x0, x1, x2, color.x, color.y, color.z);
}

void DrawTriangle(const fvec2& x0, const fvec2& x1, const fvec2& x2)
{
	DrawTriangle(x0, x1, x2, 0.0, 0.0, 0.0);
}

}
