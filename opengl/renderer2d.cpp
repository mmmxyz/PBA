#include <glad/glad.h>
#include <KHR/khrplatform.h>
#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>

#include "utils/mathfunc/mathfunc.hpp"
#include "opengl/visualizer.hpp"
#include "opengl/shader.hpp"
#include "opengl/texture.hpp"
#include "opengl/renderer2d.hpp"

namespace Renderer2D {

uint32_t persLoc;

static fvec2 center;
static float width, height;

static shader rendershader;

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
	rendershader.setprogram("../../../opengl/shadercode/vert2d.c", "../../../opengl/shadercode/frag2d.c");

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
	glEnable(GL_DEPTH_TEST);

	//vsync on
	glfwSwapInterval(1);

	rendershader.useprogram();
	persLoc = glGetUniformLocation(rendershader.program, "pers");

	width  = 5.0;
	height = 5.0;
	center = fvec2(0.0, 0.0);

	updateUniformobj();
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

		d.Va.draw();
	}
}

void Clear()
{
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}

}
