#include <glad/glad.h>
#include <KHR/khrplatform.h>
#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>

#include "utils/mathfunc/mathfunc.hpp"
#include "opengl/visualizer.hpp"
#include "opengl/shader.hpp"
#include "opengl/texture.hpp"
#include "opengl/renderer3d.hpp"
#include "opengl/drawobject.hpp"

namespace Renderer3D {

static fvec3 Lposi, Llookat, Lup, cameraposi, cameralookat, cameraup;
static float Lnear, Lfar, Lright, Lleft, Ltop, Lbottom;
static float near, far, right, left, top, bottom;
static float Lint, Lambient;

static fmat3 Affinemat;
static fvec3 Affinevec;

static uint32_t shadowpersLoc, shadoweuclidLoc;
static uint32_t edgepersLoc, edgeeuclidLoc;
static uint32_t LpersLoc, LeuclidLoc, persLoc, euclidLoc;
static uint32_t LightLoc, LintLoc, LambLoc;
static uint32_t extraeuclidLoc;

static uint32_t shadowtex;
static uint32_t fbo;

static shader shadowshader, rendershader, edgeshader;
static uint32_t currentshadernum;
// 0:shadow
// 1:renderer
// 2:edge

static uint32_t counter;

static uint32_t shadowwidth, shadowheight;

//単純な描画関数のためのdrawobject
static linevertarray Line;
static pointvertarray Point;
static linestripvertarray PolyLine;
static trianglevertarray Triangle;

static linevertarray TetraLine;
static trianglevertarray TetraTri;
//初期のときはOpenGLの関数が使えないので，それらを呼び出さないデフォルトコンストラクタを呼ぶ

fmat4 makeperspective(const float& near, const float& far, const float& right, const float& left, const float& top,
    const float& bottom)
{
	fmat4 pmat;
	pmat.m[0]  = (2 * near) / (right - left);
	pmat.m[2]  = (right + left) / (right - left);
	pmat.m[5]  = (2 * near) / (top - bottom);
	pmat.m[6]  = (top + bottom) / (top - bottom);
	pmat.m[10] = -(far + near) / (far - near);
	pmat.m[11] = -(2 * near * far) / (far - near);
	pmat.m[14] = -1.0f;
	return pmat;
}

fmat4 makeeuclid(const fvec3& cameraposition, const fvec3& cameralookat, const fvec3& cameraup)
{
	fvec3 rz = (cameraposition - cameralookat).normalize();
	fvec3 rx = (cameraup.cross(rz)).normalize();
	fvec3 ry = (rz.cross(rx)).normalize();

	fmat3 rot(rx, ry, rz);
	rot	= rot.transpose();
	fvec3 t = rot * (-cameraposition);

	return fmat4(rot, t);
}

void setLposi(const fvec3& x)
{
	Lposi = x;
}
void setLlookat(const fvec3& x)
{
	Llookat = x;
}

void setcposi(const fvec3& x)
{
	cameraposi = x;
}
void setclookat(const fvec3& x)
{
	cameralookat = x;
}

void setAffine(const fmat3& mat, const fvec3 vec)
{
	Affinemat = mat;
	Affinevec = vec;
}

void setAffine(void)
{
	Affinemat = fmat3::indentity();
	Affinevec = fvec3(0.0, 0.0, 0.0);
}

void setLightint(const float& x)
{
	Lint = x;
}

bool Init()
{

	shadowshader.setprogram("../../../opengl/shadercode/shadowvert.c", "../../../opengl/shadercode/shadowfrag.c");
	edgeshader.setprogram("../../../opengl/shadercode/edgevert.c", "../../../opengl/shadercode/edgefrag.c");
	rendershader.setprogram("../../../opengl/shadercode/vert.c", "../../../opengl/shadercode/frag.c");

	counter	     = 0;
	shadowwidth  = 4096;
	shadowheight = 4096;

	//mics OpenGL setting for 3D rendering
	glPointSize(10.0f);

	//clear color set
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	//clear depth set (same as default value)
	glClearDepth(1.0);

	//alpha blend on
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_BLEND);

	//face culling
	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);

	//depth test
	glEnable(GL_DEPTH_TEST);

	//vsync on
	glfwSwapInterval(1);

	//get Uniform Object location

	shadowshader.useprogram();
	shadowpersLoc	= glGetUniformLocation(shadowshader.program, "pers");
	shadoweuclidLoc = glGetUniformLocation(shadowshader.program, "euclid");

	edgeshader.useprogram();
	edgepersLoc   = glGetUniformLocation(edgeshader.program, "pers");
	edgeeuclidLoc = glGetUniformLocation(edgeshader.program, "euclid");

	rendershader.useprogram();
	LpersLoc       = glGetUniformLocation(rendershader.program, "Lpers");
	LeuclidLoc     = glGetUniformLocation(rendershader.program, "Leuclid");
	persLoc	       = glGetUniformLocation(rendershader.program, "pers");
	euclidLoc      = glGetUniformLocation(rendershader.program, "euclid");
	extraeuclidLoc = glGetUniformLocation(rendershader.program, "extraeuclid");

	LightLoc = glGetUniformLocation(rendershader.program, "Lpos");
	LintLoc	 = glGetUniformLocation(rendershader.program, "Lint");
	LambLoc	 = glGetUniformLocation(rendershader.program, "Lamb");

	GLuint texLoc = glGetUniformLocation(rendershader.program, "colortexture");
	glUniform1i(texLoc, 2);

	//set initial uniform value
	cameraposi   = fvec3(0.0, 5.0, 20.0);
	cameralookat = fvec3(0.0, -10.0, 0.0);
	cameraup     = fvec3(0.0, 1.0, 0.0);

	near   = 1.0f;
	far    = 100.0f;
	right  = 0.8f;
	left   = -0.8f;
	top    = 0.8f;
	bottom = -0.8f;

	Lposi	= fvec3(0.0f, 15.0f, 18.0f);
	Llookat = fvec3(0.0f, -12.0f, -6.0f);
	Lup	= fvec3(0.0f, 0.0f, -1.0f);

	Lnear	= 1.0f;
	Lfar	= 50.0f;
	Lright	= 1.0f;
	Lleft	= -1.0f;
	Ltop	= 1.0f;
	Lbottom = -1.0f;

	Lint	 = 800.0f;
	Lambient = 0.3f;

	Affinemat = fmat3::indentity();
	Affinevec = fvec3(0.0, 0.0, 0.0);

	//create shadow texture
	glGenTextures(1, &shadowtex);
	glBindTexture(GL_TEXTURE_2D, shadowtex);

	glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT, shadowwidth, shadowheight, 0, GL_DEPTH_COMPONENT, GL_FLOAT, NULL);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_FUNC, GL_LEQUAL);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_MODE, GL_COMPARE_REF_TO_TEXTURE);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);

	float bordercolor[4] = { -1.0, 0.0, 0.0, 0.0 };
	glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, bordercolor);

	glBindTexture(GL_TEXTURE_2D, 0);

	//frame buffer
	glGenFramebuffers(1, &fbo);
	glBindFramebuffer(GL_FRAMEBUFFER, fbo);
	glBindTexture(GL_TEXTURE_2D, shadowtex);
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, shadowtex, 0);
	glDrawBuffer(GL_NONE);
	glReadBuffer(GL_NONE);

	if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
		exit(-1);
	}
	glBindFramebuffer(GL_FRAMEBUFFER, 0);

	//texture mapper
	GLuint shadowtexLocation(glGetUniformLocation(rendershader.program, "shadowtex"));
	glUniform1i(shadowtexLocation, 1);

	glActiveTexture(GL_TEXTURE0 + 1);
	glBindTexture(GL_TEXTURE_2D, shadowtex);
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, 0);

	updateUniformobj();

	//setup drawobject
	Line.resetvertarray(2);
	Point.resetvertarray(1);
	PolyLine.resetvertarray(1024);
	Triangle.resetvertarray(3);

	uint32_t hoge[12];
	hoge[0]	 = 0;
	hoge[1]	 = 1;
	hoge[2]	 = 0;
	hoge[3]	 = 2;
	hoge[4]	 = 0;
	hoge[5]	 = 3;
	hoge[6]	 = 1;
	hoge[7]	 = 2;
	hoge[8]	 = 1;
	hoge[9]	 = 3;
	hoge[10] = 2;
	hoge[11] = 3;
	TetraLine.resetvertarray(4, 12, hoge);

	hoge[0] = 0;
	hoge[1] = 2;
	hoge[2] = 1;

	hoge[3] = 1;
	hoge[4] = 2;
	hoge[5] = 3;

	hoge[6] = 0;
	hoge[7] = 1;
	hoge[8] = 3;

	hoge[9]	 = 0;
	hoge[10] = 3;
	hoge[11] = 2;
	TetraTri.resetvertarray(4, 12, hoge);

	return true;
}

void useshadowshader()
{
	glBindFramebuffer(GL_FRAMEBUFFER, fbo);
	shadowshader.useprogram();
	Visualizer::setViewport(shadowwidth, shadowheight);
	currentshadernum = 0;
}

void userenderershader()
{
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
	rendershader.useprogram();
	glCullFace(GL_BACK);
	Visualizer::setViewport();
	currentshadernum = 1;
}

void useedgeshader()
{
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
	edgeshader.useprogram();
	glCullFace(GL_FRONT);
	Visualizer::setViewport();
	currentshadernum = 2;
}

void updateUniformobj()
{
	shadowshader.useprogram();
	glUniformMatrix4fv(shadowpersLoc, 1, false,
	    makeperspective(Lnear, Lfar, Lright, Lleft, Ltop, Lbottom).transpose().m);
	glUniformMatrix4fv(shadoweuclidLoc, 1, false, (makeeuclid(Lposi, Llookat, Lup) * fmat4(Affinemat, Affinevec)).transpose().m);

	rendershader.useprogram();
	glUniformMatrix4fv(LpersLoc, 1, false,
	    makeperspective(Lnear, Lfar, Lright, Lleft, Ltop, Lbottom).transpose().m);
	glUniformMatrix4fv(LeuclidLoc, 1, false, makeeuclid(Lposi, Llookat, Lup).transpose().m);
	glUniformMatrix4fv(persLoc, 1, false, makeperspective(near, far, right, left, top, bottom).transpose().m);
	glUniformMatrix4fv(euclidLoc, 1, false, makeeuclid(cameraposi, cameralookat, cameraup).transpose().m);
	glUniformMatrix4fv(extraeuclidLoc, 1, false, fmat4(Affinemat, Affinevec).transpose().m);
	glUniform3f(LightLoc, Lposi.x, Lposi.y, Lposi.z);
	glUniform1f(LintLoc, Lint);
	glUniform1f(LambLoc, Lambient);

	edgeshader.useprogram();
	glUniformMatrix4fv(edgepersLoc, 1, false, makeperspective(near, far, right, left, top, bottom).transpose().m);
	glUniformMatrix4fv(edgeeuclidLoc, 1, false, (makeeuclid(cameraposi, cameralookat, cameraup) * fmat4(Affinemat, Affinevec)).transpose().m);

	if (currentshadernum == 0)
		shadowshader.useprogram();
	if (currentshadernum == 1)
		rendershader.useprogram();
	if (currentshadernum == 2)
		edgeshader.useprogram();
}

void updateRigidTransform()
{
	shadowshader.useprogram();
	glUniformMatrix4fv(shadoweuclidLoc, 1, false, (makeeuclid(Lposi, Llookat, Lup) * fmat4(Affinemat, Affinevec)).transpose().m);

	rendershader.useprogram();
	glUniformMatrix4fv(extraeuclidLoc, 1, false, fmat4(Affinemat, Affinevec).transpose().m);

	edgeshader.useprogram();
	glUniformMatrix4fv(edgeeuclidLoc, 1, false, (makeeuclid(cameraposi, cameralookat, cameraup) * fmat4(Affinemat, Affinevec)).transpose().m);

	if (currentshadernum == 0)
		shadowshader.useprogram();
	if (currentshadernum == 1)
		rendershader.useprogram();
	if (currentshadernum == 2)
		edgeshader.useprogram();
}

void drawshadow(
    const std::vector<drawobject>& shadowlist)
{
	useshadowshader();

	for (const auto& d : shadowlist) {
		if (d.pAmat != nullptr && d.pAvec != nullptr) {
			setAffine((*d.pAmat).qtomat(), (*d.pAvec));
		} else
			setAffine();

		updateRigidTransform();

		d.Va.draw();
	}
}

void drawedge(
    const std::vector<drawobject>& edgelist)

{
	useedgeshader();

	counter++;
	if (counter >= 7) {
		for (const auto& d : edgelist) {
			if (d.pAmat != nullptr && d.pAvec != nullptr) {
				setAffine((*d.pAmat).qtomat(), (*d.pAvec));
			} else
				setAffine();

			updateRigidTransform();

			d.Va.draw();
		}
	}
	if (counter > 10)
		counter = 0;
}

void drawrender(
    const std::vector<drawobject>& renderlist)
{
	userenderershader();
	for (const auto& d : renderlist) {
		if (d.pTex != nullptr) {
			d.pTex->activate(2);
		} else
			texture::deactivate(2);

		if (d.pAmat != nullptr && d.pAvec != nullptr) {
			setAffine((*d.pAmat).qtomat(), (*d.pAvec));
		} else
			setAffine();

		updateRigidTransform();

		d.Va.draw();
	}
}

void Draw(
    const std::vector<drawobject>& shadowlist,
    const std::vector<drawobject>& edgelist,
    const std::vector<drawobject>& renderlist)
{
	drawshadow(shadowlist);
	drawedge(edgelist);
	drawrender(renderlist);
}

void Clear()
{
	//clear buffer
	glBindFramebuffer(GL_FRAMEBUFFER, fbo);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	if (currentshadernum == 0)
		glBindFramebuffer(GL_FRAMEBUFFER, fbo);
}

void DrawLine(const fvec3& x0, const fvec3& x1, const float& r, const float& g, const float& b)
{
	Line.setposition(0, x0);
	Line.setposition(1, x1);
	Line.setcolor(r, g, b, 1.0);
	Line.vboupdate();
	Line.draw();
}

void DrawLine(const fvec3& x0, const fvec3& x1, const fvec3& color)
{
	DrawLine(x0, x1, color.x, color.y, color.z);
}

void DrawLine(const fvec3& x0, const fvec3& x1)
{
	DrawLine(x0, x1, 1.0, 1.0, 1.0);
}

void DrawPoint(const fvec3& x, const float& r, const float& g, const float& b)
{
	Point.setposition(0, x);
	Point.setcolor(r, g, b, 1.0);
	Point.vboupdate();
	Point.draw();
}

void DrawPoint(const fvec3& x, const fvec3& color)
{
	DrawPoint(x, color.x, color.y, color.z);
}

void DrawPoint(const fvec3& x)
{
	DrawPoint(x, 1.0, 1.0, 1.0);
}

void DrawPolyLine(const fvec3* const X, const uint32_t size, const float& r, const float& g, const float& b)
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

void DrawPolyLine(const fvec3* const X, const uint32_t size, const fvec3& color)
{
	DrawPolyLine(X, size, color.x, color.y, color.z);
}

void DrawPolyLine(const fvec3* const X, const uint32_t size)
{
	DrawPolyLine(X, size, 1.0, 1.0, 1.0);
}

void DrawTriangle(const fvec3& x0, const fvec3& x1, const fvec3& x2, const float& r, const float& g, const float& b)
{

	Triangle.setposition(0, x0);
	Triangle.setposition(1, x1);
	Triangle.setposition(2, x2);
	Triangle.setcolor(r, g, b, 1.0);
	Triangle.vboupdate();
	Triangle.draw();
}

void DrawTriangle(const fvec3& x0, const fvec3& x1, const fvec3& x2, const fvec3& color)
{
	DrawTriangle(x0, x1, x2, color.x, color.y, color.z);
}
void DrawTriangle(const fvec3& x0, const fvec3& x1, const fvec3& x2)
{

	DrawTriangle(x0, x1, x2, 1.0, 1.0, 1.0);
}

void DrawTetrahedron(const fvec3& x0, const fvec3& x1, const fvec3& x2, const fvec3& x3, const float& r, const float& g, const float& b)
{
	TetraTri.setposition(0, x0);
	TetraTri.setposition(1, x1);
	TetraTri.setposition(2, x2);
	TetraTri.setposition(3, x3);
	TetraTri.setcolor(r, g, b, 1.0);
	TetraTri.vboupdate();
	TetraTri.draw();

	TetraLine.setposition(0, x0);
	TetraLine.setposition(1, x1);
	TetraLine.setposition(2, x2);
	TetraLine.setposition(3, x3);
	TetraLine.setcolor(r * 0.3, g * 0.3, b * 0.3, 1.0);
	TetraLine.vboupdate();
	TetraLine.draw();
}

void DrawTetrahedron(const fvec3& x0, const fvec3& x1, const fvec3& x2, const fvec3& x3, const fvec3& color)
{
	DrawTetrahedron(x0, x1, x2, x3, color.x, color.y, color.z);
}

void DrawTetrahedron(const fvec3& x0, const fvec3& x1, const fvec3& x2, const fvec3& x3)
{
	DrawTetrahedron(x0, x1, x2, x3, 1.0, 1.0, 1.0);
}

}
