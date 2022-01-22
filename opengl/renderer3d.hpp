#pragma once

#include <vector>

#include "utils/mathfunc/mathfunc.hpp"

#include "opengl/visualizer.hpp"
#include "opengl/vertarray.hpp"
#include "opengl/shader.hpp"
#include "opengl/texture.hpp"

namespace Renderer3D {

struct drawobject {
	const vertarray& Va;
	const texture* const pTex;
	const fquaternion* const pAmat;
	const fvec3* const pAvec;

	drawobject(const vertarray& va, const texture& Tex, const fquaternion& Amat, const fvec3& Avec)
	    : Va(va)
	    , pTex(&Tex)
	    , pAmat(&Amat)
	    , pAvec(&Avec)
	{
	}
	drawobject(const vertarray& va, const fquaternion& Amat, const fvec3& Avec)
	    : Va(va)
	    , pTex(nullptr)
	    , pAmat(&Amat)
	    , pAvec(&Avec)
	{
	}
	drawobject(const vertarray& va, const texture& Tex)
	    : Va(va)
	    , pTex(&Tex)
	    , pAmat(nullptr)
	    , pAvec(nullptr)
	{
	}
	drawobject(const vertarray& va)
	    : Va(va)
	    , pTex(nullptr)
	    , pAmat(nullptr)
	    , pAvec(nullptr)
	{
	}
};

bool Init();

void Draw(
    const std::vector<drawobject>& shadowlist,
    const std::vector<drawobject>& edgelist,
    const std::vector<drawobject>& renderlist);

void setLposi(const fvec3& x);
void setLlookat(const fvec3& x);
void setcposi(const fvec3& x);
void setclookat(const fvec3& x);
void setAffine(const fmat3& mat, const fvec3 vec);
void setAffine(void);
void setLightint(const float& x);

//void useshadowshader();
//void userenderershader();
//void useedgeshader();

void updateUniformobj();

void Clear();

//単純な描画関数，OpenGL側のバッファを使い回す
//Renderer3Dでは描画の情報は保持しないので，一度きりの描画になる．
//必要ならmain.cpp側で描画する場所などを記録するべし
//シェーディングはできない
void DrawLine(const fvec3& x0, const fvec3& x1, const float& r, const float& g, const float& b);
void DrawLine(const fvec3& x0, const fvec3& x1, const fvec3& color);
void DrawLine(const fvec3& x0, const fvec3& x1);

void DrawPoint(const fvec3& x, const float& r, const float& g, const float& b);
void DrawPoint(const fvec3& x, const fvec3& color);
void DrawPoint(const fvec3& x);

void DrawPolyLine(const fvec3* const X, const uint32_t size, const float& r, const float& g, const float& b);
void DrawPolyLine(const fvec3* const X, const uint32_t size, const fvec3& color);
void DrawPolyLine(const fvec3* const X, const uint32_t size);

void DrawTriangle(const fvec3& x0, const fvec3& x1, const fvec3& x2, const float& r, const float& g, const float& b);
void DrawTriangle(const fvec3& x0, const fvec3& x1, const fvec3& x2, const fvec3& color);
void DrawTriangle(const fvec3& x0, const fvec3& x1, const fvec3& x2);

void DrawTetrahedron(const fvec3& x0, const fvec3& x1, const fvec3& x2, const fvec3& x3, const float& r, const float& g, const float& b);
void DrawTetrahedron(const fvec3& x0, const fvec3& x1, const fvec3& x2, const fvec3& x3, const fvec3& color);
void DrawTetrahedron(const fvec3& x0, const fvec3& x1, const fvec3& x2, const fvec3& x3);

}
