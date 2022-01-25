#pragma once

#include <vector>

#include "utils/mathfunc/mathfunc.hpp"

#include "opengl/visualizer.hpp"
#include "opengl/vertarray.hpp"
#include "opengl/shader.hpp"
#include "opengl/texture.hpp"

namespace Renderer2D {

struct drawobject {
	const vertarray& Va;
	bool renderswitch;
	drawobject(const vertarray& va)
	    : Va(va)
	    , renderswitch(true)
	{
	}
};

bool Init();

void Draw(const std::vector<drawobject>& renderlist);

void setcenter(const fvec2& c);
void setWH(const float W, const float H);

void updateUniformobj();

void Clear();

void DrawLine(const fvec2& x0, const fvec2& x1, const float& r, const float& g, const float& b);
void DrawLine(const fvec2& x0, const fvec2& x1, const fvec3& color);
void DrawLine(const fvec2& x0, const fvec2& x1);

void DrawPoint(const fvec2& x, const float& r, const float& g, const float& b);
void DrawPoint(const fvec2& x, const fvec3& color);
void DrawPoint(const fvec2& x);

void DrawPolyLine(const fvec2* const X, const uint32_t size, const float& r, const float& g, const float& b);
void DrawPolyLine(const fvec2* const X, const uint32_t size, const fvec3& color);
void DrawPolyLine(const fvec2* const X, const uint32_t size);

}
