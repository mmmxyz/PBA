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
};

bool Init();

void Draw(const std::vector<drawobject>& renderlist);

void setcenter(const fvec2& c);
void setWH(const float W, const float H);

void updateUniformobj();

void Clear();

}
