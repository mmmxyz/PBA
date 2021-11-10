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
	texture* const pTex;
	fquaternion* const pAmat;
	fvec3* const pAvec;
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

}
