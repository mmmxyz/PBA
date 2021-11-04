#pragma once

#include "utils/mathfunc/mathfunc.hpp"

void MeshCM(
    fvec3& Cm,
    float& Mass,
    const fvec3* const MVdata,
    const uint32_t& MVsize,
    const uint32_t* const MIlist,
    const uint32_t& Misize,
    const float& rho = 1.0);

void MeshInertia(
    fmat3& Inertia,
    const fvec3& cm,
    const fvec3* const MVdata,
    const uint32_t& MVsize,
    const uint32_t* const MIlist,
    const uint32_t& Misize,
    const float& rho = 1.0);
