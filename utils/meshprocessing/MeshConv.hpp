#pragma once

#include <cstdint>

#include "utils/mathfunc/mathfunc.hpp"

#include "opengl/drawobject.hpp"

void ConvertPTMtoREM(
    const fvec3* const Pvdata,
    const uint32_t& Pvsize,
    const uint32_t* const Pidata,
    const uint32_t& pisize,
    linevertarray& varray);
