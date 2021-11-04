#pragma once

#include <cstdint>

#include "utils/mathfunc/mathfunc.hpp"

#include "opengl/vertarray.hpp"

void ConvertPTMtoREM(
    const fvec3* const Pvdata,
    const uint32_t& Pvsize,
    const uint32_t* const Pidata,
    const uint32_t& pisize,
    vertex** const Rvedata,
    uint32_t& Rvsize,
    uint32_t** const Ridata,
    uint32_t& Risize);
