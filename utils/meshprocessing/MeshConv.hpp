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

void ConvertEVtoVE(
    const uint32_t vertsize,
    const uint32_t* const elementlist,
    const uint32_t elementsize,
    uint32_t** const elsup_index,
    uint32_t** const elsup);
