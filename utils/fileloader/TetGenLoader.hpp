#pragma once

#include "utils/mathfunc/mathfunc.hpp"

bool LoadTetGentoTetrahedraMesh(
    const char* filename,
    fvec3** const vertdata,
    uint32_t& vertsize,
    uint32_t** const ilistdata,
    uint32_t& ilistsize,
    uint32_t** const facelistdata,
    uint32_t& facelistsize,
    uint32_t** const edgelistdata,
    uint32_t& edgelistsize,
    const fvec3& meshoffset,
    const float& meshscale);
