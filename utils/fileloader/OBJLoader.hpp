#pragma once

#include "opengl/vertarray.hpp"

bool LoadOBJtoRenderTriangleMesh(
    const char* filename,
    vertex** const vertdata,
    uint32_t& vertsize,
    uint32_t** const ilistdata,
    uint32_t& ilistsize,
    const fvec3& meshoffset = fvec3(0.0),
    const float& meshscale  = 1.0);

bool LoadOBJtoRenderEdgeMesh(
    const char* filename,
    vertex** const vertdata,
    uint32_t& vertsize,
    uint32_t** const ilistdata,
    uint32_t& ilistsize,
    const fvec3& meshoffset = fvec3(0.0),
    const float& meshscale  = 1.0);

bool LoadOBJtoPhysicsTriangleMesh(
    const char* filename,
    fvec3** const vertdata,
    uint32_t& vertsize,
    uint32_t** const ilistdata,
    uint32_t& ilistsize,
    const fvec3& meshoffset = fvec3(0.0),
    const float& meshscale  = 1.0);
