#pragma once

#include <cstdint>
#include "utils/mathfunc/mathfunc.hpp"

//void CubeTetrahedra(
//    const uint32_t Nx,
//    const uint32_t Ny,
//    const uint32_t Nz,
//    const float Lx,
//    const float Ly,
//    const float Lz,
//    fvec3** const vertdata,
//    uint32_t& vertsize,
//    uint32_t** const ilistdata,
//    uint32_t& ilistsize,
//    uint32_t** const facelistdata,
//    uint32_t& facelistsize,
//    uint32_t** const edgelistdata,
//    uint32_t& edgelistsize,
//    const fvec3& meshoffset,
//    const float& meshscale);

void CubeTetrahedra(
    const uint32_t N,
    const float L,
    fvec3** const vertdata,
    uint32_t& vertsize,
    uint32_t** const ilistdata,
    uint32_t& ilistsize,
    uint32_t** const facelistdata,
    uint32_t& facelistsize,
    uint32_t** const edgelistdata,
    uint32_t& edgelistsize,
    const fvec3& bias = fvec3(0.0));

void CylinderTetrahedron(
    const uint32_t N,
    const uint32_t M,
    const uint32_t L,
    const float LengthR,
    const float Lengthr,
    const float Lengthl,
    fvec3** const vertdata,
    uint32_t& vertsize,
    uint32_t** const ilistdata,
    uint32_t& ilistsize,
    uint32_t** const facelistdata,
    uint32_t& facelistsize,
    uint32_t** const edgelistdata,
    uint32_t& edgelistsize,
    const fvec3& bias = fvec3(0.0));

void RectTriangle(
    const uint32_t N,
    const uint32_t M,
    const float Lx,
    const float Ly,
    fvec2** const vertdata,
    uint32_t& vertsize,
    uint32_t** const ilistdata,
    uint32_t& ilistsize,
    uint32_t** const edgelistdata,
    uint32_t& edgelistsize,
    uint32_t** const boundarydata,
    uint32_t& boundarysize,
    const fvec2& bias = fvec2(0.0));

void RectTriangle2(
    const uint32_t N,
    const uint32_t M,
    const float Lx,
    const float Ly,
    fvec2** const vertdata,
    uint32_t& vertsize,
    uint32_t** const ilistdata,
    uint32_t& ilistsize,
    uint32_t** const edgelistdata,
    uint32_t& edgelistsize,
    uint32_t** const boundarydata,
    uint32_t& boundarysize,
    const fvec2& bias = fvec2(0.0));

void ClothFemMesh(
    const uint32_t N,
    const float L,
    fvec3** const vertdata,
    fvec2** const Restvertdata,
    uint32_t& vertsize,
    uint32_t** const ilistdata,
    uint32_t& ilistsize,
    uint32_t** const edgedata,
    uint32_t& edgesize,
    uint32_t** const InnerEdgedata,
    uint32_t& InnerEdgesize,
    fvec4** const InnerEdgeCdata,
    const fvec3& bias = fvec3(0.0));

void ClothFemMesh2(
    const uint32_t N,
    const float L,
    fvec3** const vertdata,
    fvec2** const Restvertdata,
    uint32_t& vertsize,
    uint32_t** const ilistdata,
    uint32_t& ilistsize,
    uint32_t** const edgedata,
    uint32_t& edgesize,
    uint32_t** const InnerEdgedata,
    uint32_t& InnerEdgesize,
    fvec4** const InnerEdgeCdata,
    const fvec3& bias = fvec3(0.0));
