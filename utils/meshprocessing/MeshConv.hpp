#pragma once

#include <cstdint>

#include "utils/mathfunc/mathfunc.hpp"

#include "opengl/drawobject.hpp"

void ConvertPTMtoREM(
    const fvec3* const Pvdata,
    const uint32_t& Pvsize,
    const uint32_t* const Pidata,
    const uint32_t& Pisize,
    linevertarray& varray);

void ConvertPTMtoPEM(
    const uint32_t& Pvsize,
    const uint32_t* const Pidata,
    const uint32_t& Pisize,
    uint32_t** const edgedata,
    uint32_t& edgesize);

void ConvertEVtoVE(
    const uint32_t vertsize,
    const uint32_t* const elementlist,
    const uint32_t elementsize,
    uint32_t** const elsup_index,
    uint32_t** const elsup);

void ExtractBoundaryTetrahedra(
    const uint32_t vertsize,
    const uint32_t* const elementlist,
    const uint32_t elementsize,
    const uint32_t* const VtoElist,
    const uint32_t* const VtoEind,
    const uint32_t* const trilist,
    const uint32_t trisize,
    uint32_t** const boundaryelementlist,
    uint32_t& boundaryelementsize);

void ExtractElementVertex(
    const uint32_t vertsize,
    const uint32_t* const trilist,
    const uint32_t trisize,
    const uint32_t* const VtoTlist,
    const uint32_t* const VtoTind,
    uint32_t** const TritoVertlist,
    uint32_t& TritoVertsize,
    uint32_t** const nonredTlist);
