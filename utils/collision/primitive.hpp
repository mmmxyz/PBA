#pragma once

#include "utils/mathfunc/mathfunc.hpp"

struct ClosestDV {
	float dist;
	fvec3 v;
};

struct ClosestDVV {
	float dist;
	fvec3 v0;
	fvec3 v1;
};

struct ClosestPair3D {
	float dist;
	int32_t status;
};

ClosestDV DistLinePoint(const fvec3& a0, const fvec3& a1, const fvec3& b, float* const pt = nullptr);

ClosestDV DistSegmentPoint(const fvec3& a0, const fvec3& a1, const fvec3& b);

ClosestDVV DistLineLine(const fvec3& a0, const fvec3& a1, const fvec3& b0, const fvec3& b1, float* const pt = nullptr, float* const ps = nullptr);

ClosestDVV DistLineSegment(const fvec3& a0, const fvec3& a1, const fvec3& b0, const fvec3& b1, float* const pt = nullptr, float* const ps = nullptr);

ClosestDVV DistSegmentSegment(const fvec3& a0, const fvec3& a1, const fvec3& b0, const fvec3& b1);

ClosestDV DistTrianglePoint(const fvec3& a0, const fvec3& a1, const fvec3& a2, const fvec3& b);

ClosestDVV DistTriangleLine(const fvec3& a0, const fvec3& a1, const fvec3& a2, const fvec3& b0, const fvec3& b1, float* const ps = nullptr);

ClosestDVV DistTriangleSegment(const fvec3& a0, const fvec3& a1, const fvec3& a2, const fvec3& b0, const fvec3& b1);

ClosestDVV DistTriangleTriangle(const fvec3& a0, const fvec3& a1, const fvec3& a2, const fvec3& b0, const fvec3& b1, const fvec3& b2);

bool Is_CollideTetraTetra(
    const fvec3& a0,
    const fvec3& a1,
    const fvec3& a2,
    const fvec3& a3,
    const fvec3& b0,
    const fvec3& b1,
    const fvec3& b2,
    const fvec3& b3);

ClosestPair3D DistTetraTetra(
    const fvec3& a0,
    const fvec3& a1,
    const fvec3& a2,
    const fvec3& a3,
    const fvec3& b0,
    const fvec3& b1,
    const fvec3& b2,
    const fvec3& b3);
