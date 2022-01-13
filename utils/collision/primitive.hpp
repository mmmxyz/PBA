#pragma once

#include "utils/mathfunc/mathfunc.hpp"

//todo ファイル名をprimitive collisionなどに変える

struct ClosestDV {
	float dist;
	fvec3 v;
};

struct ClosestDVV {
	float dist;
	fvec3 v0;
	fvec3 v1;
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
