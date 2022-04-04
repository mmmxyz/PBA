#pragma once

struct ClosestPair {
	float dist;
	int32_t status;
};

float DistLinePoint(
    const fvec2& a0,
    const fvec2& a1,
    const fvec2& b);

ClosestPair DistSegmentPoint(
    const fvec2& a0,
    const fvec2& a1,
    const fvec2& b);

ClosestPair DistTrianglePoint(
    const fvec2& a0,
    const fvec2& a1,
    const fvec2& a2,
    const fvec2& b);

ClosestPair DistTriangleTriangle(
    const fvec2& a0,
    const fvec2& a1,
    const fvec2& a2,
    const fvec2& b0,
    const fvec2& b1,
    const fvec2& b2);

bool Is_CollideTriangleTriangle(
    const fvec2& a0,
    const fvec2& a1,
    const fvec2& a2,
    const fvec2& b0,
    const fvec2& b1,
    const fvec2& b2);
