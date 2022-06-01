#pragma once

#include "utils/mathfunc/mathfunc.hpp"

//TODO たぶんinlineにしたほうがよさそう

//衝突の有無を返す

bool Is_CollideTetraTetra(
    const fvec3& a0,
    const fvec3& a1,
    const fvec3& a2,
    const fvec3& a3,
    const fvec3& b0,
    const fvec3& b1,
    const fvec3& b2,
    const fvec3& b3);

bool Is_CollideTriTri(
    const fvec3& a0,
    const fvec3& a1,
    const fvec3& a2,
    const fvec3& b0,
    const fvec3& b1,
    const fvec3& b2);

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

//衝突がある場合，最短距離と座標を返す．
//衝突がない場合，Is_intersect=false

struct IntersectP4 {
	bool Is_intersect;
	float p[4];
};

IntersectP4 IntersectTriangleRay(const fvec3& a0, const fvec3& a1, const fvec3& a2, const fvec3& b0, const fvec3& b1);

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

//衝突の有無によらず，最短距離を計算

struct DistP1 {
	float dist;
	float p;
};

struct DistP2 {
	float dist;
	//float p0;
	//float p1;
	float p[2];
};

struct DistP3 {
	float dist;
	//float p0;
	//float p1;
	//float p2;
	float p[3];
};

struct DistP4 {
	float dist;
	//float p0;
	//float p1;
	//float p2;
	//float p3;
	float p[4];
};

struct DistP6 {
	float dist;
	float p[6];
	//float p0;
	//float p1;
	//float p2;
	//float p3;
	//float p4;
	//float p5;
};

DistP1 DistLinePoint(const fvec3& a0, const fvec3& a1, const fvec3& v);

DistP1 DistSegmentPoint(const fvec3& a0, const fvec3& a1, const fvec3& v);

DistP2 DistLineLine(const fvec3& a0, const fvec3& a1, const fvec3& b0, const fvec3& b1);

DistP2 DistLineSegment(const fvec3& a0, const fvec3& a1, const fvec3& b0, const fvec3& b1);

DistP2 DistSegmentSegment(const fvec3& a0, const fvec3& a1, const fvec3& b0, const fvec3& b1);

DistP3 DistTPlanePoint(const fvec3& a0, const fvec3& a1, const fvec3& a2, const fvec3& b);

DistP3 DistTrianglePoint(const fvec3& a0, const fvec3& a1, const fvec3& a2, const fvec3& b);

DistP4 DistTriangleLine(const fvec3& a0, const fvec3& a1, const fvec3& a2, const fvec3& b0, const fvec3& b1);

DistP4 DistTriangleSegment(const fvec3& a0, const fvec3& a1, const fvec3& a2, const fvec3& b0, const fvec3& b1);

DistP6 DistTriangleTriangle(const fvec3& a0, const fvec3& a1, const fvec3& a2, const fvec3& b0, const fvec3& b1, const fvec3& b2);

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

//衝突がない場合は早期に計算を打ち切る

struct CollisionP6 {
	float depth;
	int32_t status;
	float p[6];
};

struct CollisionP8 {
	float depth;
	int32_t status;
	float p[8];
};

CollisionP6 CollisionTriangleTriangle(const fvec3& a0, const fvec3& a1, const fvec3& a2, const fvec3& b0, const fvec3& b1, const fvec3& b2);

CollisionP8 CollisionTetraTetra(const fvec3& a0, const fvec3& a1, const fvec3& a2, const fvec3& a3, const fvec3& b0, const fvec3& b1, const fvec3& b2, const fvec3& b3);
