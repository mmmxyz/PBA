#pragma once

#include <deque>

#include "utils/mathfunc/mathfunc.hpp"

#include "opengl/vertarray.hpp"
#include "opengl/drawobject.hpp"

class RigidObbBvh;

class RigidObbBvhNode {
    public:
	const RigidObbBvh& Root;
	uint32_t Type;
	//Inner: 0
	//Leaf : 1

	//Inter nodeのときのみ意味を持つ
	RigidObbBvhNode* RightChild;
	RigidObbBvhNode* LeftChild;

	fvec3 center; //RigidObbBvhのrotq,centerからの相対座標
	float Lengthx, Lengthy, Lengthz;

	//Leaf nodeのときのみ意味を持つ
	uint32_t Index0, Index1, Index2;

	RigidObbBvhNode(const RigidObbBvh& ROB, const fvec3* const Vdata, const uint32_t& Vsize, const uint32_t* const Ldata, const uint32_t& Lsize, std::deque<fvec3>& LVdata, std::deque<uint32_t>& LIData);
};

class RigidObbBvh {
	//剛体の重心、回転と同じ
    public:
	const fquaternion& rotq;
	const fvec3& cm;

	RigidObbBvhNode* RootNode;

	const fvec3* vertdata; //重心からの相対座標を格納している
	uint32_t vertsize;
	const uint32_t* listdata;
	uint32_t listsize;

	linevertarray lva;

	RigidObbBvh(const fquaternion& rotq, const fvec3& cm, const fvec3* const Vdata, const uint32_t& Vsize, const uint32_t* const Ldata, const uint32_t& Lsize);
	RigidObbBvh(const fquaternion& rotq, const fvec3& cm);

	void ConstructBvh(const fvec3* const Vdata, const uint32_t& Vsize, const uint32_t* const Ldata, const uint32_t& Lsize);
};

struct ContactFeature {
	fvec3 normal; // r0を動かす方向
	fvec3 r0rel, r1rel;
};

uint32_t contactontriangle(const fvec3& a0, const fvec3& a1, const fvec3& a2, const fvec3& b0, const fvec3& b1, const fvec3& b2, ContactFeature& CF1, ContactFeature& CF2);

bool Is_CollideRigidObbBvhNode(const RigidObbBvhNode* const ROBNode0, const RigidObbBvhNode* const ROBNode1, const fvec3& r0cm, const fvec3& r1cm, const fmat3& R0, const fmat3& R1, const float extend);

struct TriangleInd {
	uint32_t v00, v01, v02;
	uint32_t v10, v11, v12;
};

void RigidObbBvhPotentialCollisionPair(const RigidObbBvh& ROB0, const RigidObbBvh& ROB1, std::deque<TriangleInd>& PCList, const fvec3& r0cm, const fvec3& r1cm, const fquaternion& rotq0, const fquaternion& rotq1);
