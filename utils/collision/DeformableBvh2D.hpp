#pragma once

#include "utils/mathfunc/mathfunc.hpp"

#include <vector>

class DeformableBvh2D;

struct ContactFeature {
    public:
	uint32_t va0, va1, va2, vb0, vb1, vb2;
	ContactFeature(
	    uint32_t va0,
	    uint32_t va1,
	    uint32_t va2,
	    uint32_t vb0,
	    uint32_t vb1,
	    uint32_t vb2)
	    : va0(va0)
	    , va1(va1)
	    , va2(va2)
	    , vb0(vb0)
	    , vb1(vb1)
	    , vb2(vb2)
	{
	}
};

class DeformableBvh2DNode {
    public:
	const DeformableBvh2D& Root;
	uint32_t Type;
	//Inner : 0
	//Leaf  : 1

	fvec2 center;
	float Lengthx, Lengthy;

	//Innerの場合のみ意味を持つ
	DeformableBvh2DNode* RightChild;
	DeformableBvh2DNode* LeftChild;

	//Leafの場合のみ意味を持つ
	uint32_t index0, index1, index2;

	////////////////////

	DeformableBvh2DNode(
	    const DeformableBvh2D& root,
	    const uint32_t* const elementdata,
	    const uint32_t elementsize);

	void UpdateBvhNode();
};

void DetectCollisionNode(std::vector<ContactFeature>& ContactList, const DeformableBvh2DNode* const RNode, const DeformableBvh2DNode* const LNode);

class DeformableBvh2D {
    public:
	const fvec2* const vertdata;
	const uint32_t vertsize;
	const uint32_t* const elementdata;
	const uint32_t elementsize;
	const uint32_t* const VtoElist;
	const uint32_t* const VtoEind;

	DeformableBvh2DNode* RootNode;

	DeformableBvh2D(
	    const fvec2* const vertdata,
	    const uint32_t vertsize,
	    const uint32_t* const elementdata,
	    const uint32_t elementsize,
	    const uint32_t* const VtoElist,
	    const uint32_t* const VtoEind);

	void UpdateBvh();

	void DetectSelfCollision(std::vector<ContactFeature>& ContactList);
};
