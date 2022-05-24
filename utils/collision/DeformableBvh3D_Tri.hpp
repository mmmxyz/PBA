#pragma once

#include "utils/mathfunc/mathfunc.hpp"

#include <vector>

class DeformableBvh3DTri;

struct ContactFeature3DTri {
    public:
	uint32_t va0, va1, va2, vb0, vb1, vb2;
	ContactFeature3DTri(
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

class DeformableBvh3DTriNode {
    public:
	const DeformableBvh3DTri& Root;
	uint32_t Type;
	//Inner : 0
	//Leaf  : 1

	fvec3 center;
	float Lengthx, Lengthy, Lengthz;

	//Innerの場合のみ意味を持つ
	DeformableBvh3DTriNode* RightChild;
	DeformableBvh3DTriNode* LeftChild;

	//Leafの場合のみ意味を持つ
	uint32_t index0, index1, index2;

	////////////////////

	DeformableBvh3DTriNode(
	    const DeformableBvh3DTri& root,
	    const uint32_t* const elementdata,
	    const uint32_t elementsize);

	void UpdateBvhNode();
};

void DetectCollisionNode(std::vector<ContactFeature3DTri>& ContactList, const DeformableBvh3DTriNode* const RNode, const DeformableBvh3DTriNode* const LNode);

void DetectSemiExternalCollisionNode(std::vector<ContactFeature3DTri>& ContactList, const DeformableBvh3DTriNode* const RNode, const DeformableBvh3DTriNode* const LNode);

void DetectExternalCollisionNode(std::vector<ContactFeature3DTri>& ContactList, const DeformableBvh3DTriNode* const RNode, const DeformableBvh3DTriNode* const LNode);

class DeformableBvh3DTri {
    public:
	const fvec3* const vertdata;
	const uint32_t vertsize;
	const uint32_t* const elementdata;
	const uint32_t elementsize;
	const uint32_t* const VtoElist;
	const uint32_t* const VtoEind;

	DeformableBvh3DTriNode* RootNode;

	DeformableBvh3DTri(
	    const fvec3* const vertdata,
	    const uint32_t vertsize,
	    const uint32_t* const elementdata,
	    const uint32_t elementsize,
	    const uint32_t* const VtoElist,
	    const uint32_t* const VtoEind);

	void UpdateBvh();

	void DetectSelfCollision(std::vector<ContactFeature3DTri>& ContactList);
};

void DetectExternalCollision(const DeformableBvh3DTri& RightBvh, const DeformableBvh3DTri& LeftBvh, std::vector<ContactFeature3DTri>& ContactList);

void DetectCollision(const DeformableBvh3DTriNode* const RightBvhNode, const DeformableBvh3DTriNode* const LeftBvhNode, std::vector<ContactFeature3DTri>& ContactList);

void DetectExternalCollision(const DeformableBvh3DTriNode* const RightBvhNode, const DeformableBvh3DTriNode* const LeftBvhNode, std::vector<ContactFeature3DTri>& ContactList);

void DetectSemiExternalCollision(const DeformableBvh3DTriNode* const RightBvhNode, const DeformableBvh3DTriNode* const LeftBvhNode, std::vector<ContactFeature3DTri>& ContactList);
