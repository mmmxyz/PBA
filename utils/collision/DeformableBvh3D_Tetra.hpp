#pragma once

#include "utils/mathfunc/mathfunc.hpp"

#include <vector>

class DeformableBvh3D;

struct ContactFeature3D {
    public:
	uint32_t va0, va1, va2, va3, vb0, vb1, vb2, vb3;
	ContactFeature3D(
	    uint32_t va0,
	    uint32_t va1,
	    uint32_t va2,
	    uint32_t va3,
	    uint32_t vb0,
	    uint32_t vb1,
	    uint32_t vb2,
	    uint32_t vb3)
	    : va0(va0)
	    , va1(va1)
	    , va2(va2)
	    , va3(va3)
	    , vb0(vb0)
	    , vb1(vb1)
	    , vb2(vb2)
	    , vb3(vb3)
	{
	}
};

class DeformableBvh3DNode {
    public:
	const DeformableBvh3D& Root;
	uint32_t Type;
	//Inner : 0
	//Leaf  : 1

	fvec3 center;
	float Lengthx, Lengthy, Lengthz;

	//Innerの場合のみ意味を持つ
	DeformableBvh3DNode* RightChild;
	DeformableBvh3DNode* LeftChild;

	//Leafの場合のみ意味を持つ
	uint32_t index0, index1, index2, index3;

	////////////////////

	DeformableBvh3DNode(
	    const DeformableBvh3D& root,
	    const uint32_t* const elementdata,
	    const uint32_t elementsize);

	~DeformableBvh3DNode();

	void UpdateBvhNode();
};

void DetectCollisionNode(std::vector<ContactFeature3D>& ContactList, const DeformableBvh3DNode* const RNode, const DeformableBvh3DNode* const LNode);

void DetectSemiExternalCollisionNode(std::vector<ContactFeature3D>& ContactList, const DeformableBvh3DNode* const RNode, const DeformableBvh3DNode* const LNode);

void DetectExternalCollisionNode(std::vector<ContactFeature3D>& ContactList, const DeformableBvh3DNode* const RNode, const DeformableBvh3DNode* const LNode);

class DeformableBvh3D {
    public:
	const fvec3* const vertdata;
	const uint32_t vertsize;
	const uint32_t* const elementdata;
	const uint32_t elementsize;
	const uint32_t* const VtoElist;
	const uint32_t* const VtoEind;

	DeformableBvh3DNode* RootNode;

	DeformableBvh3D(
	    const fvec3* const vertdata,
	    const uint32_t vertsize,
	    const uint32_t* const elementdata,
	    const uint32_t elementsize,
	    const uint32_t* const VtoElist,
	    const uint32_t* const VtoEind);

	void UpdateBvh();

	void DetectSelfCollision(std::vector<ContactFeature3D>& ContactList);

	~DeformableBvh3D();
};

void DetectExternalCollision(const DeformableBvh3D& RightBvh, const DeformableBvh3D& LeftBvh, std::vector<ContactFeature3D>& ContactList);
