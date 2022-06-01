#include <vector>
#include <cmath>
#include <iostream>

#include "utils/collision/DeformableBvh3D_Tri.hpp"

#include "utils/collision/primitive2d.hpp"
#include "utils/collision/basic3d.hpp"

constexpr float Epsilon = 0.000000000001;

DeformableBvh3DTriNode::DeformableBvh3DTriNode(
    const DeformableBvh3DTri& root,
    const uint32_t* const elementdata,
    const uint32_t elementsize)
    : Root(root)
{
	if (elementsize == 3)
		this->Type = 1;
	else
		this->Type = 0;

	float maxX, minX, maxY, minY, maxZ, minZ;
	maxX = root.vertdata[elementdata[0]].x;
	minX = root.vertdata[elementdata[0]].x;
	maxY = root.vertdata[elementdata[0]].y;
	minY = root.vertdata[elementdata[0]].y;
	maxZ = root.vertdata[elementdata[0]].z;
	minZ = root.vertdata[elementdata[0]].z;

	for (uint32_t i = 0; i < elementsize; i++) {
		const fvec3& v = root.vertdata[elementdata[i]];

		if (v.x < minX)
			minX = v.x;
		if (maxX < v.x)
			maxX = v.x;

		if (v.y < minY)
			minY = v.y;
		if (maxY < v.y)
			maxY = v.y;

		if (v.z < minZ)
			minZ = v.z;
		if (maxZ < v.z)
			maxZ = v.z;
	}

	this->center  = fvec3((maxX + minX) / 2.0, (maxY + minY) / 2.0, (maxZ + minZ) / 2.0);
	this->Lengthx = std::max(maxX - minX, Epsilon);
	this->Lengthy = std::max(maxY - minY, Epsilon);
	this->Lengthz = std::max(maxZ - minZ, Epsilon);

	if (this->Type == 1) {
		this->index0 = elementdata[0];
		this->index1 = elementdata[1];
		this->index2 = elementdata[2];

		this->RightChild = nullptr;
		this->LeftChild	 = nullptr;

	} else {
		std::vector<float> OrderedCoord;

		uint32_t TetSize = elementsize / 3;
		OrderedCoord.reserve(TetSize);

		if (elementsize == 6) {
			this->RightChild = new DeformableBvh3DTriNode(root, elementdata, 3);
			this->LeftChild	 = new DeformableBvh3DTriNode(root, elementdata + 3, 3);
			return;
		} else {
			std::vector<uint32_t> RightElements;
			RightElements.reserve(elementsize);
			std::vector<uint32_t> LeftElements;
			LeftElements.reserve(elementsize);

			//lengthx == lengthy のときfalseになってやばいことになる
			if (this->Lengthx > this->Lengthy - Epsilon && this->Lengthx - Epsilon > this->Lengthz) {

				//x座標で分離
				for (uint32_t i = 0; i < TetSize; i++) {
					float cm = root.vertdata[elementdata[3 * i + 0]].x + root.vertdata[elementdata[3 * i + 1]].x + root.vertdata[elementdata[3 * i + 2]].x;
					OrderedCoord.insert(std::lower_bound(OrderedCoord.begin(), OrderedCoord.end(), cm), cm);
				}

				float SepPlane = OrderedCoord[TetSize / 2];

				for (uint32_t i = 0; i < TetSize; i++) {
					const fvec3& v0 = root.vertdata[elementdata[3 * i + 0]];
					const fvec3& v1 = root.vertdata[elementdata[3 * i + 1]];
					const fvec3& v2 = root.vertdata[elementdata[3 * i + 2]];

					float cm = v0.x + v1.x + v2.x;

					if (cm < SepPlane) {
						RightElements.emplace_back(elementdata[3 * i + 0]);
						RightElements.emplace_back(elementdata[3 * i + 1]);
						RightElements.emplace_back(elementdata[3 * i + 2]);
					} else {
						LeftElements.emplace_back(elementdata[3 * i + 0]);
						LeftElements.emplace_back(elementdata[3 * i + 1]);
						LeftElements.emplace_back(elementdata[3 * i + 2]);
					}
				}
			} else if (this->Lengthy > this->Lengthz - Epsilon && this->Lengthy > this->Lengthx - Epsilon) {

				//y座標で分離
				for (uint32_t i = 0; i < TetSize; i++) {
					float cm = root.vertdata[elementdata[3 * i + 0]].y + root.vertdata[elementdata[3 * i + 1]].y + root.vertdata[elementdata[3 * i + 2]].y;
					OrderedCoord.insert(std::lower_bound(OrderedCoord.begin(), OrderedCoord.end(), cm), cm);
				}

				float SepPlane = OrderedCoord[TetSize / 2];

				for (uint32_t i = 0; i < TetSize; i++) {
					const fvec3& v0 = root.vertdata[elementdata[3 * i + 0]];
					const fvec3& v1 = root.vertdata[elementdata[3 * i + 1]];
					const fvec3& v2 = root.vertdata[elementdata[3 * i + 2]];

					float cm = v0.y + v1.y + v2.y;

					if (cm < SepPlane) {
						RightElements.emplace_back(elementdata[3 * i + 0]);
						RightElements.emplace_back(elementdata[3 * i + 1]);
						RightElements.emplace_back(elementdata[3 * i + 2]);
					} else {
						LeftElements.emplace_back(elementdata[3 * i + 0]);
						LeftElements.emplace_back(elementdata[3 * i + 1]);
						LeftElements.emplace_back(elementdata[3 * i + 2]);
					}
				}
			} else {

				//z座標で分離
				for (uint32_t i = 0; i < TetSize; i++) {
					float cm = root.vertdata[elementdata[3 * i + 0]].z + root.vertdata[elementdata[3 * i + 1]].z + root.vertdata[elementdata[3 * i + 2]].z;
					OrderedCoord.insert(std::lower_bound(OrderedCoord.begin(), OrderedCoord.end(), cm), cm);
				}

				float SepPlane = OrderedCoord[TetSize / 2];

				for (uint32_t i = 0; i < TetSize; i++) {
					const fvec3& v0 = root.vertdata[elementdata[3 * i + 0]];
					const fvec3& v1 = root.vertdata[elementdata[3 * i + 1]];
					const fvec3& v2 = root.vertdata[elementdata[3 * i + 2]];

					float cm = v0.z + v1.z + v2.z;

					if (cm < SepPlane) {
						RightElements.emplace_back(elementdata[3 * i + 0]);
						RightElements.emplace_back(elementdata[3 * i + 1]);
						RightElements.emplace_back(elementdata[3 * i + 2]);
					} else {
						LeftElements.emplace_back(elementdata[3 * i + 0]);
						LeftElements.emplace_back(elementdata[3 * i + 1]);
						LeftElements.emplace_back(elementdata[3 * i + 2]);
					}
				}
			}

			if (RightElements.size() == 0) {
				RightElements.emplace_back(LeftElements[LeftElements.size() - 3]);
				RightElements.emplace_back(LeftElements[LeftElements.size() - 2]);
				RightElements.emplace_back(LeftElements[LeftElements.size() - 1]);
				LeftElements.pop_back();
				LeftElements.pop_back();
				LeftElements.pop_back();
			}
			if (LeftElements.size() == 0) {
				LeftElements.emplace_back(RightElements[RightElements.size() - 3]);
				LeftElements.emplace_back(RightElements[RightElements.size() - 2]);
				LeftElements.emplace_back(RightElements[RightElements.size() - 1]);
				RightElements.pop_back();
				RightElements.pop_back();
				RightElements.pop_back();
			}

			this->RightChild = new DeformableBvh3DTriNode(root, RightElements.data(), RightElements.size());
			this->LeftChild	 = new DeformableBvh3DTriNode(root, LeftElements.data(), LeftElements.size());
		}
	}
}

void DeformableBvh3DTriNode::UpdateBvhNode()
{

	if (this->Type == 1) {
		const fvec3& v0 = this->Root.vertdata[this->index0];
		const fvec3& v1 = this->Root.vertdata[this->index1];
		const fvec3& v2 = this->Root.vertdata[this->index2];

		float maxx = std::max(std::max(v0.x, v1.x), v2.x);
		float minx = std::min(std::min(v0.x, v1.x), v2.x);

		float maxy = std::max(std::max(v0.y, v1.y), v2.y);
		float miny = std::min(std::min(v0.y, v1.y), v2.y);

		float maxz = std::max(std::max(v0.z, v1.z), v2.z);
		float minz = std::min(std::min(v0.z, v1.z), v2.z);

		this->center  = fvec3((maxx + minx) / 2.0, (maxy + miny) / 2.0, (maxz + minz) / 2.0);
		this->Lengthx = std::max(maxx - minx, Epsilon);
		this->Lengthy = std::max(maxy - miny, Epsilon);
		this->Lengthz = std::max(maxz - minz, Epsilon);

	} else {
		RightChild->UpdateBvhNode();
		LeftChild->UpdateBvhNode();

		float maxx, minx, maxy, miny, maxz, minz;
		if (RightChild->center.x + (RightChild->Lengthx) / 2.0 > LeftChild->center.x + (LeftChild->Lengthx) / 2.0)
			maxx = RightChild->center.x + (RightChild->Lengthx) / 2.0;
		else
			maxx = LeftChild->center.x + (LeftChild->Lengthx) / 2.0;

		if (RightChild->center.x - (RightChild->Lengthx) / 2.0 < LeftChild->center.x - (LeftChild->Lengthx) / 2.0)
			minx = RightChild->center.x - (RightChild->Lengthx) / 2.0;
		else
			minx = LeftChild->center.x - (LeftChild->Lengthx) / 2.0;

		if (RightChild->center.y + (RightChild->Lengthy) / 2.0 > LeftChild->center.y + (LeftChild->Lengthy) / 2.0)
			maxy = RightChild->center.y + (RightChild->Lengthy) / 2.0;
		else
			maxy = LeftChild->center.y + (LeftChild->Lengthy) / 2.0;

		if (RightChild->center.y - (RightChild->Lengthy) / 2.0 < LeftChild->center.y - (LeftChild->Lengthy) / 2.0)
			miny = RightChild->center.y - (RightChild->Lengthy) / 2.0;
		else
			miny = LeftChild->center.y - (LeftChild->Lengthy) / 2.0;

		if (RightChild->center.z + (RightChild->Lengthz) / 2.0 > LeftChild->center.z + (LeftChild->Lengthz) / 2.0)
			maxz = RightChild->center.z + (RightChild->Lengthz) / 2.0;
		else
			maxz = LeftChild->center.z + (LeftChild->Lengthz) / 2.0;

		if (RightChild->center.z - (RightChild->Lengthz) / 2.0 < LeftChild->center.z - (LeftChild->Lengthz) / 2.0)
			minz = RightChild->center.z - (RightChild->Lengthz) / 2.0;
		else
			minz = LeftChild->center.z - (LeftChild->Lengthz) / 2.0;

		this->center  = fvec3((maxx + minx) / 2.0, (maxy + miny) / 2.0, (maxz + minz) / 2.0);
		this->Lengthx = std::max(maxx - minx, Epsilon);
		this->Lengthy = std::max(maxy - miny, Epsilon);
		this->Lengthz = std::max(maxz - minz, Epsilon);
	}
}

bool Is_CollideNodeAABB(const DeformableBvh3DTriNode* const RNode, const DeformableBvh3DTriNode* const LNode)
{
	//sat

	//めり込みが小さいときは無視(隣り合うやつどうしを省きたい)

	if (RNode->center.x + RNode->Lengthx / 2.0 - Epsilon < LNode->center.x - LNode->Lengthx / 2.0)
		return false;
	if (RNode->center.x - RNode->Lengthx / 2.0 + Epsilon > LNode->center.x + LNode->Lengthx / 2.0)
		return false;

	if (RNode->center.y + RNode->Lengthy / 2.0 - Epsilon < LNode->center.y - LNode->Lengthy / 2.0)
		return false;
	if (RNode->center.y - RNode->Lengthy / 2.0 + Epsilon > LNode->center.y + LNode->Lengthy / 2.0)
		return false;

	if (RNode->center.z + RNode->Lengthz / 2.0 - Epsilon < LNode->center.z - LNode->Lengthz / 2.0)
		return false;
	if (RNode->center.z - RNode->Lengthz / 2.0 + Epsilon > LNode->center.z + LNode->Lengthz / 2.0)
		return false;

	return true;
}

void DetectCollisionNode(std::vector<ContactFeature3DTri>& ContactList, const DeformableBvh3DTriNode* const RNode, const DeformableBvh3DTriNode* const LNode)
{
	//RNodeとLNodeは別のオブジェクトであることが保証される
	//RNodeとLNodeはInnerでもLeafでもありえる

	if (RNode->Type == 0 && LNode->Type == 0) {

		//RNode内部で起こる衝突
		DetectCollisionNode(ContactList, RNode->RightChild, RNode->LeftChild);
		DetectCollisionNode(ContactList, LNode->RightChild, LNode->LeftChild);

		if (Is_CollideNodeAABB(RNode, LNode)) {
			DetectCollisionNode(ContactList, RNode->RightChild, LNode->LeftChild);
			DetectCollisionNode(ContactList, RNode->RightChild, LNode->RightChild);
			DetectCollisionNode(ContactList, RNode->LeftChild, LNode->LeftChild);
			DetectCollisionNode(ContactList, RNode->LeftChild, LNode->RightChild);
		}
	} else if (RNode->Type == 0 && LNode->Type == 1) {
		//std::cout << "01" << std::endl;
		DetectCollisionNode(ContactList, RNode->RightChild, RNode->LeftChild);

		if (Is_CollideNodeAABB(RNode, LNode)) {
			DetectCollisionNode(ContactList, RNode->RightChild, LNode);
			DetectCollisionNode(ContactList, RNode->LeftChild, LNode);
		}
	} else if (RNode->Type == 1 && LNode->Type == 0) {
		//std::cout << "10" << std::endl;
		DetectCollisionNode(ContactList, LNode->RightChild, LNode->LeftChild);

		if (Is_CollideNodeAABB(RNode, LNode)) {
			DetectCollisionNode(ContactList, RNode, LNode->LeftChild);
			DetectCollisionNode(ContactList, RNode, LNode->RightChild);
		}
	} else if (RNode->Type == 1 && LNode->Type == 1) {
		//std::cout << "11" << std::endl;
		if (Is_CollideNodeAABB(RNode, LNode)) {

			//隣接する要素同士の場合を除外
			//RNodeとLNodeが同一のrootを持つことを前提にしている
			//異なるrootに属するならそもそもこの作業は必要ない

			//if (
			//    RNode->index0 == LNode->index0 || RNode->index0 == LNode->index1 || RNode->index0 == LNode->index2 || RNode->index1 == LNode->index0 || RNode->index1 == LNode->index1 || RNode->index1 == LNode->index2 || RNode->index2 == LNode->index0 || RNode->index2 == LNode->index1 || RNode->index2 == LNode->index2)
			//	return;

			if (Is_CollideTriTri(
				RNode->Root.vertdata[RNode->index0],
				RNode->Root.vertdata[RNode->index1],
				RNode->Root.vertdata[RNode->index2],
				LNode->Root.vertdata[LNode->index0],
				LNode->Root.vertdata[LNode->index1],
				LNode->Root.vertdata[LNode->index2]))
				ContactList.emplace_back(
				    RNode->index0,
				    RNode->index1,
				    RNode->index2,
				    LNode->index0,
				    LNode->index1,
				    LNode->index2);
		}
	}
}
void DetectSemiExternalCollisionNode(std::vector<ContactFeature3DTri>& ContactList, const DeformableBvh3DTriNode* const RNode, const DeformableBvh3DTriNode* const LNode)
{
	if (RNode->Type == 0 && LNode->Type == 0) {
		//std::cout << "00" << std::endl;
		if (Is_CollideNodeAABB(RNode, LNode)) {
			DetectCollisionNode(ContactList, RNode->RightChild, LNode->LeftChild);
			DetectCollisionNode(ContactList, RNode->RightChild, LNode->RightChild);
			DetectCollisionNode(ContactList, RNode->LeftChild, LNode->LeftChild);
			DetectCollisionNode(ContactList, RNode->LeftChild, LNode->RightChild);
		}
	} else if (RNode->Type == 0 && LNode->Type == 1) {
		//std::cout << "01" << std::endl;
		if (Is_CollideNodeAABB(RNode, LNode)) {
			DetectCollisionNode(ContactList, RNode->RightChild, LNode);
			DetectCollisionNode(ContactList, RNode->LeftChild, LNode);
		}
	} else if (RNode->Type == 1 && LNode->Type == 0) {
		//std::cout << "10" << std::endl;
		if (Is_CollideNodeAABB(RNode, LNode)) {
			DetectCollisionNode(ContactList, RNode, LNode->LeftChild);
			DetectCollisionNode(ContactList, RNode, LNode->RightChild);
		}
	} else if (RNode->Type == 1 && LNode->Type == 1) {
		//std::cout << "11" << std::endl;
		if (Is_CollideNodeAABB(RNode, LNode)) {

			//if (
			//    RNode->index0 == LNode->index0 || RNode->index0 == LNode->index1 || RNode->index0 == LNode->index2 || RNode->index1 == LNode->index0 || RNode->index1 == LNode->index1 || RNode->index1 == LNode->index2 || RNode->index2 == LNode->index0 || RNode->index2 == LNode->index1 || RNode->index2 == LNode->index2)
			//	return;

			if (Is_CollideTriTri(
				RNode->Root.vertdata[RNode->index0],
				RNode->Root.vertdata[RNode->index1],
				RNode->Root.vertdata[RNode->index2],
				LNode->Root.vertdata[LNode->index0],
				LNode->Root.vertdata[LNode->index1],
				LNode->Root.vertdata[LNode->index2]))
				ContactList.emplace_back(
				    RNode->index0,
				    RNode->index1,
				    RNode->index2,
				    LNode->index0,
				    LNode->index1,
				    LNode->index2);
		}
	}
}

void DetectExternalCollisionNode(std::vector<ContactFeature3DTri>& ContactList, const DeformableBvh3DTriNode* const RNode, const DeformableBvh3DTriNode* const LNode)
{

	//RNodeとLNodeは別のオブジェクトであることが保証される
	//RNodeとLNodeはInnerでもLeafでもありえる

	//ここで再帰的に呼ぶ関数をDetectCollisionNodeとしてやばいことになった．

	if (RNode->Type == 0 && LNode->Type == 0) {
		//std::cout << "00" << std::endl;
		if (Is_CollideNodeAABB(RNode, LNode)) {
			DetectExternalCollisionNode(ContactList, RNode->RightChild, LNode->LeftChild);
			DetectExternalCollisionNode(ContactList, RNode->RightChild, LNode->RightChild);
			DetectExternalCollisionNode(ContactList, RNode->LeftChild, LNode->LeftChild);
			DetectExternalCollisionNode(ContactList, RNode->LeftChild, LNode->RightChild);
		}
	} else if (RNode->Type == 0 && LNode->Type == 1) {
		//std::cout << "01" << std::endl;
		if (Is_CollideNodeAABB(RNode, LNode)) {
			DetectExternalCollisionNode(ContactList, RNode->RightChild, LNode);
			DetectExternalCollisionNode(ContactList, RNode->LeftChild, LNode);
		}
	} else if (RNode->Type == 1 && LNode->Type == 0) {
		//std::cout << "10" << std::endl;
		if (Is_CollideNodeAABB(RNode, LNode)) {
			DetectExternalCollisionNode(ContactList, RNode, LNode->LeftChild);
			DetectExternalCollisionNode(ContactList, RNode, LNode->RightChild);
		}
	} else if (RNode->Type == 1 && LNode->Type == 1) {
		//std::cout << "11" << std::endl;
		if (Is_CollideNodeAABB(RNode, LNode)) {

			//if (
			//    RNode->index0 == LNode->index0 || RNode->index0 == LNode->index1 || RNode->index0 == LNode->index2 || RNode->index1 == LNode->index0 || RNode->index1 == LNode->index1 || RNode->index1 == LNode->index2 || RNode->index2 == LNode->index0 || RNode->index2 == LNode->index1 || RNode->index2 == LNode->index2)
			//	return;

			if (Is_CollideTriTri(
				RNode->Root.vertdata[RNode->index0],
				RNode->Root.vertdata[RNode->index1],
				RNode->Root.vertdata[RNode->index2],
				LNode->Root.vertdata[LNode->index0],
				LNode->Root.vertdata[LNode->index1],
				LNode->Root.vertdata[LNode->index2]))
				ContactList.emplace_back(
				    RNode->index0,
				    RNode->index1,
				    RNode->index2,
				    LNode->index0,
				    LNode->index1,
				    LNode->index2);
		}
	}
}

//////////////////////////////////////////

bool Is_CollideAABBRay(const DeformableBvh3DTriNode* const Node, const fvec3& v0, const fvec3& v1)
{
	//v0 != v1

	//各軸において，rayはAABBとmax・minの平面と交差する
	//AABBの平面をray上に射影した際の座標を比較する

	float tmax = 999999999999.9999, tmin = -999999999999.9999;

	fvec3 dir = (v1 - v0).normalize();

	{
		float xmax = Node->center.x + Node->Lengthx / 2.0;
		float xmin = Node->center.x - Node->Lengthx / 2.0;

		if (std::abs(dir.x) < Epsilon) {
		} else if (dir.x > 0.0) {
			tmax = (xmax - v0.x) / dir.x;
			tmin = (xmin - v0.x) / dir.x;
		} else {
			tmax = (xmin - v0.x) / dir.x;
			tmin = (xmax - v0.x) / dir.x;
		}
	}

	{
		float ymax = Node->center.y + Node->Lengthy / 2.0;
		float ymin = Node->center.y - Node->Lengthy / 2.0;

		if (std::abs(dir.y) < Epsilon) {
		} else if (dir.y > 0.0) {
			tmax = std::min((ymax - v0.y) / dir.y, tmax);
			tmin = std::max((ymin - v0.y) / dir.y, tmin);
		} else {
			tmax = std::min((ymin - v0.y) / dir.y, tmax);
			tmin = std::max((ymax - v0.y) / dir.y, tmin);
		}
	}

	{
		float zmax = Node->center.z + Node->Lengthz / 2.0;
		float zmin = Node->center.z - Node->Lengthz / 2.0;

		if (std::abs(dir.z) < Epsilon) {
		} else if (dir.z > 0.0) {
			tmax = std::min((zmax - v0.z) / dir.z, tmax);
			tmin = std::max((zmin - v0.z) / dir.z, tmin);
		} else {
			tmax = std::min((zmin - v0.z) / dir.z, tmax);
			tmin = std::max((zmax - v0.z) / dir.z, tmin);
		}
	}

	if (tmax < tmin - Epsilon)
		return false;

	if (0.0 < tmax)
		return true;

	//back intersection
	return false;
}

DeformableBvh3DTri::IntersecInfo DetectBvhTriRayIntersection(const DeformableBvh3DTriNode* const Node, const fvec3& v0, const fvec3& v1)
{
	if (Is_CollideAABBRay(Node, v0, v1)) {
		if (Node->Type == 0) {
			auto IR = DetectBvhTriRayIntersection(Node->RightChild, v0, v1);
			auto IL = DetectBvhTriRayIntersection(Node->LeftChild, v0, v1);

			if (IR.elementnum != -1 && IL.elementnum != -1) {
				if (IR.pram[3] > IL.pram[3])
					return IR;
				else
					return IL;
			} else if (IR.elementnum != -1) {
				return IR;
			} else if (IL.elementnum != -1) {
				return IL;
			} else {
				return { 0.0, 0.0, 0.0, 0.0, -1 };
			}

		} else {

			auto hoge = IntersectTriangleRay(
			    Node->Root.vertdata[Node->index0],
			    Node->Root.vertdata[Node->index1],
			    Node->Root.vertdata[Node->index2],
			    v0, v1);
			if (hoge.Is_intersect) {
				uint32_t ei = (Node->index0) / 3;
				return { hoge.p[0], hoge.p[1], hoge.p[2], hoge.p[3], ei };
			} else {
				return { 0.0, 0.0, 0.0, 0.0, -1 };
			}
		}
	} else
		return { 0.0, 0.0, 0.0, 0.0, -1 };
}

//////////////////////////////////////////

DeformableBvh3DTri::DeformableBvh3DTri(
    const fvec3* const vertdata,
    const uint32_t vertsize,
    const uint32_t* const elementdata,
    const uint32_t elementsize,
    const uint32_t* const VtoElist,
    const uint32_t* const VtoEind)
    : vertdata(vertdata)
    , vertsize(vertsize)
    , elementdata(elementdata)
    , elementsize(elementsize)
    , VtoElist(VtoElist)
    , VtoEind(VtoEind)
{
	RootNode = new DeformableBvh3DTriNode((*this), elementdata, elementsize);
}

void DeformableBvh3DTri::UpdateBvh()
{
	RootNode->UpdateBvhNode();
}

void DeformableBvh3DTri::DetectSelfCollision(std::vector<ContactFeature3DTri>& ContactList)
{
	ContactList.clear();
	//ContactList.reserve(elementsize);
	DetectCollisionNode(ContactList, RootNode->RightChild, RootNode->LeftChild);
}

void DetectExternalCollision(const DeformableBvh3DTri& RightBvh, const DeformableBvh3DTri& LeftBvh, std::vector<ContactFeature3DTri>& ContactList)
{
	ContactList.clear();
	//ContactList.reserve(RightBvh.elementsize);
	DetectExternalCollisionNode(ContactList, RightBvh.RootNode, LeftBvh.RootNode);
}

DeformableBvh3DTri::IntersecInfo DeformableBvh3DTri::RayIntersection(const fvec3& v0, const fvec3& v1)
{
	return DetectBvhTriRayIntersection(RootNode, v0, v1);
}
