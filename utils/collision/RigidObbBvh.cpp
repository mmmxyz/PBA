#include <vector>
#include <deque>
#include <queue>
#include <algorithm>

//#include <iostream>

#include "utils/collision/RigidObbBvh.hpp"

#include "utils/mathfunc/mathfunc.hpp"
#include "utils/collision/geometry.hpp"

constexpr float Epsilon = 0.000001;

RigidObbBvhNode::RigidObbBvhNode(const RigidObbBvh& ROB, const fvec3* const Vdata, const uint32_t& Vsize, const uint32_t* const Ldata, const uint32_t& Lsize, std::deque<fvec3>& LVdata, std::deque<uint32_t>& LIData)
    : Root(ROB)
{

	if (Lsize == 3)
		this->Type = 1;
	else
		this->Type = 0;

	float maxX, minX, maxY, minY, maxZ, minZ;
	maxX = Vdata[Ldata[0]].x;
	minX = Vdata[Ldata[0]].x;
	maxY = Vdata[Ldata[0]].y;
	minY = Vdata[Ldata[0]].y;
	maxZ = Vdata[Ldata[0]].z;
	minZ = Vdata[Ldata[0]].z;

	for (uint32_t i = 1; i < Lsize; i++) {
		if (Vdata[Ldata[i]].x < minX)
			minX = Vdata[Ldata[i]].x;
		if (Vdata[Ldata[i]].x > maxX)
			maxX = Vdata[Ldata[i]].x;

		if (Vdata[Ldata[i]].y < minY)
			minY = Vdata[Ldata[i]].y;
		if (Vdata[Ldata[i]].y > maxY)
			maxY = Vdata[Ldata[i]].y;

		if (Vdata[Ldata[i]].z < minZ)
			minZ = Vdata[Ldata[i]].z;
		if (Vdata[Ldata[i]].z > maxZ)
			maxZ = Vdata[Ldata[i]].z;
	}

	this->center = fvec3(
	    (maxX + minX) / 2.0,
	    (maxY + minY) / 2.0,
	    (maxZ + minZ) / 2.0);

	this->Lengthx = std::max(maxX - minX, Epsilon);
	this->Lengthy = std::max(maxY - minY, Epsilon);
	this->Lengthz = std::max(maxZ - minZ, Epsilon);

	uint32_t IBase = LVdata.size();
	LVdata.emplace_back(minX, minY, maxZ);
	LVdata.emplace_back(maxX, minY, maxZ);
	LVdata.emplace_back(maxX, maxY, maxZ);
	LVdata.emplace_back(minX, maxY, maxZ);

	LVdata.emplace_back(minX, minY, minZ);
	LVdata.emplace_back(maxX, minY, minZ);
	LVdata.emplace_back(maxX, maxY, minZ);
	LVdata.emplace_back(minX, maxY, minZ);

	LIData.emplace_back(IBase + 0);
	LIData.emplace_back(IBase + 1);

	LIData.emplace_back(IBase + 1);
	LIData.emplace_back(IBase + 2);

	LIData.emplace_back(IBase + 2);
	LIData.emplace_back(IBase + 3);

	LIData.emplace_back(IBase + 3);
	LIData.emplace_back(IBase + 0);

	LIData.emplace_back(IBase + 0);
	LIData.emplace_back(IBase + 4);

	LIData.emplace_back(IBase + 1);
	LIData.emplace_back(IBase + 5);

	LIData.emplace_back(IBase + 2);
	LIData.emplace_back(IBase + 6);

	LIData.emplace_back(IBase + 3);
	LIData.emplace_back(IBase + 7);

	LIData.emplace_back(IBase + 4);
	LIData.emplace_back(IBase + 5);

	LIData.emplace_back(IBase + 5);
	LIData.emplace_back(IBase + 6);

	LIData.emplace_back(IBase + 6);
	LIData.emplace_back(IBase + 7);

	LIData.emplace_back(IBase + 7);
	LIData.emplace_back(IBase + 4);

	if (this->Type == 1) {
		this->Index0 = Ldata[0];
		this->Index1 = Ldata[1];
		this->Index2 = Ldata[2];

		this->RightChild = nullptr;
		this->LeftChild	 = nullptr;
	} else {

		//std::cout << this->Lengthx << " ";
		//std::cout << this->Lengthy << " ";
		//std::cout << this->Lengthz << std::endl;

		std::vector<float> OrderedCoord;
		uint32_t FaceSize = Lsize / 3;

		std::vector<uint32_t> RightLVec;
		std::vector<uint32_t> LeftLVec;

		//std::cout << Lsize / 3 << std::endl;
		if (Lsize == 6) {
			this->RightChild = new RigidObbBvhNode(ROB, Vdata, Vsize, Ldata, 3, LVdata, LIData);
			this->LeftChild	 = new RigidObbBvhNode(ROB, Vdata, Vsize, Ldata + 3, 3, LVdata, LIData);
		} else {

			if (this->Lengthx >= this->Lengthy && this->Lengthx >= this->Lengthz) {

				for (uint32_t i = 0; i < FaceSize; i++) {
					float cm = Vdata[Ldata[3 * i + 0]].x + Vdata[Ldata[3 * i + 1]].x + Vdata[Ldata[3 * i + 2]].x;
					OrderedCoord.insert(std::lower_bound(OrderedCoord.begin(), OrderedCoord.end(), cm), cm);
				}

				float SepPlane = OrderedCoord[FaceSize / 2];

				for (uint32_t i = 0; i < FaceSize; i++) {
					fvec3 v0 = Vdata[Ldata[3 * i + 0]];
					fvec3 v1 = Vdata[Ldata[3 * i + 1]];
					fvec3 v2 = Vdata[Ldata[3 * i + 2]];
					float cm = v0.x + v1.x + v2.x;

					if (cm < SepPlane) {
						RightLVec.emplace_back(Ldata[3 * i + 0]);
						RightLVec.emplace_back(Ldata[3 * i + 1]);
						RightLVec.emplace_back(Ldata[3 * i + 2]);
					} else {
						LeftLVec.emplace_back(Ldata[3 * i + 0]);
						LeftLVec.emplace_back(Ldata[3 * i + 1]);
						LeftLVec.emplace_back(Ldata[3 * i + 2]);
					}
				}
			} else if (this->Lengthy >= this->Lengthx && this->Lengthy >= this->Lengthz) {

				for (uint32_t i = 0; i < FaceSize; i++) {
					float cm = Vdata[Ldata[3 * i + 0]].y + Vdata[Ldata[3 * i + 1]].y + Vdata[Ldata[3 * i + 2]].y;
					OrderedCoord.insert(std::lower_bound(OrderedCoord.begin(), OrderedCoord.end(), cm), cm);
				}

				float SepPlane = OrderedCoord[FaceSize / 2];

				for (uint32_t i = 0; i < FaceSize; i++) {
					fvec3 v0 = Vdata[Ldata[3 * i + 0]];
					fvec3 v1 = Vdata[Ldata[3 * i + 1]];
					fvec3 v2 = Vdata[Ldata[3 * i + 2]];
					float cm = v0.y + v1.y + v2.y;

					if (cm < SepPlane) {
						RightLVec.emplace_back(Ldata[3 * i + 0]);
						RightLVec.emplace_back(Ldata[3 * i + 1]);
						RightLVec.emplace_back(Ldata[3 * i + 2]);
					} else {
						LeftLVec.emplace_back(Ldata[3 * i + 0]);
						LeftLVec.emplace_back(Ldata[3 * i + 1]);
						LeftLVec.emplace_back(Ldata[3 * i + 2]);
					}
				}

			} else {

				for (uint32_t i = 0; i < FaceSize; i++) {
					float cm = Vdata[Ldata[3 * i + 0]].z + Vdata[Ldata[3 * i + 1]].z + Vdata[Ldata[3 * i + 2]].z;
					OrderedCoord.insert(std::lower_bound(OrderedCoord.begin(), OrderedCoord.end(), cm), cm);
				}

				float SepPlane = OrderedCoord[FaceSize / 2];

				for (uint32_t i = 0; i < FaceSize; i++) {
					fvec3 v0 = Vdata[Ldata[3 * i + 0]];
					fvec3 v1 = Vdata[Ldata[3 * i + 1]];
					fvec3 v2 = Vdata[Ldata[3 * i + 2]];
					float cm = v0.z + v1.z + v2.z;

					if (cm < SepPlane) {
						RightLVec.emplace_back(Ldata[3 * i + 0]);
						RightLVec.emplace_back(Ldata[3 * i + 1]);
						RightLVec.emplace_back(Ldata[3 * i + 2]);
					} else {
						LeftLVec.emplace_back(Ldata[3 * i + 0]);
						LeftLVec.emplace_back(Ldata[3 * i + 1]);
						LeftLVec.emplace_back(Ldata[3 * i + 2]);
					}
				}
			}

			//std::cout << RightLVec.size() / 3 << " " << LeftLVec.size() / 3 << std::endl;

			if (RightLVec.size() == 0) {
				//ここで1,2,3の順にしてやばいことになった。
				RightLVec.emplace_back(LeftLVec[LeftLVec.size() - 3]);
				RightLVec.emplace_back(LeftLVec[LeftLVec.size() - 2]);
				RightLVec.emplace_back(LeftLVec[LeftLVec.size() - 1]);
				LeftLVec.pop_back();
				LeftLVec.pop_back();
				LeftLVec.pop_back();
			}
			if (LeftLVec.size() == 0) {
				LeftLVec.emplace_back(RightLVec[RightLVec.size() - 3]);
				LeftLVec.emplace_back(RightLVec[RightLVec.size() - 2]);
				LeftLVec.emplace_back(RightLVec[RightLVec.size() - 1]);
				RightLVec.pop_back();
				RightLVec.pop_back();
				RightLVec.pop_back();
			}

			RightChild = new RigidObbBvhNode(ROB, Vdata, Vsize, RightLVec.data(), RightLVec.size(), LVdata, LIData);
			LeftChild  = new RigidObbBvhNode(ROB, Vdata, Vsize, LeftLVec.data(), LeftLVec.size(), LVdata, LIData);
		}
	}
}

RigidObbBvh::RigidObbBvh(const fquaternion& rotq, const fvec3& cm, const fvec3* const Vdata, const uint32_t& Vsize, const uint32_t* const Ldata, const uint32_t& Lsize)
    : rotq(rotq)
    , cm(cm)
    , vertdata(Vdata)
    , vertsize(Vsize)
    , listdata(Ldata)
    , listsize(Lsize)
    , lva()
{
	std::deque<fvec3> LineVertData;
	std::deque<uint32_t> LineIndexData;

	RootNode = new RigidObbBvhNode((*this), vertdata, vertsize, listdata, listsize, LineVertData, LineIndexData);

	vertex* RVdata	 = new vertex[LineVertData.size()];
	uint32_t* RIdata = new uint32_t[LineIndexData.size()];

	//todo iteratorで回すほうが早そう
	for (uint32_t i = 0; i < LineVertData.size(); i++) {
		RVdata[i].position[0] = LineVertData[i].x;
		RVdata[i].position[1] = LineVertData[i].y;
		RVdata[i].position[2] = LineVertData[i].z;

		RVdata[i].color[0] = 1.0;
		RVdata[i].color[1] = 1.0;
		RVdata[i].color[2] = 1.0;
		RVdata[i].color[3] = 1.0;

		RVdata[i].type = 0;
	}

	//todo iteratorで回すほうが早そう
	for (uint32_t i = 0; i < LineIndexData.size(); i++) {
		RIdata[i] = LineIndexData[i];
	}

	lva.resetvertarray(LineVertData.size(), RVdata, LineIndexData.size(), RIdata);
}

RigidObbBvh::RigidObbBvh(const fquaternion& rotq, const fvec3& cm)
    : rotq(rotq)
    , cm(cm)
{
	RootNode = nullptr;
}

void RigidObbBvh::ConstructBvh(const fvec3* const Vdata, const uint32_t& Vsize, const uint32_t* const Ldata, const uint32_t& Lsize)
{

	vertdata = Vdata;
	vertsize = Vsize;
	listdata = Ldata;
	listsize = Lsize;

	std::deque<fvec3> LineVertData;
	std::deque<uint32_t> LineIndexData;

	RootNode = new RigidObbBvhNode((*this), vertdata, vertsize, listdata, listsize, LineVertData, LineIndexData);

	vertex* RVdata	 = new vertex[LineVertData.size()];
	uint32_t* RIdata = new uint32_t[LineIndexData.size()];

	//std::cout << LineVertData.size() << std::endl;
	//std::cout << LineIndexData.size() << std::endl;

	//todo iteratorで回すほうが早そう
	for (uint32_t i = 0; i < LineVertData.size(); i++) {
		RVdata[i].position[0] = LineVertData[i].x;
		RVdata[i].position[1] = LineVertData[i].y;
		RVdata[i].position[2] = LineVertData[i].z;

		RVdata[i].color[0] = 1.0;
		RVdata[i].color[1] = 1.0;
		RVdata[i].color[2] = 1.0;
		RVdata[i].color[3] = 1.0;

		RVdata[i].type = 0;
	}

	//todo iteratorで回すほうが早そう
	for (uint32_t i = 0; i < LineIndexData.size(); i++) {
		RIdata[i] = LineIndexData[i];
	}

	lva.resetvertarray(LineVertData.size(), RVdata, LineIndexData.size(), RIdata);
}

bool Is_CollideRigidObbBvhNode(const RigidObbBvhNode* const ROBNode0, const RigidObbBvhNode* const ROBNode1, const fvec3& r0cm, const fvec3& r1cm, const fmat3& R0, const fmat3& R1, const float extend = 1.0)
{
	//SAT test

	fvec3 V0[8];
	fvec3 V1[8];

	V0[0] = r0cm + R0 * (ROBNode0->center + extend * fvec3(ROBNode0->Lengthx / -2.0, ROBNode0->Lengthy / -2.0, ROBNode0->Lengthz / 2.0));
	V0[1] = r0cm + R0 * (ROBNode0->center + extend * fvec3(ROBNode0->Lengthx / 2.0, ROBNode0->Lengthy / -2.0, ROBNode0->Lengthz / 2.0));
	V0[2] = r0cm + R0 * (ROBNode0->center + extend * fvec3(ROBNode0->Lengthx / 2.0, ROBNode0->Lengthy / 2.0, ROBNode0->Lengthz / 2.0));
	V0[3] = r0cm + R0 * (ROBNode0->center + extend * fvec3(ROBNode0->Lengthx / -2.0, ROBNode0->Lengthy / 2.0, ROBNode0->Lengthz / 2.0));

	V0[4] = r0cm + R0 * (ROBNode0->center + extend * fvec3(ROBNode0->Lengthx / -2.0, ROBNode0->Lengthy / -2.0, ROBNode0->Lengthz / -2.0));
	V0[5] = r0cm + R0 * (ROBNode0->center + extend * fvec3(ROBNode0->Lengthx / 2.0, ROBNode0->Lengthy / -2.0, ROBNode0->Lengthz / -2.0));
	V0[6] = r0cm + R0 * (ROBNode0->center + extend * fvec3(ROBNode0->Lengthx / 2.0, ROBNode0->Lengthy / 2.0, ROBNode0->Lengthz / -2.0));
	V0[7] = r0cm + R0 * (ROBNode0->center + extend * fvec3(ROBNode0->Lengthx / -2.0, ROBNode0->Lengthy / 2.0, ROBNode0->Lengthz / -2.0));

	V1[0] = r1cm + R1 * (ROBNode1->center + extend * fvec3(ROBNode1->Lengthx / -2.0, ROBNode1->Lengthy / -2.0, ROBNode1->Lengthz / 2.0));
	V1[1] = r1cm + R1 * (ROBNode1->center + extend * fvec3(ROBNode1->Lengthx / 2.0, ROBNode1->Lengthy / -2.0, ROBNode1->Lengthz / 2.0));
	V1[2] = r1cm + R1 * (ROBNode1->center + extend * fvec3(ROBNode1->Lengthx / 2.0, ROBNode1->Lengthy / 2.0, ROBNode1->Lengthz / 2.0));
	V1[3] = r1cm + R1 * (ROBNode1->center + extend * fvec3(ROBNode1->Lengthx / -2.0, ROBNode1->Lengthy / 2.0, ROBNode1->Lengthz / 2.0));

	V1[4] = r1cm + R1 * (ROBNode1->center + extend * fvec3(ROBNode1->Lengthx / -2.0, ROBNode1->Lengthy / -2.0, ROBNode1->Lengthz / -2.0));
	V1[5] = r1cm + R1 * (ROBNode1->center + extend * fvec3(ROBNode1->Lengthx / 2.0, ROBNode1->Lengthy / -2.0, ROBNode1->Lengthz / -2.0));
	V1[6] = r1cm + R1 * (ROBNode1->center + extend * fvec3(ROBNode1->Lengthx / 2.0, ROBNode1->Lengthy / 2.0, ROBNode1->Lengthz / -2.0));
	V1[7] = r1cm + R1 * (ROBNode1->center + extend * fvec3(ROBNode1->Lengthx / -2.0, ROBNode1->Lengthy / 2.0, ROBNode1->Lengthz / -2.0));

	fvec3 normallist[15];

	normallist[0] = R0 * fvec3(1.0, 0.0, 0.0);
	normallist[1] = R0 * fvec3(0.0, 1.0, 0.0);
	normallist[2] = R0 * fvec3(0.0, 0.0, 1.0);

	normallist[3] = R1 * fvec3(1.0, 0.0, 0.0);
	normallist[4] = R1 * fvec3(0.0, 1.0, 0.0);
	normallist[5] = R1 * fvec3(0.0, 0.0, 1.0);

	normallist[6] = (R0 * fvec3(1.0, 0.0, 0.0)).cross(R1 * fvec3(1.0, 0.0, 0.0));
	normallist[7] = (R0 * fvec3(0.0, 1.0, 0.0)).cross(R1 * fvec3(1.0, 0.0, 0.0));
	normallist[8] = (R0 * fvec3(0.0, 0.0, 1.0)).cross(R1 * fvec3(1.0, 0.0, 0.0));

	normallist[9]  = (R0 * fvec3(1.0, 0.0, 0.0)).cross(R1 * fvec3(0.0, 1.0, 0.0));
	normallist[10] = (R0 * fvec3(0.0, 1.0, 0.0)).cross(R1 * fvec3(0.0, 1.0, 0.0));
	normallist[11] = (R0 * fvec3(0.0, 0.0, 1.0)).cross(R1 * fvec3(0.0, 1.0, 0.0));

	normallist[12] = (R0 * fvec3(1.0, 0.0, 0.0)).cross(R1 * fvec3(0.0, 0.0, 1.0));
	normallist[13] = (R0 * fvec3(0.0, 1.0, 0.0)).cross(R1 * fvec3(0.0, 0.0, 1.0));
	normallist[14] = (R0 * fvec3(0.0, 0.0, 1.0)).cross(R1 * fvec3(0.0, 0.0, 1.0));

	//for (uint32_t j = 0; j < 8; j++) {
	//	std::cout << V0[j] << std::endl;
	//}
	//std::cout << std::endl;
	//for (uint32_t j = 0; j < 8; j++) {
	//	std::cout << V1[j] << std::endl;
	//}
	//std::cout << std::endl;
	//std::cout << std::endl;

	for (uint32_t i = 0; i < 15; i++) {

		if (normallist[i].sqlength() > 0.1) {
			float max0 = V0[0].dot(normallist[i]);
			float min0 = V0[0].dot(normallist[i]);

			float max1 = V1[0].dot(normallist[i]);
			float min1 = V1[0].dot(normallist[i]);

			for (uint32_t j = 0; j < 8; j++) {
				if (V0[j].dot(normallist[i]) < min0) {
					min0 = V0[j].dot(normallist[i]);
				}
				if (V0[j].dot(normallist[i]) > max0) {
					max0 = V0[j].dot(normallist[i]);
				}
			}

			for (uint32_t j = 0; j < 8; j++) {
				if (V1[j].dot(normallist[i]) < min1) {
					min1 = V1[j].dot(normallist[i]);
				}
				if (V1[j].dot(normallist[i]) > max1) {
					max1 = V1[j].dot(normallist[i]);
				}
			}

			if (min0 > max1 || min1 > max0)
				return false;

			//if (min0 < max1 && max0 > min1) {
			//	for (uint32_t j = 0; j < 8; j++) {
			//		std::cout << V0[j].x << " ";
			//	}
			//	std::cout << std::endl;
			//	for (uint32_t j = 0; j < 8; j++) {
			//		std::cout << V1[j].x << " ";
			//	}
			//	std::cout << std::endl;

			//	std::cout << max0 << " " << min0 << " " << max1 << " " << min1 << std::endl;
			//	std::cout << normallist[i] << std::endl;
			//	return true;
			//}
			//if (min1 < max0 && max1 > min0) {
			//	std::cout << normallist[i] << std::endl;
			//	return true;
			//}
		}
	}

	return true;
}

uint32_t contactontriangle(const fvec3& a0, const fvec3& a1, const fvec3& a2, const fvec3& b0, const fvec3& b1, const fvec3& b2, ContactFeature& CF1, ContactFeature& CF2)
{

	ClosestDVV result[6];

	result[0] = DistTriangleSegment(a0, a1, a2, b0, b1);
	result[1] = DistTriangleSegment(a0, a1, a2, b1, b2);
	result[2] = DistTriangleSegment(a0, a1, a2, b2, b0);

	result[3] = DistTriangleSegment(b0, b1, b2, a0, a1);
	result[4] = DistTriangleSegment(b0, b1, b2, a1, a2);
	result[5] = DistTriangleSegment(b0, b1, b2, a2, a0);

	uint32_t flag = 0;

	for (uint32_t i = 3; i < 6; i++) {
		fvec3 hoge   = result[i].v0;
		result[i].v0 = result[i].v1;
		result[i].v1 = hoge;
	}

	fvec3 vcross;
	fvec3 normal0 = ((a1 - a0).cross(a2 - a0)).normalize();
	fvec3 normal1 = ((b1 - b0).cross(b2 - b0)).normalize();

	fvec3 v0, v1;
	fvec3 ax, ay;
	fvec3 bx, by;

	if (result[0].dist < Epsilon && result[1].dist < Epsilon) {
		flag   = 1;
		vcross = (result[0].v0 + result[1].v0) / 2.0;
	} else if (result[1].dist < Epsilon && result[2].dist < Epsilon) {
		flag   = 1;
		vcross = (result[1].v0 + result[2].v0) / 2.0;
	} else if (result[2].dist < Epsilon && result[0].dist < Epsilon) {
		flag   = 1;
		vcross = (result[2].v0 + result[0].v0) / 2.0;
	} else if (result[3].dist < Epsilon && result[4].dist < Epsilon) {
		flag   = 1;
		vcross = (result[3].v0 + result[4].v0) / 2.0;
	} else if (result[4].dist < Epsilon && result[5].dist < Epsilon) {
		flag   = 1;
		vcross = (result[4].v0 + result[5].v0) / 2.0;
	} else if (result[5].dist < Epsilon && result[3].dist < Epsilon) {
		flag   = 1;
		vcross = (result[5].v0 + result[3].v0) / 2.0;
	} else if (result[0].dist < Epsilon && result[3].dist < Epsilon) {
		flag   = 2;
		vcross = (result[0].v0 + result[3].v0) / 2.0;
		v0     = result[3].v0;
		v1     = result[0].v0;
		ax     = a0;
		ay     = a1;
		bx     = b0;
		by     = b1;
	} else if (result[0].dist < Epsilon && result[4].dist < Epsilon) {
		flag   = 2;
		vcross = (result[0].v0 + result[4].v0) / 2.0;
		v0     = result[4].v0;
		v1     = result[0].v0;
		ax     = a1;
		ay     = a2;
		bx     = b0;
		by     = b1;
	} else if (result[0].dist < Epsilon && result[5].dist < Epsilon) {
		flag   = 2;
		vcross = (result[0].v0 + result[5].v0) / 2.0;
		v0     = result[5].v0;
		v1     = result[0].v0;
		ax     = a2;
		ay     = a0;
		bx     = b0;
		by     = b1;
	} else if (result[1].dist < Epsilon && result[3].dist < Epsilon) {
		flag   = 2;
		vcross = (result[1].v0 + result[3].v0) / 2.0;
		v0     = result[3].v0;
		v1     = result[1].v0;
		ax     = a0;
		ay     = a1;
		bx     = b1;
		by     = b2;
	} else if (result[1].dist < Epsilon && result[4].dist < Epsilon) {
		flag   = 2;
		vcross = (result[1].v0 + result[4].v0) / 2.0;
		v0     = result[4].v0;
		v1     = result[1].v0;
		ax     = a1;
		ay     = a2;
		bx     = b1;
		by     = b2;
	} else if (result[1].dist < Epsilon && result[5].dist < Epsilon) {
		flag   = 2;
		vcross = (result[1].v0 + result[5].v0) / 2.0;
		v0     = result[5].v0;
		v1     = result[1].v0;
		ax     = a2;
		ay     = a0;
		bx     = b1;
		by     = b2;
	} else if (result[2].dist < Epsilon && result[3].dist < Epsilon) {
		flag   = 2;
		vcross = (result[2].v0 + result[3].v0) / 2.0;
		v0     = result[3].v0;
		v1     = result[2].v0;
		ax     = a0;
		ay     = a1;
		bx     = b2;
		by     = b0;
	} else if (result[2].dist < Epsilon && result[4].dist < Epsilon) {
		flag   = 2;
		vcross = (result[2].v0 + result[4].v0) / 2.0;
		v0     = result[4].v0;
		v1     = result[2].v0;
		ax     = a1;
		ay     = a2;
		bx     = b2;
		by     = b0;
	} else if (result[2].dist < Epsilon && result[5].dist < Epsilon) {
		flag   = 2;
		vcross = (result[2].v0 + result[5].v0) / 2.0;
		v0     = result[5].v0;
		v1     = result[2].v0;
		ax     = a2;
		ay     = a0;
		bx     = b2;
		by     = b0;
	}

	if (flag == 1) {
		fvec3 v01;
		if (b0.dot(normal0) < b1.dot(normal0) + Epsilon && b0.dot(normal0) < b2.dot(normal0) + Epsilon)
			v01 = b0;
		else if (b1.dot(normal0) < b0.dot(normal0) + Epsilon && b1.dot(normal0) < b2.dot(normal0) + Epsilon)
			v01 = b1;
		else
			v01 = b2;

		fvec3 v10;
		if (a0.dot(normal1) < a1.dot(normal1) + Epsilon && a0.dot(normal1) < a2.dot(normal1) + Epsilon)
			v10 = a0;
		else if (a1.dot(normal1) < a0.dot(normal1) + Epsilon && a1.dot(normal1) < a2.dot(normal1) + Epsilon)
			v10 = a1;
		else
			v10 = a2;

		//for (uint32_t i = 0; i < 3; i++) {
		//	fvec3 Ax, Ay;
		//	if (i == 0) {
		//		Ax = a0;
		//		Ay = a1;
		//	} else if (i == 1) {
		//		Ax = a1;
		//		Ay = a2;
		//	} else {
		//		Ax = a2;
		//		Ay = a0;
		//	}

		//	fvec3 na = ((Ay - Ax).cross(normal0)).normalize();
		//	if (std::abs((vcross - v01).dot(na)) > Epsilon) {
		//		float t = ((vcross - Ax).dot(na) / (vcross - v01).dot(na));
		//		if (0.0 < t && t < 1.0)
		//			v01 = t * (v01 - vcross) + vcross;
		//	}
		//}

		//for (uint32_t i = 0; i < 3; i++) {
		//	fvec3 Bx, By;
		//	if (i == 0) {
		//		Bx = b0;
		//		By = b1;
		//	} else if (i == 1) {
		//		Bx = b1;
		//		By = b2;
		//	} else {
		//		Bx = b2;
		//		By = b0;
		//	}

		//	fvec3 nb = ((By - Bx).cross(normal1)).normalize();
		//	if (std::abs((vcross - v10).dot(nb)) > Epsilon) {
		//		float s = ((vcross - Bx).dot(nb) / (vcross - v10).dot(nb));
		//		if (0.0 < s && s < 1.0)
		//			v10 = s * (v10 - vcross) + vcross;
		//	}
		//}

		//面同士が同じ方向を向いているときは無視してよい
		//このときに修正を加えるとさらにめり込む方向に押される可能性がある
		//manifoldなメッシュを考える場合は反対向きのメッシュも衝突している。
		if (normal0.dot(normal1) >= 0.0)
			return 0;

		if ((v01 - vcross).dot(-normal0) < (vcross - v10).dot(normal1)) {
			CF1 = ContactFeature { -normal0, vcross, v01 };
		} else {
			CF1 = ContactFeature { normal1, v10, vcross };
		}

		return 1;
	} else if (flag == 2) {
		fvec3 v01;
		if ((bx - vcross).dot(normal0) < (by - vcross).dot(normal0))
			v01 = bx;
		else
			v01 = by;

		fvec3 na = ((ay - ax).cross(normal0)).normalize();
		if (std::abs((v1 - v01).dot(na)) > Epsilon) {
			float t = ((v1 - ax).dot(na) / (v1 - v01).dot(na));
			if (0.0 < t && t < 1.0)
				v01 = t * (v01 - v1) + v1;
		}

		fvec3 v10;
		if ((ax - vcross).dot(normal1) < (ay - vcross).dot(normal1))
			v10 = ax;
		else
			v10 = ay;

		fvec3 nb = ((by - bx).cross(normal1)).normalize();
		if (std::abs((v0 - v10).dot(nb)) > Epsilon) {
			float s = ((v0 - bx).dot(nb) / (v0 - v10).dot(nb));
			if (0.0 < s && s < 1.0)
				v10 = s * (v10 - v0) + v0;
		}

		//std::cout << nb << std::endl;
		//std::cout << bx << std::endl;
		//std::cout << by << std::endl;
		//std::cout << v10 << std::endl;
		//std::cout << v0 << std::endl;
		//std::cout << s << std::endl;

		if ((v1 - v0).sqlength() < Epsilon)
			return 0;

		fvec3 normal = (v1 - v0).normalize();

		if ((v01 - vcross).dot(-normal0) < (vcross - v10).dot(normal1)) {
			if ((v01 - vcross).dot(-normal0) < (v1 - v0).dot(normal)) {
				CF1 = ContactFeature { -normal0, vcross, v01 };
			} else {
				CF1 = ContactFeature { normal, v0, v1 };
			}
		} else {
			if ((vcross - v10).dot(normal1) < (v1 - v0).dot(normal)) {
				CF1 = ContactFeature { normal1, v10, vcross };
			} else {
				CF1 = ContactFeature { normal, v0, v1 };
			}
		}

		return 1;
	}

	return 0;
}

struct PairROBNode {
	const RigidObbBvhNode* const ROBa;
	const RigidObbBvhNode* const ROBb;
};

void RigidObbBvhNodePotentialCollisionCheck(const RigidObbBvhNode* const ROBNode0, const RigidObbBvhNode* const ROBNode1, std::deque<TriangleInd>& PClist, std::deque<PairROBNode>& PROBList, const fvec3 r0cm, const fvec3 r1cm, const fmat3& R0, const fmat3& R1)
{
	//自己衝突は考えないため、ROBNode0 == ROBNode1 の場合は考えない

	if (ROBNode0->Type == 0 && ROBNode1->Type == 0) {

		//0RightChild vs 1RightChild
		if (Is_CollideRigidObbBvhNode(ROBNode0->RightChild, ROBNode1->RightChild, r0cm, r1cm, R0, R1))
			PROBList.emplace_back(PairROBNode { ROBNode0->RightChild, ROBNode1->RightChild });
		//0LeftChild vs 1RightChild
		if (Is_CollideRigidObbBvhNode(ROBNode0->LeftChild, ROBNode1->RightChild, r0cm, r1cm, R0, R1))
			PROBList.emplace_back(PairROBNode { ROBNode0->LeftChild, ROBNode1->RightChild });
		//0LeftChild vs 1RightChild
		if (Is_CollideRigidObbBvhNode(ROBNode0->RightChild, ROBNode1->LeftChild, r0cm, r1cm, R0, R1))
			PROBList.emplace_back(PairROBNode { ROBNode0->RightChild, ROBNode1->LeftChild });
		//0LeftChild vs 1LeftChild
		if (Is_CollideRigidObbBvhNode(ROBNode0->LeftChild, ROBNode1->LeftChild, r0cm, r1cm, R0, R1))
			PROBList.emplace_back(PairROBNode { ROBNode0->LeftChild, ROBNode1->LeftChild });

	} else if (ROBNode0->Type == 1 && ROBNode1->Type == 0) {

		//0 vs 1RightChild
		if (Is_CollideRigidObbBvhNode(ROBNode0, ROBNode1->RightChild, r0cm, r1cm, R0, R1))
			PROBList.emplace_back(PairROBNode { ROBNode0, ROBNode1->RightChild });
		//0 vs 1LeftChild
		if (Is_CollideRigidObbBvhNode(ROBNode0, ROBNode1->LeftChild, r0cm, r1cm, R0, R1))
			PROBList.emplace_back(PairROBNode { ROBNode0, ROBNode1->LeftChild });

	} else if (ROBNode0->Type == 0 && ROBNode1->Type == 1) {

		//0RightChild vs 1
		if (Is_CollideRigidObbBvhNode(ROBNode0->RightChild, ROBNode1, r0cm, r1cm, R0, R1))
			PROBList.emplace_back(PairROBNode { ROBNode0->RightChild, ROBNode1 });
		//0LeftChild vs 1
		if (Is_CollideRigidObbBvhNode(ROBNode0->LeftChild, ROBNode1, r0cm, r1cm, R0, R1))
			PROBList.emplace_back(PairROBNode { ROBNode0->LeftChild, ROBNode1 });

	} else {
		if (Is_CollideRigidObbBvhNode(ROBNode0, ROBNode1, r0cm, r1cm, R0, R1, 2.0)) {

			//ここでは荒く候補を決めるのでtriangle vs trianleはしない。
			PClist.emplace_back(
			    TriangleInd {
				ROBNode0->Index0,
				ROBNode0->Index1,
				ROBNode0->Index2,
				ROBNode1->Index0,
				ROBNode1->Index1,
				ROBNode1->Index2,
			    });
		}
	}

	if (!PROBList.empty()) {
		auto hoge = PROBList.back();
		PROBList.pop_back();
		RigidObbBvhNodePotentialCollisionCheck(hoge.ROBa, hoge.ROBb, PClist, PROBList, r0cm, r1cm, R0, R1);
	}
}

void RigidObbBvhPotentialCollisionPair(const RigidObbBvh& ROB0, const RigidObbBvh& ROB1, std::deque<TriangleInd>& PCList, const fvec3& r0cm, const fvec3& r1cm, const fquaternion& rotq0, const fquaternion& rotq1)
{
	std::deque<PairROBNode> PROBList;

	const fmat3 R0 = rotq0.qtomat();
	const fmat3 R1 = rotq1.qtomat();

	//ここで、ROB0,ROB1らに紐付けられたrigidbodyの座標、つまりexternal forceや前ステップのvelocityによる前進eulerでの更新をする前の座標で計算するとやばいことになる。

	if (Is_CollideRigidObbBvhNode(ROB0.RootNode, ROB1.RootNode, r0cm, r1cm, R0, R1)) {
		RigidObbBvhNodePotentialCollisionCheck(ROB0.RootNode, ROB1.RootNode, PCList, PROBList, r0cm, r1cm, R0, R1);
	}
}
