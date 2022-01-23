#include <iostream>
#include <vector>
#include <set>
#include <fstream>
#include <deque>
#include <sstream>
#include <string>
#include <algorithm>

#include "utils/mathfunc/mathfunc.hpp"

#include "opengl/drawobject.hpp"

#include "utils/fileloader/OBJLoader.hpp"

bool LoadOBJtoRenderTriangleMesh(
    const char* filename,
    trianglevertarray& varray,
    const fvec3& meshoffset,
    const float& meshscale)
{
	std::ifstream file(filename);
	if (file.fail()) {
		std::cout << "fail to open file: " << filename << std::endl;
		return false;
	}

	uint32_t FaceSize = 0;

	std::deque<fvec3> PositionSet;
	PositionSet.emplace_back(fvec3(0.0));
	std::deque<fvec3> NormalSet;
	NormalSet.emplace_back(fvec3(0.0));
	std::deque<fvec2> UVSet;
	UVSet.emplace_back(fvec2(0.0));

	std::string line;

	bool is_Normal = true;

	while (std::getline(file, line)) {
		std::istringstream iss(line);
		std::string linetype;
		iss >> linetype;

		if (linetype == "v") {
			float x, y, z;
			iss >> x >> y >> z;
			x *= meshscale;
			y *= meshscale;
			z *= meshscale;
			PositionSet.emplace_back(fvec3(x, y, z) + meshoffset);
		} else if (linetype == "vn") {
			float x, y, z;
			iss >> x >> y >> z;
			NormalSet.emplace_back(fvec3(x, y, z));
		} else if (linetype == "vt") {
			float u, v;
			iss >> u >> v;
			UVSet.emplace_back(fvec2(u, v));
		} else if (linetype == "f") {
			std::string stv[4];
			iss >> stv[0] >> stv[1] >> stv[2] >> stv[3];
			if (stv[3] == "")
				FaceSize += 1;
			else
				FaceSize += 2;
		}
	}

	file.clear();
	file.seekg(0L, std::ios::beg);

	std::deque<uint32_t> FaceIndexList(3 * FaceSize);
	std::deque<vec3<uint32_t>> VertIndexList(PositionSet.size() + 1);
	// x: index in PositionSet
	// y: index in UVSet
	// z: index in NormalSet

	for (auto& Vert : VertIndexList)
		Vert.x = 0;

	uint32_t CurrentFaceIndex = 0;
	while (std::getline(file, line)) {
		std::istringstream iss(line);
		std::string linetype;
		iss >> linetype;

		if (linetype == "f") {
			std::string stv[4];
			iss >> stv[0] >> stv[1] >> stv[2] >> stv[3];
			if (stv[3] == "") {
				vec3<uint32_t> VertIndexes[3];
				for (uint32_t i = 0; i < 3; i++) {
					VertIndexes[i].x = 0;
					VertIndexes[i].y = 0;
					VertIndexes[i].z = 0;

					uint32_t slashcounter = std::count(stv[i].begin(), stv[i].end(), '/');
					std::istringstream viss(stv[i]);
					std::string Index;

					std::getline(viss, Index, '/');
					VertIndexes[i].x = std::stoi(Index);

					std::getline(viss, Index, '/');
					if (Index != "" && slashcounter == 1) {
						// Position/UV

						VertIndexes[i].y = std::stoi(Index);

						is_Normal = false;
					} else if (Index != "" && slashcounter == 2) {
						// Position/UV/Normal

						VertIndexes[i].y = std::stoi(Index);

						std::getline(viss, Index, '/');
						VertIndexes[i].z = std::stoi(Index);
					} else if (Index == "" && slashcounter == 2) {
						// Position//Normal

						std::getline(viss, Index, '/');
						VertIndexes[i].z = std::stoi(Index);
					} else if ((viss.rdstate() & std::ios_base::eofbit) && slashcounter == 0) {
						// Position

						is_Normal = false;
					}
				}

				for (uint32_t i = 0, OrderArray[3] = { 0, 1, 2 }; i < 3; i++) {
					uint32_t j = OrderArray[i];
					if (VertIndexList[VertIndexes[j].x].x == 0
					    || (VertIndexList[VertIndexes[j].x].y == VertIndexes[j].y
						&& VertIndexList[VertIndexes[j].x].z == VertIndexes[j].z)) {

						VertIndexList[VertIndexes[j].x].x = VertIndexes[j].x;
						VertIndexList[VertIndexes[j].x].y = VertIndexes[j].y;
						VertIndexList[VertIndexes[j].x].z = VertIndexes[j].z;

						FaceIndexList[3 * CurrentFaceIndex + i] = VertIndexes[j].x;
					} else {
						VertIndexList.emplace_back(VertIndexes[j]);
						FaceIndexList[3 * CurrentFaceIndex + i] = VertIndexList.size() - 1;
					}
				}
				CurrentFaceIndex++;

			} else {

				vec3<uint32_t> VertIndexes[4];
				for (uint32_t i = 0; i < 4; i++) {
					VertIndexes[i].x = 0;
					VertIndexes[i].y = 0;
					VertIndexes[i].z = 0;

					uint32_t slashcounter = std::count(stv[i].begin(), stv[i].end(), '/');
					std::istringstream viss(stv[i]);
					std::string Index;

					std::getline(viss, Index, '/');
					VertIndexes[i].x = std::stoi(Index);

					std::getline(viss, Index, '/');
					if (Index != "" && slashcounter == 1) {
						// Position/UV

						VertIndexes[i].y = std::stoi(Index);

						is_Normal = false;
					} else if (Index != "" && slashcounter == 2) {
						// Position/UV/Normal

						VertIndexes[i].y = std::stoi(Index);

						std::getline(viss, Index, '/');
						VertIndexes[i].z = std::stoi(Index);
					} else if (Index == "" && slashcounter == 2) {
						// Position//Normal

						std::getline(viss, Index, '/');
						VertIndexes[i].z = std::stoi(Index);
					} else if ((viss.rdstate() & std::ios_base::eofbit) && slashcounter == 0) {
						// Position

						is_Normal = false;
					}
				}
				for (uint32_t i = 0, OrderArray[3] = { 0, 1, 2 }; i < 3; i++) {
					uint32_t j = OrderArray[i];
					if (VertIndexList[VertIndexes[j].x].x == 0
					    || (VertIndexList[VertIndexes[j].x].y == VertIndexes[j].y
						&& VertIndexList[VertIndexes[j].x].z == VertIndexes[j].z)) {

						VertIndexList[VertIndexes[j].x].x = VertIndexes[j].x;
						VertIndexList[VertIndexes[j].x].y = VertIndexes[j].y;
						VertIndexList[VertIndexes[j].x].z = VertIndexes[j].z;

						FaceIndexList[3 * CurrentFaceIndex + i] = VertIndexes[j].x;
					} else {
						VertIndexList.emplace_back(VertIndexes[j]);
						FaceIndexList[3 * CurrentFaceIndex + i] = VertIndexList.size() - 1;
					}
				}
				CurrentFaceIndex++;

				for (uint32_t i = 0, OrderArray[3] = { 0, 2, 3 }; i < 3; i++) {
					uint32_t j = OrderArray[i];
					if (VertIndexList[VertIndexes[j].x].x == 0
					    || (VertIndexList[VertIndexes[j].x].y == VertIndexes[j].y
						&& VertIndexList[VertIndexes[j].x].z == VertIndexes[j].z)) {

						VertIndexList[VertIndexes[j].x].x = VertIndexes[j].x;
						VertIndexList[VertIndexes[j].x].y = VertIndexes[j].y;
						VertIndexList[VertIndexes[j].x].z = VertIndexes[j].z;

						FaceIndexList[3 * CurrentFaceIndex + i] = VertIndexes[j].x;
					} else {
						VertIndexList.emplace_back(VertIndexes[j]);
						FaceIndexList[3 * CurrentFaceIndex + i] = VertIndexList.size() - 1;
					}
				}
				CurrentFaceIndex++;
			}
		}
	}

	if (!is_Normal) {
		NormalSet.clear();
		NormalSet.resize(VertIndexList.size());

		for (auto& normal : NormalSet)
			normal = fvec3(0.0);

		for (uint32_t i = 0; i < FaceSize; i++) {
			fvec3 v0 = PositionSet[VertIndexList[FaceIndexList[i * 3 + 0]].x];
			fvec3 v1 = PositionSet[VertIndexList[FaceIndexList[i * 3 + 1]].x];
			fvec3 v2 = PositionSet[VertIndexList[FaceIndexList[i * 3 + 2]].x];

			fvec3 normal = ((v1 - v0).cross(v2 - v0)).normalize();

			NormalSet[FaceIndexList[3 * i + 0]] = NormalSet[FaceIndexList[3 * i + 0]] + normal;
			NormalSet[FaceIndexList[3 * i + 1]] = NormalSet[FaceIndexList[3 * i + 1]] + normal;
			NormalSet[FaceIndexList[3 * i + 2]] = NormalSet[FaceIndexList[3 * i + 2]] + normal;
		}

		for (uint32_t i = 0; i < VertIndexList.size(); i++) {
			NormalSet[i]	   = NormalSet[i].normalize();
			VertIndexList[i].z = i;
		}
	}

	uint32_t* ilist = new uint32_t[FaceSize * 3];
	for (uint32_t i = 0; i < FaceSize * 3; i++)
		ilist[i] = FaceIndexList[i];

	varray.resetvertarray(VertIndexList.size(), FaceSize * 3, ilist);

	delete[] ilist;

	for (uint32_t i = 0; i < VertIndexList.size(); i++) {
		varray[i].position[0] = PositionSet[VertIndexList[i].x].x;
		varray[i].position[1] = PositionSet[VertIndexList[i].x].y;
		varray[i].position[2] = PositionSet[VertIndexList[i].x].z;

		varray[i].color[0] = 1.0f;
		varray[i].color[1] = 1.0f;
		varray[i].color[2] = 1.0f;
		varray[i].color[3] = 1.0f;

		varray[i].uv[0] = UVSet[VertIndexList[i].y].x;
		varray[i].uv[1] = UVSet[VertIndexList[i].y].y;

		fvec3 normal(
		    NormalSet[VertIndexList[i].z].x,
		    NormalSet[VertIndexList[i].z].y,
		    NormalSet[VertIndexList[i].z].z);
		normal = normal.normalize();

		varray[i].normal[0] = normal.x;
		varray[i].normal[1] = normal.y;
		varray[i].normal[2] = normal.z;

		varray[i].type = 0;
	}

	return true;
}

bool LoadOBJtoRenderEdgeMesh(
    const char* filename,
    linevertarray& varray,
    const fvec3& meshoffset,
    const float& meshscale)
{
	std::ifstream file(filename);
	if (file.fail()) {
		std::cout << "fail to open file: " << filename << std::endl;
		return false;
	}

	uint32_t FaceSize = 0;
	uint32_t VertSize = 0;

	std::deque<fvec3> PositionSet;
	PositionSet.emplace_back(fvec3(0.0));

	std::string line;

	while (std::getline(file, line)) {
		std::istringstream iss(line);
		std::string linetype;
		iss >> linetype;

		if (linetype == "v") {
			float x, y, z;
			iss >> x >> y >> z;
			x *= meshscale;
			y *= meshscale;
			z *= meshscale;
			PositionSet.emplace_back(fvec3(x, y, z) + meshoffset);
			VertSize++;
		} else if (linetype == "f") {
			std::string stv[4];
			iss >> stv[0] >> stv[1] >> stv[2] >> stv[3];
			if (stv[3] == "")
				FaceSize += 1;
			else
				FaceSize += 2;
		}
	}

	file.clear();
	file.seekg(0L, std::ios::beg);

	std::vector<std::set<uint32_t>> AdjacencyList(VertSize + 1);
	//AdjacencyList[0]は使わない。

	uint32_t CurrentFaceIndex = 0;
	while (std::getline(file, line)) {
		std::istringstream iss(line);
		std::string linetype;
		iss >> linetype;
		if (linetype == "f") {
			std::string stv[4];
			iss >> stv[0] >> stv[1] >> stv[2] >> stv[3];
			if (stv[3] == "") {
				uint32_t VertInFace[3];
				for (uint32_t i = 0; i < 3; i++) {
					std::istringstream viss(stv[i]);
					std::string Index;

					std::getline(viss, Index, '/');
					VertInFace[i] = std::stoi(Index);
				}

				for (uint32_t i = 0, EdgeStart[3] = { 0, 1, 2 }, EdgeEnd[3] = { 1, 2, 0 }; i < 3; i++) {
					uint32_t ES = VertInFace[EdgeStart[i]];
					uint32_t EE = VertInFace[EdgeEnd[i]];
					if (ES > EE) {
						uint32_t temp = EE;
						EE	      = ES;
						ES	      = temp;
					}

					AdjacencyList[ES].insert(EE);
				}

			} else {
				uint32_t VertInFace[4];
				for (uint32_t i = 0; i < 4; i++) {
					std::istringstream viss(stv[i]);
					std::string Index;

					std::getline(viss, Index, '/');
					VertInFace[i] = std::stoi(Index);
				}

				for (uint32_t i = 0, EdgeStart[3] = { 0, 1, 2 }, EdgeEnd[3] = { 1, 2, 0 }; i < 3; i++) {
					uint32_t ES = VertInFace[EdgeStart[i]];
					uint32_t EE = VertInFace[EdgeEnd[i]];
					if (ES > EE) {
						uint32_t temp = EE;
						EE	      = ES;
						ES	      = temp;
					}

					AdjacencyList[ES].insert(EE);
				}

				for (uint32_t i = 0, EdgeStart[3] = { 0, 2, 3 }, EdgeEnd[3] = { 2, 3, 0 }; i < 3; i++) {
					uint32_t ES = VertInFace[EdgeStart[i]];
					uint32_t EE = VertInFace[EdgeEnd[i]];
					if (ES > EE) {
						uint32_t temp = EE;
						EE	      = ES;
						ES	      = temp;
					}

					AdjacencyList[ES].insert(EE);
				}
			}
		}
	}

	uint32_t ilistsize = 0;
	for (auto x : AdjacencyList) {
		ilistsize += 2 * x.size();
	}

	uint32_t* ilistdata = new uint32_t[ilistsize];

	uint32_t counter = 0;
	for (uint32_t i = 1; i < AdjacencyList.size(); i++) {
		uint32_t counterinlist = 0;
		for (const auto& x : AdjacencyList[i]) {
			ilistdata[counter + 2 * counterinlist]	   = i;
			ilistdata[counter + 2 * counterinlist + 1] = x;
			counterinlist++;
		}
		counter += 2 * AdjacencyList[i].size();
	}

	varray.resetvertarray(VertSize + 1, ilistsize, ilistdata);

	delete[] ilistdata;

	for (uint32_t i = 0; i < VertSize + 1; i++) {
		varray[i].position[0] = PositionSet[i].x;
		varray[i].position[1] = PositionSet[i].y;
		varray[i].position[2] = PositionSet[i].z;

		varray[i].color[0] = 1.0f;
		varray[i].color[1] = 1.0f;
		varray[i].color[2] = 1.0f;
		varray[i].color[3] = 1.0f;

		varray[i].normal[0] = 0.0;
		varray[i].normal[1] = 0.0;
		varray[i].normal[2] = 0.0;

		varray[i].type = 0;
	}

	return true;
}

bool LoadOBJtoPhysicsTriangleMesh(
    const char* filename,
    fvec3** const vpdata,
    uint32_t& vertsize,
    uint32_t** const ilistdata,
    uint32_t& ilistsize,
    const fvec3& meshoffset,
    const float& meshscale)
{
	std::ifstream file(filename);
	if (file.fail()) {
		std::cout << "fail to open file: " << filename << std::endl;
		return false;
	}

	uint32_t FaceSize = 0;
	uint32_t VertSize = 0;

	std::deque<fvec3> PositionSet;

	std::string line;

	while (std::getline(file, line)) {
		std::istringstream iss(line);
		std::string linetype;
		iss >> linetype;

		if (linetype == "v") {
			float x, y, z;
			iss >> x >> y >> z;
			x *= meshscale;
			y *= meshscale;
			z *= meshscale;
			PositionSet.emplace_back(fvec3(x, y, z) + meshoffset);
			VertSize++;
		} else if (linetype == "f") {
			std::string stv[4];
			iss >> stv[0] >> stv[1] >> stv[2] >> stv[3];
			if (stv[3] == "")
				FaceSize += 1;
			else
				FaceSize += 2;
		}
	}

	file.clear();
	file.seekg(0L, std::ios::beg);

	(*vpdata) = new fvec3[VertSize];
	vertsize  = VertSize;
	for (uint32_t i = 0; i < VertSize; i++)
		(*vpdata)[i] = PositionSet[i];

	(*ilistdata) = new uint32_t[3 * FaceSize];
	ilistsize    = 3 * FaceSize;

	uint32_t CurrentFaceIndex = 0;
	while (std::getline(file, line)) {
		std::istringstream iss(line);
		std::string linetype;
		iss >> linetype;

		if (linetype == "f") {
			std::string stv[4];
			iss >> stv[0] >> stv[1] >> stv[2] >> stv[3];
			if (stv[3] == "") {
				uint32_t VertInFace[3];
				for (uint32_t i = 0; i < 3; i++) {
					std::istringstream viss(stv[i]);
					std::string Index;

					std::getline(viss, Index, '/');
					VertInFace[i] = std::stoi(Index);
				}

				for (uint32_t i = 0, OrderArray[3] = { 0, 1, 2 }; i < 3; i++) {
					(*ilistdata)[3 * CurrentFaceIndex + i] = VertInFace[OrderArray[i]] - 1;
				}
				CurrentFaceIndex++;
			} else {
				uint32_t VertInFace[4];
				for (uint32_t i = 0; i < 4; i++) {
					std::istringstream viss(stv[i]);
					std::string Index;

					std::getline(viss, Index, '/');
					VertInFace[i] = std::stoi(Index);
				}
				for (uint32_t i = 0, OrderArray[3] = { 0, 1, 2 }; i < 3; i++) {
					(*ilistdata)[3 * CurrentFaceIndex + i] = VertInFace[OrderArray[i]] - 1;
				}
				CurrentFaceIndex++;
				for (uint32_t i = 0, OrderArray[3] = { 0, 2, 3 }; i < 3; i++) {
					(*ilistdata)[3 * CurrentFaceIndex + i] = VertInFace[OrderArray[i]] - 1;
				}
				CurrentFaceIndex++;
			}
		}
	}

	return true;
}
