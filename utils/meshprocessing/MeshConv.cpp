#include <cstdint>
#include <vector>
#include <set>

#include "opengl/drawobject.hpp"

void ConvertPTMtoREM(
    const fvec3* const Pvdata,
    const uint32_t& Pvsize,
    const uint32_t* const Pidata,
    const uint32_t& Pisize,
    linevertarray& varray)
{

	uint32_t Rvsize = Pvsize;

	std::vector<std::set<uint32_t>> AdjacencyList(Pvsize);

	for (uint32_t i = 0; i < Pisize / 3; i++) {
		uint32_t VI0 = Pidata[3 * i + 0];
		uint32_t VI1 = Pidata[3 * i + 1];
		uint32_t VI2 = Pidata[3 * i + 2];

		if (VI0 < VI1)
			AdjacencyList[VI0].insert(VI1);
		else
			AdjacencyList[VI1].insert(VI0);

		if (VI1 < VI2)
			AdjacencyList[VI1].insert(VI2);
		else
			AdjacencyList[VI2].insert(VI1);

		if (VI2 < VI0)
			AdjacencyList[VI2].insert(VI0);
		else
			AdjacencyList[VI0].insert(VI2);
	}

	uint32_t Risize = 0;
	for (auto x : AdjacencyList) {
		Risize += 2 * x.size();
	}

	uint32_t* Ridata = new uint32_t[Risize];

	uint32_t counter = 0;
	for (uint32_t i = 0; i < AdjacencyList.size(); i++) {
		uint32_t counterinlist = 0;
		for (const auto& x : AdjacencyList[i]) {
			Ridata[counter + 2 * counterinlist]	= i;
			Ridata[counter + 2 * counterinlist + 1] = x;
			counterinlist++;
		}
		counter += 2 * AdjacencyList[i].size();
	}

	varray.resetvertarray(Rvsize, Risize, Ridata);

	delete[] Ridata;

	for (uint32_t i = 0; i < Rvsize; i++) {
		varray[i].position[0] = Pvdata[i].x;
		varray[i].position[1] = Pvdata[i].y;
		varray[i].position[2] = Pvdata[i].z;

		varray[i].color[0] = 1.0;
		varray[i].color[1] = 1.0;
		varray[i].color[2] = 1.0;
		varray[i].color[3] = 1.0;

		varray[i].type = 0;
	}
}

void ConvertPTMtoPEM(
    const uint32_t& Pvsize,
    const uint32_t* const Pidata,
    const uint32_t& Pisize,
    uint32_t** const edgedata,
    uint32_t& edgesize)
{

	std::vector<std::set<uint32_t>> AdjacencyList(Pvsize);

	for (uint32_t i = 0; i < Pisize / 3; i++) {
		uint32_t VI0 = Pidata[3 * i + 0];
		uint32_t VI1 = Pidata[3 * i + 1];
		uint32_t VI2 = Pidata[3 * i + 2];

		if (VI0 < VI1)
			AdjacencyList[VI0].insert(VI1);
		else
			AdjacencyList[VI1].insert(VI0);

		if (VI1 < VI2)
			AdjacencyList[VI1].insert(VI2);
		else
			AdjacencyList[VI2].insert(VI1);

		if (VI2 < VI0)
			AdjacencyList[VI2].insert(VI0);
		else
			AdjacencyList[VI0].insert(VI2);
	}

	edgesize = 0;
	for (auto x : AdjacencyList) {
		edgesize += 2 * x.size();
	}

	(*edgedata) = new uint32_t[edgesize];

	uint32_t counter = 0;
	for (uint32_t i = 0; i < AdjacencyList.size(); i++) {
		uint32_t counterinlist = 0;
		for (const auto& x : AdjacencyList[i]) {
			(*edgedata)[counter + 2 * counterinlist]     = i;
			(*edgedata)[counter + 2 * counterinlist + 1] = x;
			counterinlist++;
		}
		counter += 2 * AdjacencyList[i].size();
	}
}

void ConvertEVtoVE(
    const uint32_t vertsize,
    const uint32_t* const elementlist,
    const uint32_t elementsize,
    uint32_t** const elsup_index,
    uint32_t** const elsup)
{
	(*elsup_index) = new uint32_t[vertsize + 1];
	(*elsup)       = new uint32_t[elementsize];

	for (uint32_t i = 0; i < vertsize + 1; i++) {
		(*elsup_index)[i] = 0;
	}

	for (uint32_t i = 0; i < elementsize; i++) {
		(*elsup_index)[elementlist[i]]++;
	}

	uint32_t counter = 0;
	for (uint32_t i = 0; i < vertsize + 1; i++) {
		counter += (*elsup_index)[i];
		(*elsup_index)[i] = counter - (*elsup_index)[i];
	}

	for (uint32_t i = 0; i < elementsize; i++) {
		(*elsup)[(*elsup_index)[elementlist[i]]] = i;
		(*elsup_index)[elementlist[i]]++;
	}

	for (uint32_t i = vertsize; i > 0; i--) {
		(*elsup_index)[i] = (*elsup_index)[i - 1];
	}
	(*elsup_index)[0] = 0;
}

bool Is_SameTriangle(
    const uint32_t& a0,
    const uint32_t& a1,
    const uint32_t& a2,
    const uint32_t& b0,
    const uint32_t& b1,
    const uint32_t& b2)
{
	if (
	    (a0 == b0 && a1 == b1 && a2 == b2)
	    || (a0 == b1 && a1 == b2 && a2 == b0)
	    || (a0 == b2 && a1 == b0 && a2 == b1))
		return true;

	return false;
}

void ExtractBoundaryTetrahedra(
    const uint32_t vertsize,
    const uint32_t* const elementlist,
    const uint32_t elementsize,
    const uint32_t* const VtoElist,
    const uint32_t* const VtoEind,
    const uint32_t* const trilist,
    const uint32_t trisize,
    uint32_t** const boundaryelementlist,
    uint32_t& boundaryelementsize)
{
	int32_t* is_boundary = new int32_t[elementsize / 4];

	for (uint32_t i = 0; i < elementsize / 4; i++) {
		is_boundary[i] = 0;
	}

	for (uint32_t i = 0; i < trisize; i++) {
		for (uint32_t j = VtoEind[trilist[i]]; j < VtoEind[trilist[i] + 1]; j++) {
			uint32_t elementind = VtoElist[j] / 4;

			uint32_t a0 = trilist[3 * (i / 3) + 0];
			uint32_t a1 = trilist[3 * (i / 3) + 1];
			uint32_t a2 = trilist[3 * (i / 3) + 2];

			uint32_t b0 = elementlist[4 * elementind + 0];
			uint32_t b1 = elementlist[4 * elementind + 1];
			uint32_t b2 = elementlist[4 * elementind + 2];
			uint32_t b3 = elementlist[4 * elementind + 3];

			//if (Is_SameTriangle(a0, a1, a2, b0, b2, b1) || Is_SameTriangle(a0, a1, a2, b0, b1, b3) || Is_SameTriangle(a0, a1, a2, b1, b2, b3) || Is_SameTriangle(a0, a1, a2, b2, b0, b3))
			is_boundary[elementind] = 1;
		}
	}

	boundaryelementsize = 0;
	for (uint32_t i = 0; i < elementsize / 4; i++) {
		if (is_boundary[i] != 0) {
			boundaryelementsize += 4;
		}
	}

	(*boundaryelementlist) = new uint32_t[boundaryelementsize];

	uint32_t ind = 0;
	for (uint32_t i = 0; i < elementsize / 4; i++) {
		if (is_boundary[i] != 0) {
			(*boundaryelementlist)[4 * ind + 0] = elementlist[4 * i + 0];
			(*boundaryelementlist)[4 * ind + 1] = elementlist[4 * i + 1];
			(*boundaryelementlist)[4 * ind + 2] = elementlist[4 * i + 2];
			(*boundaryelementlist)[4 * ind + 3] = elementlist[4 * i + 3];
			ind++;
		}
	}

	delete[] is_boundary;
}

void ExtractElementVertex(
    const uint32_t vertsize,
    const uint32_t* const trilist,
    const uint32_t trisize,
    const uint32_t* const VtoTlist,
    const uint32_t* const VtoTind,
    uint32_t** const TritoVertlist,
    uint32_t& TritoVertsize,
    uint32_t** const nonredTlist)
{
	bool* is_boundary = new bool[vertsize];

	for (uint32_t i = 0; i < vertsize; i++) {
		is_boundary[i] = false;
	}

	for (uint32_t i = 0; i < trisize; i++) {
		is_boundary[trilist[i]] = true;
	}

	TritoVertsize = 0;
	for (uint32_t i = 0; i < vertsize; i++) {
		if (is_boundary[i] == true)
			TritoVertsize++;
	}

	*TritoVertlist = new uint32_t[TritoVertsize];

	uint32_t counter = 0;
	for (uint32_t i = 0; i < vertsize; i++) {
		if (is_boundary[i] == true) {
			(*TritoVertlist)[counter] = i;
			counter++;
		}
	}

	*nonredTlist = new uint32_t[trisize];

	for (uint32_t i = 0; i < TritoVertsize; i++) {
		for (uint32_t j = VtoTind[(*TritoVertlist)[i]]; j < VtoTind[(*TritoVertlist)[i] + 1]; j++) {
			(*nonredTlist)[VtoTlist[j]] = i;
		}
	}

	delete[] is_boundary;
}
