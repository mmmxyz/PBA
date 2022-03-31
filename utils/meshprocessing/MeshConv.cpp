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
