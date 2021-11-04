#include <cstdint>
#include <vector>
#include <set>

#include "utils/meshprocessing/MeshConv.hpp"

void ConvertPTMtoREM(
    const fvec3* const Pvdata,
    const uint32_t& Pvsize,
    const uint32_t* const Pidata,
    const uint32_t& Pisize,
    vertex** const Rvdata,
    uint32_t& Rvsize,
    uint32_t** const Ridata,
    uint32_t& Risize)
{

	Rvsize	= Pvsize;
	*Rvdata = new vertex[Rvsize];

	for (uint32_t i = 0; i < Rvsize; i++) {
		(*Rvdata)[i].position[0] = Pvdata[i].x;
		(*Rvdata)[i].position[1] = Pvdata[i].y;
		(*Rvdata)[i].position[2] = Pvdata[i].z;

		(*Rvdata)[i].color[0] = 1.0;
		(*Rvdata)[i].color[1] = 1.0;
		(*Rvdata)[i].color[2] = 1.0;
		(*Rvdata)[i].color[3] = 1.0;

		(*Rvdata)[i].type = 0;
	}

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

	Risize = 0;
	for (auto x : AdjacencyList) {
		Risize += 2 * x.size();
	}

	*Ridata = new uint32_t[Risize];

	uint32_t counter = 0;
	for (uint32_t i = 0; i < AdjacencyList.size(); i++) {
		uint32_t counterinlist = 0;
		for (const auto& x : AdjacencyList[i]) {
			(*Ridata)[counter + 2 * counterinlist]	   = i;
			(*Ridata)[counter + 2 * counterinlist + 1] = x;
			counterinlist++;
		}
		counter += 2 * AdjacencyList[i].size();
	}
}
