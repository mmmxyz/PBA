#pragma once

#include "utils/mathfunc/mathfunc.hpp"

struct MeshInfo {
	uint32_t vertsize;
	fvec2* Restvertdata;
	uint32_t* tilist;
	uint32_t tisize;
	uint32_t* edgelist;
	uint32_t edgesize;
	uint32_t* InnerEdgelist;
	fvec4* InnerEdgeClist;
	uint32_t InnerEdgesize;
	fmat2* Alist;
	float* Vlist;
	float mass;
	float dt;
	MeshInfo(
	    uint32_t vertsize,
	    fvec2* Restvertdata,
	    uint32_t* tilist,
	    uint32_t tisize,
	    uint32_t* edgelist,
	    uint32_t edgesize,
	    uint32_t* InnerEdgelist,
	    fvec4* InnerEdgeClist,
	    uint32_t InnerEdgesize,
	    fmat2* Alist,
	    float* Vlist,
	    float mass,
	    float dt)
	    : vertsize(vertsize)
	    , Restvertdata(Restvertdata)
	    , tilist(tilist)
	    , tisize(tisize)
	    , edgelist(edgelist)
	    , edgesize(edgesize)
	    , InnerEdgelist(InnerEdgelist)
	    , InnerEdgeClist(InnerEdgeClist)
	    , InnerEdgesize(InnerEdgesize)
	    , Alist(Alist)
	    , Vlist(Vlist)
	    , dt(dt)
	    , mass(mass)
	{
	}
};

//void Init(const fvec2* const RestPositionList, const uint32_t* const TriIndList, const uint32_t* const edgelist, const uint32_t* const InnerEdgeIndList, const fvec4* const InnerEdgeCList, const fmat2* const AList, const float* const VList, const uint32_t N);

void Init(MeshInfo& minfo);

void FemElasticProjectGPU(fvec3* const tempp, const float lambda, const float mu);

void FemBendProjectGPU(fvec3* const tempp, const float bendCof);

void FemAreaProjectGPU(fvec3* const tempp);

void ClearLambdaGPU();
