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
	uint32_t* VtoTlist;
	uint32_t* VtoTind;
	uint32_t* VtoElist;
	uint32_t* VtoEind;
	uint32_t* VtoIlist;
	uint32_t* VtoIind;
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
	    uint32_t* VtoTlist,
	    uint32_t* VtoTind,
	    uint32_t* VtoElist,
	    uint32_t* VtoEind,
	    uint32_t* VtoIlist,
	    uint32_t* VtoIind,
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
	    , VtoTlist(VtoTlist)
	    , VtoTind(VtoTind)
	    , VtoElist(VtoElist)
	    , VtoEind(VtoEind)
	    , VtoIlist(VtoIlist)
	    , VtoIind(VtoIind)
	    , Alist(Alist)
	    , Vlist(Vlist)
	    , dt(dt)
	    , mass(mass)
	{
	}
};

void Init(MeshInfo& minfo);

void MassSpringIterationGPU(fvec3* const tempp, const float SpringCof, const float edgedist, const bool isfix, const uint32_t N, const float bendCof, const uint32_t iternum);

void ElasticIterationGPU(fvec3* const tempp, const float lambda, const float mu, const float edgedist, const bool isfix, const uint32_t N, const float bendCof, const uint32_t iternum);

void MassSpringProjectGPU(fvec3* const tempp, const float SpringCof);

void FixedProjectionGPU(fvec3* const tempp, const float edgedist, const bool isfix, const uint32_t N);

void FemElasticProjectGPU(fvec3* const tempp, const float lambda, const float mu);

void FemBendProjectGPU(fvec3* const tempp, const float bendCof);

void FemAreaProjectGPU(fvec3* const tempp);

void ClearLambdaGPU();
