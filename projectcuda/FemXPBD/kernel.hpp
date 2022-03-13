#pragma once

#include "utils/mathfunc/mathfunc.hpp"

struct MeshInfo {
	uint32_t vertsize;
	fvec3* Restvertdata;
	uint32_t* elementlist;
	uint32_t elementsize;
	uint32_t* VtoElist;
	uint32_t* VtoEind;
	fmat3* Alist;
	float* Vlist;
	fquaternion* qlist;
	float mass;
	float dt;
	MeshInfo(
	    uint32_t vertsize,
	    fvec3* Restvertdata,
	    uint32_t* elementlist,
	    uint32_t elementsize,
	    uint32_t* VtoElist,
	    uint32_t* VtoEind,
	    fmat3* Alist,
	    float* Vlist,
	    fquaternion* qlist,
	    float mass,
	    float dt)
	    : vertsize(vertsize)
	    , Restvertdata(Restvertdata)
	    , elementlist(elementlist)
	    , elementsize(elementsize)
	    , VtoElist(VtoElist)
	    , VtoEind(VtoEind)
	    , Alist(Alist)
	    , Vlist(Vlist)
	    , qlist(qlist)
	    , dt(dt)
	    , mass(mass)
	{
	}
};

void Init(MeshInfo& minfo);

void ElasticIterationGPU(fvec3* const tempp, const float lambda, const float mu, const int32_t MaterialInd, const uint32_t iternum, const bool isinv);

void Clearqlist();

void ClearLambdaGPU();
