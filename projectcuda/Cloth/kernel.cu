
#include "utilscuda/mathfunc/mathfunc.cuh"

#include <cuda_runtime.h>

#include <iostream>

#include "kernel.hpp"

__constant__ fvec2* d_RestPositionList	  = nullptr;
__constant__ uint32_t* d_TriIndList	  = nullptr;
__constant__ uint32_t* d_InnerEdgeIndList = nullptr;
__constant__ uint32_t* d_EdgeList	  = nullptr;
__constant__ fvec4* d_InnerEdgeCList	  = nullptr;

__constant__ uint32_t* d_VtoTlist = nullptr;
__constant__ uint32_t* d_VtoTind  = nullptr;
__constant__ uint32_t* d_VtoElist = nullptr;
__constant__ uint32_t* d_VtoEind  = nullptr;
__constant__ uint32_t* d_VtoIlist = nullptr;
__constant__ uint32_t* d_VtoIind  = nullptr;

__constant__ fmat2* d_AList = nullptr;
__constant__ float* d_VList = nullptr;

__constant__ fvec3* d_dx    = nullptr;
__constant__ fvec3* d_tempp = nullptr;

__constant__ float* d_LambdaContinuumList = nullptr;
__constant__ float* d_LambdaBendList	  = nullptr;
__constant__ float* d_LambdaSpringList	  = nullptr;

__constant__ uint32_t d_vertsize;
__constant__ uint32_t d_tisize;
__constant__ uint32_t d_edgesize;
__constant__ uint32_t d_InnerEdgesize;

__constant__ float d_mass;
__constant__ float d_dt;

static fvec2* cpu_d_RestPositionList	= nullptr;
static uint32_t* cpu_d_TriIndList	= nullptr;
static uint32_t* cpu_d_InnerEdgeIndList = nullptr;
static uint32_t* cpu_d_EdgeList		= nullptr;
static fvec4* cpu_d_InnerEdgeCList	= nullptr;
static uint32_t* cpu_d_VtoTlist		= nullptr;
static uint32_t* cpu_d_VtoTind		= nullptr;
static uint32_t* cpu_d_VtoElist		= nullptr;
static uint32_t* cpu_d_VtoEind		= nullptr;
static uint32_t* cpu_d_VtoIlist		= nullptr;
static uint32_t* cpu_d_VtoIind		= nullptr;
static fmat2* cpu_d_AList		= nullptr;
static float* cpu_d_VList		= nullptr;
static fvec3* cpu_d_dx			= nullptr;
static fvec3* cpu_d_tempp		= nullptr;
static float* cpu_d_LambdaContinuumList = nullptr;
static float* cpu_d_LambdaBendList	= nullptr;
static float* cpu_d_LambdaSpringList	= nullptr;

static fvec2* cpu_RestPositionList    = nullptr;
static uint32_t* cpu_TriIndList	      = nullptr;
static uint32_t* cpu_InnerEdgeIndList = nullptr;
static uint32_t* cpu_EdgeList	      = nullptr;
static fvec4* cpu_InnerEdgeCList      = nullptr;
static uint32_t* cpu_VtoTlist	      = nullptr;
static uint32_t* cpu_VtoTind	      = nullptr;
static uint32_t* cpu_VtoElist	      = nullptr;
static uint32_t* cpu_VtoEind	      = nullptr;
static uint32_t* cpu_VtoIlist	      = nullptr;
static uint32_t* cpu_VtoIind	      = nullptr;
static fmat2* cpu_AList		      = nullptr;
static float* cpu_VList		      = nullptr;
static float* cpu_LambdaContinuumList = nullptr;
static float* cpu_LambdaBendList      = nullptr;
static float* cpu_LambdaSpringList    = nullptr;

static uint32_t cpu_vertsize;
static uint32_t cpu_tisize;
static uint32_t cpu_edgesize;
static uint32_t cpu_InnerEdgesize;

uint32_t dxsize;
uint32_t lambdasize;

__global__ void
FemElasticProjectGPU_Kernel(const float lambda, const float mu)
{

	uint32_t i = blockDim.x * blockIdx.x + threadIdx.x;

	if (i < d_tisize / 3) {
		fvec3 x0 = d_tempp[d_TriIndList[3 * i + 0]];
		fvec3 x1 = d_tempp[d_TriIndList[3 * i + 1]];
		fvec3 x2 = d_tempp[d_TriIndList[3 * i + 2]];

		fmat2 A	 = d_AList[i];
		fmat32 F = fmat32(x1 - x0, x2 - x0) * A;

		float C00 = F.m[2 * 0 + 0] * F.m[2 * 0 + 0] + F.m[2 * 1 + 0] * F.m[2 * 1 + 0] + F.m[2 * 2 + 0] * F.m[2 * 2 + 0];
		float C01 = F.m[2 * 0 + 0] * F.m[2 * 0 + 1] + F.m[2 * 1 + 0] * F.m[2 * 1 + 1] + F.m[2 * 2 + 0] * F.m[2 * 2 + 1];
		float C11 = F.m[2 * 0 + 1] * F.m[2 * 0 + 1] + F.m[2 * 1 + 1] * F.m[2 * 1 + 1] + F.m[2 * 2 + 1] * F.m[2 * 2 + 1];

		fmat2 E;
		E.m[0] = 0.5f * (C00 - 1);
		E.m[1] = 0.5f * C01;
		E.m[2] = 0.5f * C01;
		E.m[3] = 0.5f * (C11 - 1);

		float V = d_VList[i];
		V *= 0.01f;

		float W	 = V * (mu * E.sqlength() + 0.5f * lambda * E.trace() * E.trace());
		fmat32 B = V * (2.0f * mu * F * E + lambda * E.trace() * F);

		//printf("%f\n", C00);
		//printf("%f\n", C01);
		//printf("%f\n", C11);
		//if (i == 0) {
		//	printf("%f ", x0.x);
		//	printf("%f ", x0.y);
		//	printf("%f   ", x0.z);

		//	printf("%f ", X0.x);
		//	printf("%f \n", X0.y);
		//}

		if (W > 0.0) {
			float C = std::sqrt(2.0 * W);

			fmat32 BAt = B * A.transpose();

			fvec3 dC1 = (1.0f / C) * fvec3(BAt.m[0], BAt.m[2], BAt.m[4]);
			fvec3 dC2 = (1.0f / C) * fvec3(BAt.m[1], BAt.m[3], BAt.m[5]);
			fvec3 dC0 = -(dC1 + dC2);

			float dtdtdlambda = (-C - d_LambdaContinuumList[i]) / ((dC0.sqlength() + dC1.sqlength() + dC2.sqlength()) / d_mass + 1.0 / (d_dt * d_dt));
			dtdtdlambda *= 0.30;

			d_dx[3 * i + 0] = dtdtdlambda * (1.0f / d_mass) * dC0;
			d_dx[3 * i + 1] = dtdtdlambda * (1.0f / d_mass) * dC1;
			d_dx[3 * i + 2] = dtdtdlambda * (1.0f / d_mass) * dC2;

			d_LambdaContinuumList[i] += dtdtdlambda / (d_dt * d_dt);
		} else {

			d_dx[3 * i + 0] = fvec3(0.0);
			d_dx[3 * i + 1] = fvec3(0.0);
			d_dx[3 * i + 2] = fvec3(0.0);
		}
	}
}

__global__ void
FemBendProjectGPU_Kernel(const float bendCof)
{

	uint32_t i = blockDim.x * blockIdx.x + threadIdx.x;

	if (i < d_InnerEdgesize / 4) {

		fvec3 x0 = d_tempp[d_InnerEdgeIndList[4 * i + 0]];
		fvec3 x1 = d_tempp[d_InnerEdgeIndList[4 * i + 1]];
		fvec3 x2 = d_tempp[d_InnerEdgeIndList[4 * i + 2]];
		fvec3 x3 = d_tempp[d_InnerEdgeIndList[4 * i + 3]];

		fvec4 Cot = d_InnerEdgeCList[i];

		fmat4 X	   = fmat4(fvec4(x0), fvec4(x1), fvec4(x2), fvec4(x3));
		fvec4 XCot = X * Cot;

		float Q = bendCof * XCot.sqlength();

		if (Q > 0.0) {
			float C = std::sqrt(2.0 * Q);

			fvec3 dC0 = (1.0f / C) * bendCof * fvec3(XCot.x * Cot.x, XCot.y * Cot.x, XCot.z * Cot.x);
			fvec3 dC1 = (1.0f / C) * bendCof * fvec3(XCot.x * Cot.y, XCot.y * Cot.y, XCot.z * Cot.y);
			fvec3 dC2 = (1.0f / C) * bendCof * fvec3(XCot.x * Cot.z, XCot.y * Cot.z, XCot.z * Cot.z);
			fvec3 dC3 = (1.0f / C) * bendCof * fvec3(XCot.x * Cot.w, XCot.y * Cot.w, XCot.z * Cot.w);

			float dtdtdlambda = (-C - d_LambdaBendList[i]) / ((dC0.sqlength() + dC1.sqlength() + dC2.sqlength() + dC3.sqlength()) / d_mass + 1.0 / (d_dt * d_dt));
			dtdtdlambda *= 0.8;

			d_dx[4 * i + 0] = dtdtdlambda * (1.0f / d_mass) * dC0;
			d_dx[4 * i + 1] = dtdtdlambda * (1.0f / d_mass) * dC1;
			d_dx[4 * i + 2] = dtdtdlambda * (1.0f / d_mass) * dC2;
			d_dx[4 * i + 3] = dtdtdlambda * (1.0f / d_mass) * dC3;

			d_LambdaBendList[i] += dtdtdlambda / (d_dt * d_dt);
		} else {
			d_dx[4 * i + 0] = fvec3(0.0);
			d_dx[4 * i + 1] = fvec3(0.0);
			d_dx[4 * i + 2] = fvec3(0.0);
			d_dx[4 * i + 3] = fvec3(0.0);
		}
	}
}

__global__ void
FemAreaProjectGPU_Kernel()
{

	uint32_t i = blockDim.x * blockIdx.x + threadIdx.x;

	if (i < d_InnerEdgesize / 4) {
		fvec3 x0 = d_tempp[d_TriIndList[3 * i + 0]];
		fvec3 x1 = d_tempp[d_TriIndList[3 * i + 1]];
		fvec3 x2 = d_tempp[d_TriIndList[3 * i + 2]];

		float V = d_VList[i];

		float C = ((x1 - x0).cross(x2 - x0)).sqlength() - 4.0 * V * V;

		if (abs(C) > 0.001) {

			fvec3 dC1 = 2.0f * (x2 - x0).cross((x1 - x0).cross(x2 - x0));
			fvec3 dC2 = 2.0f * (x1 - x0).cross((x2 - x0).cross(x1 - x0));
			fvec3 dC0 = -(dC1 + dC2);

			float dtdtdlambda = (-C) / ((dC0.sqlength() + dC1.sqlength() + dC2.sqlength()) / d_mass);

			d_dx[3 * i + 0] = dtdtdlambda * (1.0f / d_mass) * dC0;
			d_dx[3 * i + 1] = dtdtdlambda * (1.0f / d_mass) * dC1;
			d_dx[3 * i + 2] = dtdtdlambda * (1.0f / d_mass) * dC2;

		} else {

			d_dx[3 * i + 0] = fvec3(0.0);
			d_dx[3 * i + 1] = fvec3(0.0);
			d_dx[3 * i + 2] = fvec3(0.0);
		}
	}
}

__global__ void
MassSpringProjectGPU_Kernel(const float SpringCof)
{

	uint32_t i = blockDim.x * blockIdx.x + threadIdx.x;

	if (i < d_edgesize / 2) {

		fvec2 X0 = d_RestPositionList[d_EdgeList[2 * i + 0]];
		fvec2 X1 = d_RestPositionList[d_EdgeList[2 * i + 1]];

		fvec3 x0 = d_tempp[d_EdgeList[2 * i + 0]];
		fvec3 x1 = d_tempp[d_EdgeList[2 * i + 1]];

		float C = SpringCof * ((x1 - x0).sqlength() - (X1 - X0).sqlength());

		if (abs(C) > 0.001) {

			fvec3 dC0 = SpringCof * (x0 - x1);
			fvec3 dC1 = SpringCof * (x1 - x0);

			float dtdtlambda = (-C - d_LambdaSpringList[i]) / ((dC0.sqlength() + dC1.sqlength()) / d_mass + 1.0 / (d_dt * d_dt));
			dtdtlambda *= 0.3;

			d_dx[2 * i + 0] = dtdtlambda * (1.0f / d_mass) * dC0;
			d_dx[2 * i + 1] = dtdtlambda * (1.0f / d_mass) * dC1;

			d_LambdaSpringList[i] += dtdtlambda / (d_dt * d_dt);
		} else {
			d_dx[2 * i + 0] = fvec3(0.0);
			d_dx[2 * i + 1] = fvec3(0.0);
		}
	}
}

__global__ void
dxtotemppTriangle()
{
	uint32_t i = blockDim.x * blockIdx.x + threadIdx.x;

	if (i < d_vertsize) {
		for (uint32_t j = d_VtoTind[i]; j < d_VtoTind[i + 1]; j++) {
			d_tempp[i] = d_tempp[i] + d_dx[d_VtoTlist[j]];
		}
	}
}

__global__ void
dxtotemppEdge()
{
	uint32_t i = blockDim.x * blockIdx.x + threadIdx.x;

	if (i < d_vertsize) {
		for (uint32_t j = d_VtoEind[i]; j < d_VtoEind[i + 1]; j++) {
			d_tempp[i] = d_tempp[i] + d_dx[d_VtoElist[j]];
		}
	}
}

__global__ void
dxtotemppInnerEdge()
{
	uint32_t i = blockDim.x * blockIdx.x + threadIdx.x;

	if (i < d_vertsize) {
		for (uint32_t j = d_VtoIind[i]; j < d_VtoIind[i + 1]; j++) {
			d_tempp[i] = d_tempp[i] + d_dx[d_VtoIlist[j]];
		}
	}
}

__global__ void
ClearLambda()
{
	uint32_t i = blockDim.x * blockIdx.x + threadIdx.x;

	if (i < d_tisize / 3) {
		d_LambdaContinuumList[i] = 0.0;
	}
	if (i < d_edgesize / 2) {
		d_LambdaSpringList[i] = 0.0;
	}
	if (i < d_InnerEdgesize / 4) {
		d_LambdaBendList[i] = 0.0;
	}
}

__global__ void
FixedProjectionGPU_Kernel(const float edgedist, const bool isfix, const uint32_t N)
{
	uint32_t i = blockDim.x * blockIdx.x + threadIdx.x;

	if (i < d_vertsize) {

		if (isfix) {
			if ((d_RestPositionList[i] - d_RestPositionList[N * (N - 1)]).length() < 0.5) {
				d_tempp[i].x = d_RestPositionList[i].x + 0.5 * edgedist;
				d_tempp[i].y = d_RestPositionList[i].y;
				d_tempp[i].z = 0.0f;
			}
			if ((d_RestPositionList[i] - d_RestPositionList[N * N - 1]).length() < 0.5) {
				d_tempp[i].x = d_RestPositionList[i].x - 0.5 * edgedist;
				d_tempp[i].y = d_RestPositionList[i].y;
				d_tempp[i].z = 0.0f;
			}

		} else {
			if (d_tempp[i].y < -12.0)
				d_tempp[i].y = -11.99;
		}
	}
}

__global__ void
PrintInfo()
{
	printf("vertsize %d \n", d_vertsize);
	printf("tisize %d \n", d_tisize);
	printf("edgesize %d \n", d_edgesize);
	printf("InnerEdgesize %d \n", d_InnerEdgesize);

	//for (int i = 0; i < d_InnerEdgesize / 4; i++) {
	//	printf("%f  ", d_InnerEdgeCList[i].x);
	//	printf("%f  ", d_InnerEdgeCList[i].y);
	//	printf("%f  ", d_InnerEdgeCList[i].z);
	//	printf("%f  ", d_InnerEdgeCList[i].w);
	//	printf("%d ", d_InnerEdgeIndList[4 * i + 0]);
	//	printf("%d ", d_InnerEdgeIndList[4 * i + 1]);
	//	printf("%d ", d_InnerEdgeIndList[4 * i + 2]);
	//	printf("%d \n", d_InnerEdgeIndList[4 * i + 3]);
	//}

	for (int i = 0; i < d_InnerEdgesize; i++) {
		printf("%d %d \n", d_VtoIlist[i], d_InnerEdgeIndList[d_VtoIlist[i]]);
	}
}

void Init(MeshInfo& mInfo)
{

	cpu_RestPositionList = mInfo.Restvertdata;
	cpu_TriIndList	     = mInfo.tilist;
	cpu_InnerEdgeIndList = mInfo.InnerEdgelist;
	cpu_EdgeList	     = mInfo.edgelist;
	cpu_InnerEdgeCList   = mInfo.InnerEdgeClist;

	cpu_VtoTlist = mInfo.VtoTlist;
	cpu_VtoTind  = mInfo.VtoTind;
	cpu_VtoElist = mInfo.VtoElist;
	cpu_VtoEind  = mInfo.VtoEind;
	cpu_VtoIlist = mInfo.VtoIlist;
	cpu_VtoIind  = mInfo.VtoIind;

	cpu_AList = mInfo.Alist;
	cpu_VList = mInfo.Vlist;

	cpu_vertsize	  = mInfo.vertsize;
	cpu_tisize	  = mInfo.tisize;
	cpu_edgesize	  = mInfo.edgesize;
	cpu_InnerEdgesize = mInfo.InnerEdgesize;

	dxsize	   = max(max(cpu_tisize, cpu_edgesize), cpu_InnerEdgesize);
	lambdasize = max(max(cpu_tisize / 3, cpu_edgesize / 2), cpu_InnerEdgesize / 4);

	//memory allocation
	cudaMalloc(&cpu_d_RestPositionList, cpu_vertsize * sizeof(fvec2));
	cudaMalloc(&cpu_d_TriIndList, cpu_tisize * sizeof(uint32_t));
	cudaMalloc(&cpu_d_InnerEdgeIndList, cpu_InnerEdgesize * sizeof(uint32_t));
	cudaMalloc(&cpu_d_InnerEdgeCList, (cpu_InnerEdgesize / 4) * sizeof(fvec4));
	cudaMalloc(&cpu_d_EdgeList, cpu_edgesize * sizeof(uint32_t));

	cudaMalloc(&cpu_d_VtoTlist, cpu_tisize * sizeof(uint32_t));
	cudaMalloc(&cpu_d_VtoElist, cpu_edgesize * sizeof(uint32_t));
	cudaMalloc(&cpu_d_VtoIlist, cpu_InnerEdgesize * sizeof(uint32_t));
	cudaMalloc(&cpu_d_VtoTind, (cpu_vertsize + 1) * sizeof(uint32_t));
	cudaMalloc(&cpu_d_VtoEind, (cpu_vertsize + 1) * sizeof(uint32_t));
	cudaMalloc(&cpu_d_VtoIind, (cpu_vertsize + 1) * sizeof(uint32_t));

	cudaMalloc(&cpu_d_AList, (cpu_tisize / 3) * sizeof(fmat2));
	cudaMalloc(&cpu_d_VList, (cpu_tisize / 3) * sizeof(float));
	cudaMalloc(&cpu_d_dx, dxsize * sizeof(fvec3));
	cudaMalloc(&cpu_d_tempp, cpu_vertsize * sizeof(fvec3));

	cudaMalloc(&cpu_d_LambdaContinuumList, (cpu_tisize / 3) * sizeof(float));
	cudaMalloc(&cpu_d_LambdaBendList, (cpu_InnerEdgesize / 4) * sizeof(float));
	cudaMalloc(&cpu_d_LambdaSpringList, (cpu_edgesize / 2) * sizeof(float));

	cudaDeviceSynchronize();

	//move pointer

	cudaMemcpyToSymbol(d_RestPositionList, &cpu_d_RestPositionList, sizeof(fvec2*));
	cudaMemcpyToSymbol(d_TriIndList, &cpu_d_TriIndList, sizeof(uint32_t*));
	cudaMemcpyToSymbol(d_InnerEdgeIndList, &cpu_d_InnerEdgeIndList, sizeof(uint32_t*));
	cudaMemcpyToSymbol(d_InnerEdgeCList, &cpu_d_InnerEdgeCList, sizeof(fvec4*));
	cudaMemcpyToSymbol(d_EdgeList, &cpu_d_EdgeList, sizeof(uint32_t*));

	cudaMemcpyToSymbol(d_VtoTlist, &cpu_d_VtoTlist, sizeof(uint32_t*));
	cudaMemcpyToSymbol(d_VtoElist, &cpu_d_VtoElist, sizeof(uint32_t*));
	cudaMemcpyToSymbol(d_VtoIlist, &cpu_d_VtoIlist, sizeof(uint32_t*));
	cudaMemcpyToSymbol(d_VtoTind, &cpu_d_VtoTind, sizeof(uint32_t*));
	cudaMemcpyToSymbol(d_VtoEind, &cpu_d_VtoEind, sizeof(uint32_t*));
	cudaMemcpyToSymbol(d_VtoIind, &cpu_d_VtoIind, sizeof(uint32_t*));

	cudaMemcpyToSymbol(d_AList, &cpu_d_AList, sizeof(fmat2*));
	cudaMemcpyToSymbol(d_VList, &cpu_d_VList, sizeof(float*));

	cudaMemcpyToSymbol(d_dx, &cpu_d_dx, sizeof(fvec3*));
	cudaMemcpyToSymbol(d_tempp, &cpu_d_tempp, sizeof(fvec3*));

	cudaMemcpyToSymbol(d_LambdaContinuumList, &cpu_d_LambdaContinuumList, sizeof(float*));
	cudaMemcpyToSymbol(d_LambdaBendList, &cpu_d_LambdaBendList, sizeof(float*));
	cudaMemcpyToSymbol(d_LambdaSpringList, &cpu_d_LambdaSpringList, sizeof(float*));

	cudaMemcpyToSymbol(d_vertsize, &cpu_vertsize, sizeof(uint32_t));
	cudaMemcpyToSymbol(d_tisize, &cpu_tisize, sizeof(uint32_t));
	cudaMemcpyToSymbol(d_edgesize, &cpu_edgesize, sizeof(uint32_t));
	cudaMemcpyToSymbol(d_InnerEdgesize, &cpu_InnerEdgesize, sizeof(uint32_t));

	cudaMemcpyToSymbol(d_mass, &mInfo.mass, sizeof(float));
	cudaMemcpyToSymbol(d_dt, &mInfo.dt, sizeof(float));

	//move static information
	cudaMemcpy(cpu_d_RestPositionList, cpu_RestPositionList, cpu_vertsize * sizeof(fvec2), cudaMemcpyHostToDevice);
	cudaMemcpy(cpu_d_TriIndList, cpu_TriIndList, cpu_tisize * sizeof(uint32_t), cudaMemcpyHostToDevice);
	cudaMemcpy(cpu_d_InnerEdgeIndList, cpu_InnerEdgeIndList, cpu_InnerEdgesize * sizeof(uint32_t), cudaMemcpyHostToDevice);
	cudaMemcpy(cpu_d_InnerEdgeCList, cpu_InnerEdgeCList, (cpu_InnerEdgesize / 4) * sizeof(fvec4), cudaMemcpyHostToDevice);
	cudaMemcpy(cpu_d_EdgeList, cpu_EdgeList, cpu_edgesize * sizeof(uint32_t), cudaMemcpyHostToDevice);

	cudaMemcpy(cpu_d_VtoTlist, cpu_VtoTlist, cpu_tisize * sizeof(uint32_t), cudaMemcpyHostToDevice);
	cudaMemcpy(cpu_d_VtoElist, cpu_VtoElist, cpu_edgesize * sizeof(uint32_t), cudaMemcpyHostToDevice);
	cudaMemcpy(cpu_d_VtoIlist, cpu_VtoIlist, cpu_InnerEdgesize * sizeof(uint32_t), cudaMemcpyHostToDevice);
	cudaMemcpy(cpu_d_VtoTind, cpu_VtoTind, (cpu_vertsize + 1) * sizeof(uint32_t), cudaMemcpyHostToDevice);
	cudaMemcpy(cpu_d_VtoEind, cpu_VtoEind, (cpu_vertsize + 1) * sizeof(uint32_t), cudaMemcpyHostToDevice);
	cudaMemcpy(cpu_d_VtoIind, cpu_VtoIind, (cpu_vertsize + 1) * sizeof(uint32_t), cudaMemcpyHostToDevice);

	cudaMemcpy(cpu_d_AList, cpu_AList, (cpu_tisize / 3) * sizeof(fmat2), cudaMemcpyHostToDevice);
	cudaMemcpy(cpu_d_VList, cpu_VList, (cpu_tisize / 3) * sizeof(float), cudaMemcpyHostToDevice);
	cudaDeviceSynchronize();

	//printf("vertsize %d \n", cpu_vertsize);
	//printf("tisize %d \n", cpu_tisize);
	//printf("edgesize %d \n", cpu_edgesize);
	//printf("InnerEdgesize %d \n", cpu_InnerEdgesize);

	//PrintInfo<<<1, 1>>>();

	//for (int i = 0; i < cpu_InnerEdgesize / 4; i++) {
	//	printf("%f  ", cpu_InnerEdgeCList[i].x);
	//	printf("%f  ", cpu_InnerEdgeCList[i].y);
	//	printf("%f  ", cpu_InnerEdgeCList[i].z);
	//	printf("%f  ", cpu_InnerEdgeCList[i].w);
	//	printf("%d ", cpu_InnerEdgeIndList[4 * i + 0]);
	//	printf("%d ", cpu_InnerEdgeIndList[4 * i + 1]);
	//	printf("%d ", cpu_InnerEdgeIndList[4 * i + 2]);
	//	printf("%d \n", cpu_InnerEdgeIndList[4 * i + 3]);
	//}
}

void ClearLambdaGPU()
{
	ClearLambda<<<lambdasize / 32 + 1, 32>>>();
}

void ElasticIterationGPU(fvec3* const tempp, const float lambda, const float mu, const float edgedist, const bool isfix, const uint32_t N, const float bendCof, const uint32_t iternum)
{

	cudaMemcpy(cpu_d_tempp, tempp, cpu_vertsize * sizeof(fvec3), cudaMemcpyHostToDevice);
	cudaDeviceSynchronize();

	uint32_t IEhsize = cpu_InnerEdgesize / 4;
	uint32_t IEDs	 = IEhsize / 320;

	uint32_t Thsize = cpu_tisize / 3;
	uint32_t TDs	= Thsize / 320;

	for (uint32_t x = 0; x < iternum; x++) {

		FixedProjectionGPU_Kernel<<<cpu_vertsize / 320 + 1, 320>>>(edgedist, isfix, N);
		cudaDeviceSynchronize();

		FemElasticProjectGPU_Kernel<<<TDs + 1, 320>>>(lambda, mu);
		cudaDeviceSynchronize();
		dxtotemppTriangle<<<cpu_vertsize / 320 + 1, 320>>>();
		cudaDeviceSynchronize();

		if (x % 10 == 0) {
			FemBendProjectGPU_Kernel<<<IEDs + 1, 320>>>(bendCof);
			cudaDeviceSynchronize();
			dxtotemppInnerEdge<<<cpu_vertsize / 320 + 1, 320>>>();
			cudaDeviceSynchronize();

			FemAreaProjectGPU_Kernel<<<TDs + 1, 320>>>();
			cudaDeviceSynchronize();
			dxtotemppTriangle<<<cpu_vertsize / 320 + 1, 320>>>();
			cudaDeviceSynchronize();
		}
	}

	cudaMemcpy(tempp, cpu_d_tempp, cpu_vertsize * sizeof(fvec3), cudaMemcpyDeviceToHost);
	cudaDeviceSynchronize();
}

void MassSpringIterationGPU(fvec3* const tempp, const float SpringCof, const float edgedist, const bool isfix, const uint32_t N, const float bendCof, const uint32_t iternum)
{

	cudaMemcpy(cpu_d_tempp, tempp, cpu_vertsize * sizeof(fvec3), cudaMemcpyHostToDevice);
	cudaDeviceSynchronize();

	uint32_t Ehsize = cpu_edgesize / 2;
	uint32_t EDS	= Ehsize / 640;

	uint32_t IEhsize = cpu_InnerEdgesize / 4;
	uint32_t IEDs	 = IEhsize / 640;

	uint32_t Thsize = cpu_tisize / 3;
	uint32_t TDs	= Thsize / 640;

	for (uint32_t x = 0; x < iternum; x++) {

		FixedProjectionGPU_Kernel<<<cpu_vertsize / 640 + 1, 640>>>(edgedist, isfix, N);
		cudaDeviceSynchronize();

		MassSpringProjectGPU_Kernel<<<EDS + 1, 640>>>(SpringCof);
		cudaDeviceSynchronize();
		dxtotemppEdge<<<cpu_vertsize / 640 + 1, 640>>>();
		cudaDeviceSynchronize();

		if (x % 10 == 0) {
			FemBendProjectGPU_Kernel<<<IEDs + 1, 640>>>(bendCof);
			cudaDeviceSynchronize();
			dxtotemppInnerEdge<<<cpu_vertsize / 640 + 1, 640>>>();
			cudaDeviceSynchronize();

			FemAreaProjectGPU_Kernel<<<TDs + 1, 640>>>();
			cudaDeviceSynchronize();
			dxtotemppTriangle<<<cpu_vertsize / 640 + 1, 640>>>();
			cudaDeviceSynchronize();
		}
	}

	cudaMemcpy(tempp, cpu_d_tempp, cpu_vertsize * sizeof(fvec3), cudaMemcpyDeviceToHost);
	cudaDeviceSynchronize();
}

void MassSpringProjectGPU(fvec3* const tempp, const float SpringCof)
{
	cudaMemcpy(cpu_d_tempp, tempp, cpu_vertsize * sizeof(fvec3), cudaMemcpyHostToDevice);
	cudaDeviceSynchronize();

	uint32_t Ehsize = cpu_edgesize / 2;
	uint32_t Ds	= Ehsize / 320;
	MassSpringProjectGPU_Kernel<<<Ds + 1, 320>>>(SpringCof);
	cudaDeviceSynchronize();

	dxtotemppEdge<<<cpu_vertsize / 320 + 1, 320>>>();
	cudaDeviceSynchronize();

	cudaMemcpy(tempp, cpu_d_tempp, cpu_vertsize * sizeof(fvec3), cudaMemcpyDeviceToHost);
	cudaDeviceSynchronize();
}

void FixedProjectionGPU(fvec3* const tempp, const float edgedist, const bool isfix, const uint32_t N)
{

	cudaMemcpy(cpu_d_tempp, tempp, cpu_vertsize * sizeof(fvec3), cudaMemcpyHostToDevice);
	cudaDeviceSynchronize();

	uint32_t Ds = cpu_vertsize / 320;
	FixedProjectionGPU_Kernel<<<Ds + 1, 320>>>(edgedist, isfix, N);
	cudaDeviceSynchronize();

	cudaMemcpy(tempp, cpu_d_tempp, cpu_vertsize * sizeof(fvec3), cudaMemcpyDeviceToHost);
	cudaDeviceSynchronize();
}

void FemBendProjectGPU(fvec3* const tempp, const float bendCof)
{

	cudaMemcpy(cpu_d_tempp, tempp, cpu_vertsize * sizeof(fvec3), cudaMemcpyHostToDevice);
	cudaDeviceSynchronize();

	uint32_t IEhsize = cpu_InnerEdgesize / 4;
	uint32_t Ds	 = IEhsize / 320;
	FemBendProjectGPU_Kernel<<<Ds + 1, 320>>>(bendCof);
	cudaDeviceSynchronize();

	dxtotemppInnerEdge<<<cpu_vertsize / 320 + 1, 320>>>();
	cudaDeviceSynchronize();

	cudaMemcpy(tempp, cpu_d_tempp, cpu_vertsize * sizeof(fvec3), cudaMemcpyDeviceToHost);
	cudaDeviceSynchronize();
}

void FemAreaProjectGPU(fvec3* const tempp)
{
	cudaMemcpy(cpu_d_tempp, tempp, cpu_vertsize * sizeof(fvec3), cudaMemcpyHostToDevice);
	cudaDeviceSynchronize();
	uint32_t Thsize = cpu_tisize / 3;
	uint32_t Ds	= Thsize / 320;
	FemAreaProjectGPU_Kernel<<<Ds + 1, 320>>>();
	cudaDeviceSynchronize();

	dxtotemppTriangle<<<cpu_vertsize / 320 + 1, 320>>>();
	cudaDeviceSynchronize();

	cudaMemcpy(tempp, cpu_d_tempp, cpu_vertsize * sizeof(fvec3), cudaMemcpyDeviceToHost);
	cudaDeviceSynchronize();
}

void FemElasticProjectGPU(fvec3* const tempp, const float lambda, const float mu)
{
	cudaMemcpy(cpu_d_tempp, tempp, cpu_vertsize * sizeof(fvec3), cudaMemcpyHostToDevice);
	cudaDeviceSynchronize();
	uint32_t Thsize = cpu_tisize / 3;
	uint32_t Ds	= Thsize / 320;
	FemElasticProjectGPU_Kernel<<<Ds + 1, 320>>>(lambda, mu);
	cudaDeviceSynchronize();

	dxtotemppTriangle<<<cpu_vertsize / 320 + 1, 320>>>();
	cudaDeviceSynchronize();

	cudaMemcpy(tempp, cpu_d_tempp, cpu_vertsize * sizeof(fvec3), cudaMemcpyDeviceToHost);
	cudaDeviceSynchronize();
}
