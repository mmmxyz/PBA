
#include "utilscuda/mathfunc/mathfunc.cuh"

#include <cuda_runtime.h>

#include <iostream>

static fvec2* d_RestPositionList    = nullptr;
static uint32_t* d_TriIndList	    = nullptr;
static uint32_t* d_InnerEdgeIndList = nullptr;
static fvec4* d_InnerEdgeCList	    = nullptr;
static fmat2* d_AList		    = nullptr;
static float* d_VList		    = nullptr;

static fvec3* d_dx	    = nullptr;
static float* d_ELambdaList = nullptr;
static fvec3* d_tempp	    = nullptr;

static uint32_t* TriIndListptr = nullptr;

__global__ void
FemElasticProjectGPU_Kernel(fvec3* const tempp, float* const ELambdaList, fvec3* const dx, const uint32_t N, const fvec2* const RestPositionList, const uint32_t* const TriIndList, const uint32_t* const InnerEdgeIndList, const fvec4* const InnerEdgeCList, const fmat2* const AList, const float* const VList, const float lambda, const float mu)
{

	const float mass = 7.68e-05;
	const float dt	 = 1.0 / 60.0;

	uint32_t i = blockDim.x * blockIdx.x + threadIdx.x;

	if (i < 2 * (N - 1) * (N - 1)) {

		fvec2 X0 = RestPositionList[TriIndList[3 * i + 0]];
		fvec2 X1 = RestPositionList[TriIndList[3 * i + 1]];
		fvec2 X2 = RestPositionList[TriIndList[3 * i + 2]];

		fvec3 x0 = tempp[TriIndList[3 * i + 0]];
		fvec3 x1 = tempp[TriIndList[3 * i + 1]];
		fvec3 x2 = tempp[TriIndList[3 * i + 2]];

		fmat2 A	 = AList[i];
		fmat32 F = fmat32(x1 - x0, x2 - x0) * A;

		float C00 = F.m[2 * 0 + 0] * F.m[2 * 0 + 0] + F.m[2 * 1 + 0] * F.m[2 * 1 + 0] + F.m[2 * 2 + 0] * F.m[2 * 2 + 0];
		float C01 = F.m[2 * 0 + 0] * F.m[2 * 0 + 1] + F.m[2 * 1 + 0] * F.m[2 * 1 + 1] + F.m[2 * 2 + 0] * F.m[2 * 2 + 1];
		float C11 = F.m[2 * 0 + 1] * F.m[2 * 0 + 1] + F.m[2 * 1 + 1] * F.m[2 * 1 + 1] + F.m[2 * 2 + 1] * F.m[2 * 2 + 1];

		fmat2 E;
		E.m[0] = 0.5f * (C00 - 1);
		E.m[1] = 0.5f * C01;
		E.m[2] = 0.5f * C01;
		E.m[3] = 0.5f * (C11 - 1);

		float V = VList[i];
		V *= 0.01f;

		float W	 = V * (mu * E.sqlength() + 0.5f * lambda * E.trace() * E.trace());
		fmat32 B = V * (2.0f * mu * F * E + lambda * E.trace() * F);

		//printf("%f\n", E.sqlength());

		if (W > 0.0) {
			float C = std::sqrt(2.0 * W);

			fmat32 BAt = B * A.transpose();

			fvec3 dC1 = (1.0f / C) * fvec3(BAt.m[0], BAt.m[2], BAt.m[4]);
			fvec3 dC2 = (1.0f / C) * fvec3(BAt.m[1], BAt.m[3], BAt.m[5]);
			fvec3 dC0 = -(dC1 + dC2);

			float dtdtdlambda = (-C - ELambdaList[i]) / ((dC0.sqlength() + dC1.sqlength() + dC2.sqlength()) / mass + 1.0 / (dt * dt));
			dtdtdlambda /= 2.0f;

			dx[3 * i + 0] = dtdtdlambda * (1.0f / mass) * dC0;
			dx[3 * i + 1] = dtdtdlambda * (1.0f / mass) * dC1;
			dx[3 * i + 2] = dtdtdlambda * (1.0f / mass) * dC2;

			ELambdaList[i] += dtdtdlambda / (dt * dt);
		} else {

			dx[3 * i + 0] = fvec3(0.0);
			dx[3 * i + 1] = fvec3(0.0);
			dx[3 * i + 2] = fvec3(0.0);
		}
	}
}

void Init(const fvec2* const RestPositionList, uint32_t* const TriIndList, const uint32_t* const InnerEdgeIndList, const fvec4* const InnerEdgeCList, const fmat2* const AList, const float* const VList, const uint32_t N)
{

	TriIndListptr = TriIndList;

	cudaMalloc(&d_RestPositionList, N * N * sizeof(fvec2));
	cudaMalloc(&d_TriIndList, 3 * 2 * (N - 1) * (N - 1) * sizeof(uint32_t));
	cudaMalloc(&d_InnerEdgeIndList, 4 * ((N - 1) * (N - 1) + 2 * (N - 2) * (N - 2)) * sizeof(uint32_t));
	cudaMalloc(&d_InnerEdgeCList, ((N - 1) * (N - 1) + 2 * (N - 2) * (N - 2)) * sizeof(uint32_t));
	cudaMalloc(&d_AList, 2 * (N - 1) * (N - 1) * sizeof(fmat2));
	cudaMalloc(&d_VList, 2 * (N - 1) * (N - 1) * sizeof(float));

	cudaMalloc(&d_dx, 3 * 2 * (N - 1) * (N - 1) * sizeof(fvec3));
	cudaMalloc(&d_ELambdaList, 3 * 2 * (N - 1) * (N - 1) * sizeof(float));
	cudaMalloc(&d_tempp, N * N * sizeof(fvec3));
	cudaDeviceSynchronize();

	cudaMemcpy(d_RestPositionList, RestPositionList, N * N * sizeof(fvec2), cudaMemcpyHostToDevice);
	cudaMemcpy(d_TriIndList, TriIndList, 3 * 2 * (N - 1) * (N - 1) * sizeof(uint32_t), cudaMemcpyHostToDevice);
	cudaMemcpy(d_InnerEdgeIndList, InnerEdgeIndList, 4 * ((N - 1) * (N - 1) + 2 * (N - 2) * (N - 2)) * sizeof(uint32_t), cudaMemcpyHostToDevice);
	cudaMemcpy(d_InnerEdgeCList, InnerEdgeCList, ((N - 1) * (N - 1) + 2 * (N - 2) * (N - 2)) * sizeof(uint32_t), cudaMemcpyHostToDevice);
	cudaMemcpy(d_AList, AList, 2 * (N - 1) * (N - 1) * sizeof(fmat2), cudaMemcpyHostToDevice);
	cudaMemcpy(d_VList, VList, 2 * (N - 1) * (N - 1) * sizeof(float), cudaMemcpyHostToDevice);
	cudaDeviceSynchronize();
}

void FemElasticProjectGPU(fvec3* const tempp, float* const ELambdaList, const uint32_t N, const float lambda, const float mu)
{
	cudaMemcpy(d_ELambdaList, ELambdaList, 3 * 2 * (N - 1) * (N - 1) * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_tempp, tempp, N * N * sizeof(fvec3), cudaMemcpyHostToDevice);
	cudaDeviceSynchronize();

	uint32_t Thsize = 2 * (N - 1) * (N - 1);
	uint32_t Ds	= Thsize / 320;
	FemElasticProjectGPU_Kernel<<<Ds + 1, 320>>>(d_tempp, d_ELambdaList, d_dx, N, d_RestPositionList, d_TriIndList, d_InnerEdgeIndList, d_InnerEdgeCList, d_AList, d_VList, lambda, mu);
	cudaDeviceSynchronize();

	cudaMemcpy(ELambdaList, d_ELambdaList, 3 * 2 * (N - 1) * (N - 1) * sizeof(float), cudaMemcpyDeviceToHost);
	fvec3* dx = new fvec3[3 * 2 * (N - 1) * (N - 1)];
	cudaMemcpy(dx, d_dx, 3 * 2 * (N - 1) * (N - 1) * sizeof(fvec3), cudaMemcpyDeviceToHost);
	cudaDeviceSynchronize();

	for (uint32_t i = 0; i < 3 * 2 * (N - 1) * (N - 1); i++) {
		tempp[TriIndListptr[i]] = tempp[TriIndListptr[i]] + dx[i];
	}

	delete[] dx;
}
