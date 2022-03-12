#include "utilscuda/mathfunc/polardecompose.cuh"

#include "utilscuda/mathfunc/mathfunc.cuh"

#include "utilscuda/fem/fem.cuh"

#include <cuda_runtime.h>
#include <iostream>

#include "kernel.hpp"

__constant__ fvec3* d_RestPositionList = nullptr;
__constant__ uint32_t* d_elementlist   = nullptr;

__constant__ uint32_t* d_VtoElist = nullptr;
__constant__ uint32_t* d_VtoEind  = nullptr;

__constant__ fmat3* d_AList = nullptr;
__constant__ float* d_VList = nullptr;

__constant__ fquaternion* d_qlist = nullptr;

__constant__ fvec3* d_dx    = nullptr;
__constant__ fvec3* d_tempp = nullptr;

__constant__ float* d_LambdaContinuumList = nullptr;

__constant__ uint32_t d_vertsize;
__constant__ uint32_t d_elementsize;

__constant__ float d_mass;
__constant__ float d_dt;

static fvec3* cpu_d_RestPositionList = nullptr;

static uint32_t* cpu_d_elementlist = nullptr;

static uint32_t* cpu_d_VtoElist = nullptr;
static uint32_t* cpu_d_VtoEind	= nullptr;

static fmat3* cpu_d_AList = nullptr;
static float* cpu_d_VList = nullptr;

static fquaternion* cpu_d_qlist = nullptr;

static fvec3* cpu_d_dx	  = nullptr;
static fvec3* cpu_d_tempp = nullptr;

static float* cpu_d_LambdaContinuumList = nullptr;

static fvec3* cpu_RestPositionList = nullptr;

static uint32_t* cpu_elementlist = nullptr;

static uint32_t* cpu_VtoElist = nullptr;
static uint32_t* cpu_VtoEind  = nullptr;

static fmat3* cpu_AList = nullptr;
static float* cpu_VList = nullptr;

static fquaternion* cpu_qlist = nullptr;

static uint32_t cpu_vertsize;
static uint32_t cpu_elementsize;

uint32_t dxsize;
uint32_t lambdasize;

__global__ void
FemElasticProjectGPU_Kernel(const float lambda, const float mu, const int32_t MaterialInd)
{

	uint32_t i = blockDim.x * blockIdx.x + threadIdx.x;

	if (i < d_elementsize / 4) {
		const fvec3& x0 = d_tempp[d_elementlist[4 * i + 0]];
		const fvec3& x1 = d_tempp[d_elementlist[4 * i + 1]];
		const fvec3& x2 = d_tempp[d_elementlist[4 * i + 2]];
		const fvec3& x3 = d_tempp[d_elementlist[4 * i + 3]];

		const fmat3& A = d_AList[i];
		const fmat3 F  = fmat3(x1 - x0, x2 - x0, x3 - x0) * A;
		const fmat3 E  = 0.5f * (F.transpose() * F - fmat3::indentity());

		const float& V = d_VList[i];

		float W;
		fvec3 dx0, dx1, dx2, dx3;
		bool ValidEnergy;

		if (MaterialInd == 0) {
			ValidEnergy = CUDAutils::FemElasticDxStVenant(F, E, A, V, lambda, mu, W, dx0, dx1, dx2, dx3);
		} else if (MaterialInd == 1) {
			ValidEnergy = CUDAutils::FemElasticDxNeoHookean(F, E, A, V, lambda, mu, W, dx0, dx1, dx2, dx3);
		} else if (MaterialInd == 2) {
			ValidEnergy = CUDAutils::FemElasticDxCoRotational(F, E, A, d_qlist[i], V, lambda, mu, W, dx0, dx1, dx2, dx3);
		}

		if (ValidEnergy && W > 0.0001) {

			float C = sqrt(2.0 * W);

			fvec3 dC1 = (1.0f / C) * dx1;
			fvec3 dC2 = (1.0f / C) * dx2;
			fvec3 dC3 = (1.0f / C) * dx3;
			fvec3 dC0 = (1.0f / C) * dx0;

			float dtdtdlambda = (-C - d_LambdaContinuumList[i]) / ((dC0.sqlength() + dC1.sqlength() + dC2.sqlength() + dC3.sqlength()) / d_mass + 1.0 / (d_dt * d_dt));
			dtdtdlambda *= 0.10;

			d_dx[4 * i + 0] = dtdtdlambda * (1.0f / d_mass) * dC0;
			d_dx[4 * i + 1] = dtdtdlambda * (1.0f / d_mass) * dC1;
			d_dx[4 * i + 2] = dtdtdlambda * (1.0f / d_mass) * dC2;
			d_dx[4 * i + 3] = dtdtdlambda * (1.0f / d_mass) * dC3;

			d_LambdaContinuumList[i] += dtdtdlambda / (d_dt * d_dt);
		} else {
			d_dx[4 * i + 0] = fvec3(0.0);
			d_dx[4 * i + 1] = fvec3(0.0);
			d_dx[4 * i + 2] = fvec3(0.0);
			d_dx[4 * i + 3] = fvec3(0.0);
		}
	}
}

__global__ void
FemFixInversionGPU_Kernel(const float lambda, const float mu, const int32_t MaterialInd)
{
	uint32_t i = blockDim.x * blockIdx.x + threadIdx.x;

	if (i < d_elementsize / 4) {
		const fvec3& x0 = d_tempp[d_elementlist[4 * i + 0]];
		const fvec3& x1 = d_tempp[d_elementlist[4 * i + 1]];
		const fvec3& x2 = d_tempp[d_elementlist[4 * i + 2]];
		const fvec3& x3 = d_tempp[d_elementlist[4 * i + 3]];

		const fmat3& A = d_AList[i];
		const fmat3 F  = fmat3(x1 - x0, x2 - x0, x3 - x0) * A;
		const fmat3 E  = 0.5f * (F.transpose() * F - fmat3::indentity());

		if (F.det() < 0.00001 && lambda > 0.00001 && mu > 0.00001 && MaterialInd != 2) {
			fquaternion q = CUDAutils::ExtractRotation(F, 3, d_qlist[i]);
			d_qlist[i]    = q;
			fmat3 R	      = q.qtomat();

			fmat3 S	  = R.transpose() * F;
			float trS = S.trace();

			float W = 0.5 * (F.sqlength() - 2.0 * trS + 3);
			fmat3 B = (F - R);

			float C = std::sqrt(2.0 * W);

			fmat3 BAt = B * A.transpose();

			fvec3 dC1 = fvec3(BAt.m[0], BAt.m[3], BAt.m[6]);
			fvec3 dC2 = fvec3(BAt.m[1], BAt.m[4], BAt.m[7]);
			fvec3 dC3 = fvec3(BAt.m[2], BAt.m[5], BAt.m[8]);
			fvec3 dC0 = -(dC1 + dC2 + dC3);

			float dtdtdlambda = -C / ((dC0.sqlength() + dC1.sqlength() + dC2.sqlength() + dC3.sqlength()));
			dtdtdlambda *= 0.35;

			d_dx[4 * i + 0] = dtdtdlambda * dC0;
			d_dx[4 * i + 1] = dtdtdlambda * dC1;
			d_dx[4 * i + 2] = dtdtdlambda * dC2;
			d_dx[4 * i + 3] = dtdtdlambda * dC3;
		} else {
			d_dx[4 * i + 0] = fvec3(0.0);
			d_dx[4 * i + 1] = fvec3(0.0);
			d_dx[4 * i + 2] = fvec3(0.0);
			d_dx[4 * i + 3] = fvec3(0.0);
		}
	}
}

__global__ void
dxtotemppElement()
{
	uint32_t i = blockDim.x * blockIdx.x + threadIdx.x;

	if (i < d_vertsize) {
		for (uint32_t j = d_VtoEind[i]; j < d_VtoEind[i + 1]; j++) {
			d_tempp[i] = d_tempp[i] + d_dx[d_VtoElist[j]];
		}
	}
}

__global__ void
ClearLambda()
{
	uint32_t i = blockDim.x * blockIdx.x + threadIdx.x;

	if (i < d_elementsize / 4) {
		d_LambdaContinuumList[i] = 0.0;
	}
}

__global__ void
ClearQ()
{
	uint32_t i = blockDim.x * blockIdx.x + threadIdx.x;

	if (i < d_elementsize / 4) {
		d_qlist[i] = fquaternion(0.0, 0.0, 0.0, 1.0);
	}
}

__global__ void
FixedProjectionGPU_Kernel()
{
	uint32_t i = blockDim.x * blockIdx.x + threadIdx.x;

	if (i < d_vertsize) {
		fvec3& position = d_tempp[i];
		if (position.y < -12.0)
			position.y = -11.99;
		if (position.y > 12.0)
			position.y = 11.99;

		if (position.x < -12.0)
			position.x = -11.99;
		if (position.x > 12.0)
			position.x = 11.99;

		if (position.z < -12.0)
			position.z = -11.99;
		if (position.z > 12.0)
			position.z = 11.99;
	}
}

__global__ void
PrintInfo()
{
	//printf("Esize %d \n", d_elementsize);

	//for (int i = 0; i < d_elementsize; i++) {
	//	printf("%d %d \n", d_VtoElist[i], d_elementlist[d_VtoElist[i]]);
	//}

	for (int i = 0; i < d_elementsize / 4; i++) {
		//printf("%f,%f,%f,%f \n", d_qlist[i].x, d_qlist[i].y, d_qlist[i].z, d_qlist[i].w);
		//printf("%f \n", d_AList[i].m[0]);
	}
}

void Init(MeshInfo& mInfo)
{

	cpu_RestPositionList = mInfo.Restvertdata;
	cpu_elementlist	     = mInfo.elementlist;

	cpu_VtoElist = mInfo.VtoElist;
	cpu_VtoEind  = mInfo.VtoEind;

	cpu_AList = mInfo.Alist;
	cpu_VList = mInfo.Vlist;

	cpu_qlist = mInfo.qlist;

	cpu_vertsize	= mInfo.vertsize;
	cpu_elementsize = mInfo.elementsize;

	dxsize	   = cpu_elementsize;
	lambdasize = cpu_elementsize / 4;

	//memory allocation
	cudaMalloc(&cpu_d_RestPositionList, cpu_vertsize * sizeof(fvec3));
	cudaMalloc(&cpu_d_elementlist, cpu_elementsize * sizeof(uint32_t));

	cudaMalloc(&cpu_d_VtoElist, cpu_elementsize * sizeof(uint32_t));
	cudaMalloc(&cpu_d_VtoEind, (cpu_vertsize + 1) * sizeof(uint32_t));

	cudaMalloc(&cpu_d_AList, (cpu_elementsize / 4) * sizeof(fmat3));
	cudaMalloc(&cpu_d_VList, (cpu_elementsize / 4) * sizeof(float));
	cudaMalloc(&cpu_d_qlist, (cpu_elementsize / 4) * sizeof(fquaternion));

	cudaMalloc(&cpu_d_dx, dxsize * sizeof(fvec3));
	cudaMalloc(&cpu_d_tempp, cpu_vertsize * sizeof(fvec3));

	cudaMalloc(&cpu_d_LambdaContinuumList, (cpu_elementsize / 4) * sizeof(float));

	cudaDeviceSynchronize();

	//move pointer

	cudaMemcpyToSymbol(d_RestPositionList, &cpu_d_RestPositionList, sizeof(fvec3*));
	cudaMemcpyToSymbol(d_elementlist, &cpu_d_elementlist, sizeof(uint32_t*));

	cudaMemcpyToSymbol(d_VtoElist, &cpu_d_VtoElist, sizeof(uint32_t*));
	cudaMemcpyToSymbol(d_VtoEind, &cpu_d_VtoEind, sizeof(uint32_t*));

	cudaMemcpyToSymbol(d_AList, &cpu_d_AList, sizeof(fmat3*));
	cudaMemcpyToSymbol(d_VList, &cpu_d_VList, sizeof(float*));

	cudaMemcpyToSymbol(d_qlist, &cpu_d_qlist, sizeof(fquaternion*));

	cudaMemcpyToSymbol(d_dx, &cpu_d_dx, sizeof(fvec3*));
	cudaMemcpyToSymbol(d_tempp, &cpu_d_tempp, sizeof(fvec3*));

	cudaMemcpyToSymbol(d_LambdaContinuumList, &cpu_d_LambdaContinuumList, sizeof(float*));

	cudaMemcpyToSymbol(d_vertsize, &cpu_vertsize, sizeof(uint32_t));
	cudaMemcpyToSymbol(d_elementsize, &cpu_elementsize, sizeof(uint32_t));

	cudaMemcpyToSymbol(d_mass, &mInfo.mass, sizeof(float));
	cudaMemcpyToSymbol(d_dt, &mInfo.dt, sizeof(float));

	//move static information
	cudaMemcpy(cpu_d_RestPositionList, cpu_RestPositionList, cpu_vertsize * sizeof(fvec3), cudaMemcpyHostToDevice);
	cudaMemcpy(cpu_d_elementlist, cpu_elementlist, cpu_elementsize * sizeof(uint32_t), cudaMemcpyHostToDevice);

	cudaMemcpy(cpu_d_VtoElist, cpu_VtoElist, cpu_elementsize * sizeof(uint32_t), cudaMemcpyHostToDevice);
	cudaMemcpy(cpu_d_VtoEind, cpu_VtoEind, (cpu_vertsize + 1) * sizeof(uint32_t), cudaMemcpyHostToDevice);

	cudaMemcpy(cpu_d_AList, cpu_AList, (cpu_elementsize / 4) * sizeof(fmat3), cudaMemcpyHostToDevice);
	cudaMemcpy(cpu_d_VList, cpu_VList, (cpu_elementsize / 4) * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(cpu_d_qlist, cpu_qlist, (cpu_elementsize / 4) * sizeof(fquaternion), cudaMemcpyHostToDevice);

	cudaDeviceSynchronize();

	//printf("vertsize %d \n", cpu_vertsize);
	//printf("tisize %d \n", cpu_tisize);
	//printf("edgesize %d \n", cpu_edgesize);
	//printf("elementsize %d \n", cpu_elementsize / 4);

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

void Clearqlist()
{
	ClearQ<<<(cpu_elementsize / 4) / 32 + 1, 32>>>();
}

void ClearLambdaGPU()
{
	ClearLambda<<<lambdasize / 32 + 1, 32>>>();
}

void ElasticIterationGPU(fvec3* const tempp, const float lambda, const float mu, const int32_t MaterialInd, const uint32_t iternum)
{

	cudaMemcpy(cpu_d_tempp, tempp, cpu_vertsize * sizeof(fvec3), cudaMemcpyHostToDevice);
	cudaDeviceSynchronize();

	uint32_t Ehsize = cpu_elementsize / 4;
	uint32_t EDs	= Ehsize / 320;

	for (uint32_t x = 0; x < iternum; x++) {

		FixedProjectionGPU_Kernel<<<cpu_vertsize / 320 + 1, 320>>>();
		cudaDeviceSynchronize();

		FemElasticProjectGPU_Kernel<<<EDs + 1, 320>>>(lambda, mu, MaterialInd);
		cudaDeviceSynchronize();
		dxtotemppElement<<<cpu_vertsize / 320 + 1, 320>>>();
		cudaDeviceSynchronize();

		if (x % 10 == 0 && MaterialInd != 2) {
			FemFixInversionGPU_Kernel<<<EDs + 1, 320>>>(lambda, mu, MaterialInd);
			cudaDeviceSynchronize();
			dxtotemppElement<<<cpu_vertsize / 320 + 1, 320>>>();
			cudaDeviceSynchronize();
		}
	}

	cudaMemcpy(tempp, cpu_d_tempp, cpu_vertsize * sizeof(fvec3), cudaMemcpyDeviceToHost);
	cudaDeviceSynchronize();
}
