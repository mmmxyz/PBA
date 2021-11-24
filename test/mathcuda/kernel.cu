#include <cstdint>
#include <cuda_runtime.h>

#include "utils/mathfunc/mathfunc.hpp"
#include "utilscuda/mathfunc/mathfunc.cuh"

__global__ void detarray_kernel(const cudamat3* const array, float* const retptr, const uint32_t N)
{
	uint32_t i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i < N)
		retptr[i] = array[i].det();
}

void detarray(fmat3* const matptr, float* const retptr, const uint32_t& N)
{

	cudamat3* d_cmatarray;
	cudaMalloc(&d_cmatarray, N * sizeof(cudamat3));
	cudaDeviceSynchronize();

	float* d_retptr;
	cudaMalloc(&d_retptr, N * sizeof(float));
	cudaDeviceSynchronize();

	//cudamat3 cmatarray[N];
	//for (uint32_t i = 0; i < N; i++)
	//	for (uint32_t x = 0; x < 9; x++)
	//		cmatarray[i].m[x] = matptr[i].m[x];

	cudaMemcpy(d_cmatarray, matptr, N * sizeof(cudamat3), cudaMemcpyHostToDevice);
	cudaDeviceSynchronize();

	uint32_t Ds = N / 320;
	detarray_kernel<<<Ds + 1, 320>>>(d_cmatarray, d_retptr, N);
	cudaDeviceSynchronize();

	cudaMemcpy(retptr, d_retptr, N * sizeof(float), cudaMemcpyDeviceToHost);
	cudaDeviceSynchronize();

	cudaFree(d_cmatarray);
	cudaFree(d_retptr);
}
