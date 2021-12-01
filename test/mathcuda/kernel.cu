#include <cstdint>
#include <cuda_runtime.h>

#include "utilscuda/mathfunc/mathfunc.cuh"

//template <class T>
//T func(T x)
//{
//	return x + 1;
//}
//
//template <>
//__host__ __device__ float func(float x)
//{
//	printf("hello");
//	return x + 1;
//}

__global__ void detarray_kernel(const fmat3* const array0, const fmat3* const array1, fmat3* const retptr, const uint32_t N)
{
	uint32_t i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i < N)
		retptr[i] = array0[i] * array1[i];
}

void detarray(const fmat3* const mata, const fmat3* const matb, fmat3* const retptr, const uint32_t& N)
{

	//func(10.0f);

	fmat3* d_cmatarray0;
	cudaMalloc(&d_cmatarray0, N * sizeof(fmat3));
	cudaDeviceSynchronize();

	fmat3* d_cmatarray1;
	cudaMalloc(&d_cmatarray1, N * sizeof(fmat3));
	cudaDeviceSynchronize();

	fmat3* d_retptr;
	cudaMalloc(&d_retptr, N * sizeof(fmat3));
	cudaDeviceSynchronize();

	cudaMemcpy(d_cmatarray0, mata, N * sizeof(fmat3), cudaMemcpyHostToDevice);
	cudaDeviceSynchronize();
	cudaMemcpy(d_cmatarray1, matb, N * sizeof(fmat3), cudaMemcpyHostToDevice);
	cudaDeviceSynchronize();

	uint32_t Ds = N / 320;
	detarray_kernel<<<Ds + 1, 320>>>(d_cmatarray0, d_cmatarray1, d_retptr, N);
	cudaDeviceSynchronize();

	cudaMemcpy(retptr, d_retptr, N * sizeof(fmat3), cudaMemcpyDeviceToHost);
	cudaDeviceSynchronize();

	cudaFree(d_cmatarray0);
	cudaFree(d_cmatarray1);
	cudaFree(d_retptr);
}
