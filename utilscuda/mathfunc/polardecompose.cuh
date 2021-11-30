#pragma once

#include "utilscuda/mathfunc/mathfunc.cuh"

__host__ __device__ cudaquaternion ExtractRotation(const cudamat3& A, const uint32_t& numitr, const cudaquaternion initq)
{

	cudaquaternion q = initq;

	cudavec3 a0 = cudavec3(A.m[0], A.m[3], A.m[6]);
	cudavec3 a1 = cudavec3(A.m[1], A.m[4], A.m[7]);
	cudavec3 a2 = cudavec3(A.m[2], A.m[5], A.m[8]);

	for (uint32_t i = 0; i < numitr; i++) {

		cudamat3 R = q.qtomat();

		cudavec3 r0 = cudavec3(R.m[0], R.m[3], R.m[6]);
		cudavec3 r1 = cudavec3(R.m[1], R.m[4], R.m[7]);
		cudavec3 r2 = cudavec3(R.m[2], R.m[5], R.m[8]);

		cudavec3 omega = r0.cross(a0) + r1.cross(a1) + r2.cross(a2);
		omega	       = omega / (std::abs(a0.dot(r0) + a1.dot(r1) + a2.dot(r2)) + 0.000000001);

		q = cudaquaternion(omega) * q;
		q = q.normalize();
	}
	return q;
}
