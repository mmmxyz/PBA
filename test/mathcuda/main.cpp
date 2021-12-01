#include <iomanip>
#include <iostream>
#include <cstdint>
#include <cmath>
#include <random>

#include "utils/mathfunc/mathfunc.hpp"
#include "utilscuda/mathfunc/mathfunc.cuh"
#include "kernel.hpp"

int main(int argc, char const* argv[])
{

	uint32_t N = 1000;

	std::random_device s;
	std::mt19937 e(s());
	std::normal_distribution<> d(1.0, 1.0);

	fmat3* A = new fmat3[N];
	fmat3* B = new fmat3[N];

	for (uint32_t i = 0; i < N; i++)
		for (uint32_t x = 0; x < 9; x++) {
			A[i].m[x] = d(e);
			B[i].m[x] = d(e);
		}

	fmat3* mulList = new fmat3[N];

	detarray(A, B, mulList, N);

	for (uint32_t i = 0; i < N; i++) {
		std::cout << i << std::endl
			  << mulList[i] << A[i] * B[i] << (mulList[i] - A[i] * B[i]).sqlength() << std::endl;

		if ((mulList[i] - A[i] * B[i]).sqlength() > 0.000001)
			std::cout << "error" << std::endl;
	}

	delete[] A;
	delete[] B;
	delete[] mulList;

	return 0;
}
