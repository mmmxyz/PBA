#include <iomanip>
#include <iostream>
#include <cstdint>
#include <cmath>
#include <random>

#include "utils/mathfunc/mathfunc.hpp"
#include "kernel.hpp"

int main(int argc, char const* argv[])
{

	uint32_t N = 10000;

	std::random_device s;
	std::mt19937 e(s());
	std::normal_distribution<> d(1.0, 1.0);

	fmat3* A = new fmat3[N];

	for (uint32_t i = 0; i < N; i++)
		for (uint32_t x = 0; x < 9; x++)
			A[i].m[x] = d(e);

	float* detList = new float[N];

	detarray(A, detList, N);

	for (uint32_t i = 0; i < N; i++) {
		std::cout << std::setprecision(10) << detList[i] << " " << A[i].det() << std::endl;
		//if (std::abs(detList[i] - A[i].det()) > 0.00001)
		//	std::cout << "erroe" << std::endl;
	}

	delete[] A;
	delete[] detList;

	return 0;
}
