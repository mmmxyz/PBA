#include <iostream>
#include <cstdint>
#include <omp.h>

#include "utils/mathfunc/mathfunc.hpp"

int main(int argc, char const* argv[])
{

	fvec3 a(1.0);

	std::cout << a << std::endl;

	fmat3 b;

	std::cout << b * a << std::endl;

	return 0;
}
