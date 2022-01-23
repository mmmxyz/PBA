#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include <cstdint>

#include "utils/mathfunc/mathfunc.hpp"
#include "utils/fileloader/VTKLoader.hpp"

bool LoadVTKtoTetrahedraMesh(
    const char* filename,
    fvec3** const vertdata,
    uint32_t& vertsize,
    uint32_t** const ilistdata,
    uint32_t& ilistsize,
    const fvec3& meshoffset,
    const float& meshscale)
{
	std::ifstream file(filename);
	if (file.fail()) {
		std::cout << "fail to open file: " << filename << std::endl;
		return false;
	}

	std::string line;

	std::getline(file, line);
	std::getline(file, line);
	std::getline(file, line);
	std::getline(file, line);

	{
		std::getline(file, line);
		//std::cout << line << std::endl;
		std::istringstream iss(line);
		std::string nsize;
		iss >> nsize;
		iss >> nsize;

		vertsize = std::stoi(nsize);
	}

	*vertdata = new fvec3[vertsize];

	for (uint32_t i = 0; i < vertsize; i++) {
		std::getline(file, line);
		std::istringstream iss(line);
		float x, y, z;

		iss >> x >> y >> z;

		(*vertdata)[i].x = (x * meshscale) + meshoffset.x;
		(*vertdata)[i].y = (y * meshscale) + meshoffset.y;
		(*vertdata)[i].z = (z * meshscale) + meshoffset.z;

		//std::cout << x << " " << y << " " << z << std::endl;
	}

	std::getline(file, line);

	{
		std::getline(file, line);
		//std::cout << line << std::endl;
		std::istringstream iss(line);
		std::string nsize;
		iss >> nsize;
		iss >> nsize;

		ilistsize = 4 * std::stoi(nsize);
	}

	*ilistdata = new uint32_t[ilistsize];

	for (uint32_t i = 0; i < ilistsize / 4; i++) {
		std::getline(file, line);
		std::istringstream iss(line);

		uint32_t hoge, x, y, z, w;

		iss >> hoge >> x >> y >> z >> w;

		(*ilistdata)[4 * i + 0] = x;
		(*ilistdata)[4 * i + 1] = y;
		(*ilistdata)[4 * i + 2] = z;
		(*ilistdata)[4 * i + 3] = w;
	}
}
