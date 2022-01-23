
#include "utils/fileloader/TetGenLoader.hpp"

#include <iostream>
#include <fstream>
#include <cstdint>
#include <string>
#include <sstream>

bool LoadTetGentoTetrahedraMesh(
    const char* filename,
    fvec3** const vertdata,
    uint32_t& vertsize,
    uint32_t** const ilistdata,
    uint32_t& ilistsize,
    uint32_t** const facelistdata,
    uint32_t& facelistsize,
    uint32_t** const edgelistdata,
    uint32_t& edgelistsize,
    const fvec3& meshoffset,
    const float& meshscale)
{
	std::string filenamenode = filename;
	filenamenode += ".node";

	std::ifstream file(filenamenode);

	if (file.fail()) {
		std::cout << "fail to open file: " << filenamenode << std::endl;
		return false;
	}

	std::string line;

	{
		std::getline(file, line);
		//std::cout << line << std::endl;
		std::istringstream iss(line);
		std::string nsize;
		iss >> nsize;

		vertsize = std::stoi(nsize);
	}

	//std::cout << vertsize << std::endl;

	*vertdata = new fvec3[vertsize];

	for (uint32_t i = 0; i < vertsize; i++) {
		std::getline(file, line);

		std::istringstream iss(line);

		float hoge, x, y, z;
		iss >> hoge >> x >> y >> z;

		(*vertdata)[i].x = (x * meshscale) + meshoffset.x;
		(*vertdata)[i].y = (y * meshscale) + meshoffset.y;
		(*vertdata)[i].z = (z * meshscale) + meshoffset.z;

		//std::cout << x << " " << y << " " << z << std::endl;
	}

	file.close();

	std::string filenameele = filename;
	filenameele += ".ele";

	file.open(filenameele);

	if (file.fail()) {
		std::cout << "fail to open file: " << filenameele << std::endl;
		return false;
	}

	{
		std::getline(file, line);
		//std::cout << line << std::endl;
		std::istringstream iss(line);
		std::string nsize;
		iss >> nsize;

		ilistsize = 4 * std::stoi(nsize);
	}

	//std::cout << ilistsize << std::endl;

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

	file.close();

	std::string filenameface = filename;
	filenameface += ".face";

	file.open(filenameface);

	if (file.fail()) {
		std::cout << "fail to open file: " << filenameface << std::endl;
		return false;
	}

	{
		std::getline(file, line);
		//std::cout << line << std::endl;
		std::istringstream iss(line);
		std::string nsize;
		iss >> nsize;

		facelistsize = 3 * std::stoi(nsize);
	}

	//std::cout << facelistsize << std::endl;

	*facelistdata = new uint32_t[facelistsize];

	for (uint32_t i = 0; i < facelistsize / 3; i++) {
		std::getline(file, line);

		std::istringstream iss(line);

		uint32_t hoge, x, y, z, w;
		iss >> hoge >> x >> y >> z >> w;

		(*facelistdata)[3 * i + 0] = x;
		(*facelistdata)[3 * i + 1] = z;
		(*facelistdata)[3 * i + 2] = y;
	}

	file.close();

	std::string filenameedge = filename;
	filenameedge += ".edge";

	file.open(filenameedge);

	if (file.fail()) {
		std::cout << "fail to open file: " << filenameedge << std::endl;
		return false;
	}

	{
		std::getline(file, line);
		//std::cout << line << std::endl;
		std::istringstream iss(line);
		std::string nsize;
		iss >> nsize;

		edgelistsize = 2 * std::stoi(nsize);
	}

	//std::cout << edgelistsize << std::endl;

	*edgelistdata = new uint32_t[edgelistsize];

	for (uint32_t i = 0; i < edgelistsize / 2; i++) {
		std::getline(file, line);

		std::istringstream iss(line);

		uint32_t x, y, z, w;
		iss >> x >> y >> z >> w;

		(*edgelistdata)[2 * i + 0] = y;
		(*edgelistdata)[2 * i + 1] = z;
	}
}
