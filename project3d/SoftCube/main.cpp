#include <iostream>
#include <cmath>
#include <cstdint>
#include <vector>
#include <deque>
#include <algorithm>
#include <random>
#include <omp.h>

#include "opengl/visualizer.hpp"
#include "opengl/vertarray.hpp"
#include "opengl/drawobject.hpp"
#include "opengl/renderer3d.hpp"

#include "utils/mathfunc/mathfunc.hpp"

#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>

constexpr uint32_t FPS = 60;
constexpr float dt     = 1.0 / FPS;
//glfwSwapInterval(1); なので60FPS以下になる。
//これはモニタやGPUに依るかもしれない。

float mu     = 40;
float lambda = 80;
float rho    = 0.1;

float rightedge;
float rightrot = 0.001 * 3.1415f;
bool rightwall = false;

class CubeMesh {
    public:
	std::vector<fvec3> PositionList;
	std::vector<fvec3> RestPositionList;
	std::vector<fvec3> VelocityList;
	std::vector<uint32_t> TetIndList;

	const fvec3 bias;

	linevertarray lva;
	vertex* lvdata;
	uint32_t lvsize;
	uint32_t* lilist;
	uint32_t lisize;

	trianglevertarray tva;
	vertex* tvdata;
	uint32_t tvsize;
	uint32_t* tilist;
	uint32_t tisize;

	std::vector<float> Lamdalist;

	std::vector<fmat3> AList;
	std::vector<float> VList;

	const uint32_t N;

	int32_t MaterialInd = 0;

	float mass;
	//mass of nodepoint

	CubeMesh(uint32_t N, float length, const fvec3& bias)
	    : N(N)
	    , lva()
	    , tva()
	    , bias(bias)
	{

		PositionList.reserve(N * N * N);
		RestPositionList.reserve(N * N * N);
		VelocityList.reserve(N * N * N);
		TetIndList.reserve(4 * 6 * (N - 1) * (N - 1) * (N - 1));

		AList.reserve(6 * (N - 1) * (N - 1) * (N - 1));
		VList.reserve(6 * (N - 1) * (N - 1) * (N - 1));

		Lamdalist.resize(6 * (N - 1) * (N - 1) * (N - 1));

		for (uint32_t y = 0; y < N; y++) {
			float vy = -0.5 * length + length * (y / float(N - 1));
			for (uint32_t z = 0; z < N; z++) {
				float vz = -0.5 * length + length * (z / float(N - 1));
				for (uint32_t x = 0; x < N; x++) {
					float vx = -0.5 * length + length * (x / float(N - 1));
					PositionList.emplace_back(fvec3(vx, vy, vz) + bias);
					RestPositionList.emplace_back(fvec3(vx, vy, vz) + bias);
					VelocityList.emplace_back(fvec3(0.0));
				}
			}
		}

		rightedge = 0.5 * length;

		mass = length * length * length * rho / (N * N * N);

		for (uint32_t y = 0; y < N - 1; y++) {
			for (uint32_t z = 0; z < N - 1; z++) {
				for (uint32_t x = 0; x < N - 1; x++) {
					uint32_t Ind1 = N * N * y + N * z + x;
					uint32_t Ind0 = Ind1 + 1;
					uint32_t Ind2 = Ind0 + N;
					uint32_t Ind3 = Ind1 + N;

					uint32_t Ind4 = N * N + Ind0;
					uint32_t Ind5 = N * N + Ind1;
					uint32_t Ind6 = N * N + Ind2;
					uint32_t Ind7 = N * N + Ind3;

					TetIndList.emplace_back(Ind0);
					TetIndList.emplace_back(Ind1);
					TetIndList.emplace_back(Ind2);
					TetIndList.emplace_back(Ind4);

					TetIndList.emplace_back(Ind6);
					TetIndList.emplace_back(Ind4);
					TetIndList.emplace_back(Ind2);
					TetIndList.emplace_back(Ind5);

					TetIndList.emplace_back(Ind4);
					TetIndList.emplace_back(Ind5);
					TetIndList.emplace_back(Ind1);
					TetIndList.emplace_back(Ind2);

					TetIndList.emplace_back(Ind1);
					TetIndList.emplace_back(Ind3);
					TetIndList.emplace_back(Ind2);
					TetIndList.emplace_back(Ind7);

					TetIndList.emplace_back(Ind5);
					TetIndList.emplace_back(Ind7);
					TetIndList.emplace_back(Ind1);
					TetIndList.emplace_back(Ind2);

					TetIndList.emplace_back(Ind5);
					TetIndList.emplace_back(Ind6);
					TetIndList.emplace_back(Ind7);
					TetIndList.emplace_back(Ind2);
				}
			}
		}

		for (uint32_t i = 0; i < 6 * (N - 1) * (N - 1) * (N - 1); i++) {
			fvec3 X0 = RestPositionList[TetIndList[4 * i + 0]];
			fvec3 X1 = RestPositionList[TetIndList[4 * i + 1]];
			fvec3 X2 = RestPositionList[TetIndList[4 * i + 2]];
			fvec3 X3 = RestPositionList[TetIndList[4 * i + 3]];

			AList.emplace_back(fmat3(X1 - X0, X2 - X0, X3 - X0).inverse());
			VList.emplace_back(fvec3::STP(X1 - X0, X2 - X0, X3 - X0) / 6.0);
		}

		tvsize = 6 * N * N;
		tvdata = new vertex[tvsize];
		tisize = 6 * 6 * (N - 1) * (N - 1);
		tilist = new uint32_t[tisize];

		//face counter [0,6(N-1)^2]

		for (uint32_t i = 0; i < 6; i++) {
			for (uint32_t x = 0; x < N - 1; x++) {
				for (uint32_t y = 0; y < N - 1; y++) {
					uint32_t Ind0						      = N * y + x;
					tilist[6 * (N - 1) * (N - 1) * i + 6 * ((N - 1) * y + x) + 0] = i * N * N + Ind0;
					tilist[6 * (N - 1) * (N - 1) * i + 6 * ((N - 1) * y + x) + 1] = i * N * N + Ind0 + 1;
					tilist[6 * (N - 1) * (N - 1) * i + 6 * ((N - 1) * y + x) + 2] = i * N * N + Ind0 + 0 + N;

					tilist[6 * (N - 1) * (N - 1) * i + 6 * ((N - 1) * y + x) + 3] = i * N * N + Ind0 + 0 + N;
					tilist[6 * (N - 1) * (N - 1) * i + 6 * ((N - 1) * y + x) + 4] = i * N * N + Ind0 + 1;
					tilist[6 * (N - 1) * (N - 1) * i + 6 * ((N - 1) * y + x) + 5] = i * N * N + Ind0 + 1 + N;
				}
			}
		}

		for (uint32_t i = 0; i < 6 * N * N; i++) {
			tvdata[i].color[0] = 1.0;
			tvdata[i].color[1] = 1.0;
			tvdata[i].color[2] = 1.0;
			tvdata[i].color[3] = 1.0;
			tvdata[i].type	   = 1;
		}

		lvsize = 6 * N * N;
		lvdata = new vertex[lvsize];
		lisize = 6 * (6 * (N - 1) + 4) * (N - 1);
		lilist = new uint32_t[lisize];

		for (uint32_t i = 0; i < 6; i++) {
			for (uint32_t y = 0; y < N - 1; y++) {
				for (uint32_t x = 0; x < N - 1; x++) {
					uint32_t Ind0								    = N * y + x;
					lilist[(6 * (N - 1) + 4) * (N - 1) * i + (6 * (N - 1) + 2) * y + 6 * x + 0] = i * N * N + Ind0;
					lilist[(6 * (N - 1) + 4) * (N - 1) * i + (6 * (N - 1) + 2) * y + 6 * x + 1] = i * N * N + Ind0 + 1;

					lilist[(6 * (N - 1) + 4) * (N - 1) * i + (6 * (N - 1) + 2) * y + 6 * x + 2] = i * N * N + Ind0 + 1;
					lilist[(6 * (N - 1) + 4) * (N - 1) * i + (6 * (N - 1) + 2) * y + 6 * x + 3] = i * N * N + Ind0 + 0 + N;

					lilist[(6 * (N - 1) + 4) * (N - 1) * i + (6 * (N - 1) + 2) * y + 6 * x + 4] = i * N * N + Ind0 + 0 + N;
					lilist[(6 * (N - 1) + 4) * (N - 1) * i + (6 * (N - 1) + 2) * y + 6 * x + 5] = i * N * N + Ind0;
				}
			}
		}

		for (uint32_t i = 0; i < 6 * N * N; i++) {
			lvdata[i].color[0] = 0.0;
			lvdata[i].color[1] = 0.0;
			lvdata[i].color[2] = 0.0;
			lvdata[i].color[3] = 1.0;
			lvdata[i].type	   = 0;
		}

		this->UpdataVertarray();

		tva.resetvertarray(tvsize, tvdata, tisize, tilist);
		lva.resetvertarray(lvsize, lvdata, lisize, lilist);
	}

	void UpdataVertarray()
	{

		//z+
#pragma omp parallel
		for (uint32_t y = 0; y < N; y++) {
			for (uint32_t x = 0; x < N; x++) {
				fvec3 v					  = PositionList[N * N * y + N * (N - 1) + x];
				tvdata[0 * N * N + N * y + x].position[0] = v.x;
				tvdata[0 * N * N + N * y + x].position[1] = v.y;
				tvdata[0 * N * N + N * y + x].position[2] = v.z;
				lvdata[0 * N * N + N * y + x].position[0] = v.x;
				lvdata[0 * N * N + N * y + x].position[1] = v.y;
				lvdata[0 * N * N + N * y + x].position[2] = v.z;
			}
		}
		//z-
#pragma omp parallel
		for (uint32_t y = 0; y < N; y++) {
			for (uint32_t x = 0; x < N; x++) {
				fvec3 v					  = PositionList[N * N * y + N - 1 - x];
				tvdata[1 * N * N + N * y + x].position[0] = v.x;
				tvdata[1 * N * N + N * y + x].position[1] = v.y;
				tvdata[1 * N * N + N * y + x].position[2] = v.z;
				lvdata[1 * N * N + N * y + x].position[0] = v.x;
				lvdata[1 * N * N + N * y + x].position[1] = v.y;
				lvdata[1 * N * N + N * y + x].position[2] = v.z;
			}
		}

		//x+
#pragma omp parallel
		for (uint32_t y = 0; y < N; y++) {
			for (uint32_t z = 0; z < N; z++) {
				fvec3 v					  = PositionList[N * N * y + N * N - 1 - N * z];
				tvdata[2 * N * N + N * y + z].position[0] = v.x;
				tvdata[2 * N * N + N * y + z].position[1] = v.y;
				tvdata[2 * N * N + N * y + z].position[2] = v.z;
				lvdata[2 * N * N + N * y + z].position[0] = v.x;
				lvdata[2 * N * N + N * y + z].position[1] = v.y;
				lvdata[2 * N * N + N * y + z].position[2] = v.z;
			}
		}

		//x-
#pragma omp parallel
		for (uint32_t y = 0; y < N; y++) {
			for (uint32_t z = 0; z < N; z++) {
				fvec3 v					  = PositionList[N * N * y + N * z];
				tvdata[3 * N * N + N * y + z].position[0] = v.x;
				tvdata[3 * N * N + N * y + z].position[1] = v.y;
				tvdata[3 * N * N + N * y + z].position[2] = v.z;
				lvdata[3 * N * N + N * y + z].position[0] = v.x;
				lvdata[3 * N * N + N * y + z].position[1] = v.y;
				lvdata[3 * N * N + N * y + z].position[2] = v.z;
			}
		}

		//y+
#pragma omp parallel
		for (uint32_t z = 0; z < N; z++) {
			for (uint32_t x = 0; x < N; x++) {
				fvec3 v					  = PositionList[N * N * (N - 1) + N * z + (N - 1) - x];
				tvdata[4 * N * N + N * z + x].position[0] = v.x;
				tvdata[4 * N * N + N * z + x].position[1] = v.y;
				tvdata[4 * N * N + N * z + x].position[2] = v.z;
				lvdata[4 * N * N + N * z + x].position[0] = v.x;
				lvdata[4 * N * N + N * z + x].position[1] = v.y;
				lvdata[4 * N * N + N * z + x].position[2] = v.z;
			}
		}

		//y-
#pragma omp parallel
		for (uint32_t z = 0; z < N; z++) {
			for (uint32_t x = 0; x < N; x++) {
				fvec3 v					  = PositionList[N * N * 0 + N * z + x];
				tvdata[5 * N * N + N * z + x].position[0] = v.x;
				tvdata[5 * N * N + N * z + x].position[1] = v.y;
				tvdata[5 * N * N + N * z + x].position[2] = v.z;
				lvdata[5 * N * N + N * z + x].position[0] = v.x;
				lvdata[5 * N * N + N * z + x].position[1] = v.y;
				lvdata[5 * N * N + N * z + x].position[2] = v.z;
			}
		}

		//normal

		std::vector<fvec3> NormalSet(6 * N * N);

#pragma omp parallel
		for (auto& x : NormalSet)
			x = fvec3(0.0);

#pragma omp parallel
		for (uint32_t i = 0; i < 12 * (N - 1) * (N - 1); i++) {
			fvec3 v0 = fvec3(tvdata[tilist[3 * i + 0]].position);
			fvec3 v1 = fvec3(tvdata[tilist[3 * i + 1]].position);
			fvec3 v2 = fvec3(tvdata[tilist[3 * i + 2]].position);

			fvec3 normal = (v1 - v0).cross(v2 - v0);
			if (normal.sqlength() > 0.000001)
				normal = normal.normalize();

			NormalSet[tilist[3 * i + 0]] = NormalSet[tilist[3 * i + 0]] + normal;
			NormalSet[tilist[3 * i + 1]] = NormalSet[tilist[3 * i + 1]] + normal;
			NormalSet[tilist[3 * i + 2]] = NormalSet[tilist[3 * i + 2]] + normal;
		}

#pragma omp parallel
		for (uint32_t i = 0; i < 6 * N * N; i++) {
			tvdata[i].normal[0] = NormalSet[i].x;
			tvdata[i].normal[1] = NormalSet[i].y;
			tvdata[i].normal[2] = NormalSet[i].z;
		}
	}

	void Setdata()
	{
		tva.vboupdate();
		lva.vboupdate();
	}

	void ClearLamda()
	{
#pragma omp parallel
		for (auto& x : Lamdalist) {
			x = 0.0;
		}
	}

	void FemProjectGS(std::vector<fvec3>& TempPosition)
	{
		for (uint32_t i = 0; i < 6 * (N - 1) * (N - 1) * (N - 1); i++) {
			fvec3 X0 = RestPositionList[TetIndList[4 * i + 0]];
			fvec3 X1 = RestPositionList[TetIndList[4 * i + 1]];
			fvec3 X2 = RestPositionList[TetIndList[4 * i + 2]];
			fvec3 X3 = RestPositionList[TetIndList[4 * i + 3]];

			fvec3 x0 = TempPosition[TetIndList[4 * i + 0]];
			fvec3 x1 = TempPosition[TetIndList[4 * i + 1]];
			fvec3 x2 = TempPosition[TetIndList[4 * i + 2]];
			fvec3 x3 = TempPosition[TetIndList[4 * i + 3]];

			fmat3 A = AList[i];
			fmat3 F = mat3(x1 - x0, x2 - x0, x3 - x0) * A;
			fmat3 E = 0.5 * (F.transpose() * F - fmat3::indentity());

			float V = VList[i];

			float W;
			fmat3 B;

			if (MaterialInd == 0) {
				W = 5.0 * V * (mu * E.sqlength() + 0.5 * lambda * E.trace() * E.trace());
				B = 5.0 * V * (2 * mu * F * E + lambda * E.trace() * F);
			} else {
				float J	   = F.det();
				float logJ = std::log(J);
				if (J < 0.0)
					logJ = 0.0;
				W = V * (0.5 * mu * (F.sqlength() - 3) - mu * logJ + 0.5 * lambda * logJ * logJ);
				B = V * (mu * F - mu * (F.inverse()).transpose() + lambda * logJ * (F.inverse()).transpose());
			}

			//std::cout << W << std::endl;

			if (W > 0.0) {
				float C = std::sqrt(2.0 * W);

				fmat3 BAt = B * A.transpose();

				fvec3 dC1 = (1.0 / C) * fvec3(BAt.m[0], BAt.m[3], BAt.m[6]);
				fvec3 dC2 = (1.0 / C) * fvec3(BAt.m[1], BAt.m[4], BAt.m[7]);
				fvec3 dC3 = (1.0 / C) * fvec3(BAt.m[2], BAt.m[5], BAt.m[8]);
				fvec3 dC0 = -(dC1 + dC2 + dC3);

				float dtdtdlambda = (-C - Lamdalist[i]) / ((dC0.sqlength() + dC1.sqlength() + dC2.sqlength() + dC3.sqlength()) / mass + 1.0 / (dt * dt));

				TempPosition[TetIndList[4 * i + 0]] = TempPosition[TetIndList[4 * i + 0]] + dtdtdlambda * (1.0 / mass) * dC0;
				TempPosition[TetIndList[4 * i + 1]] = TempPosition[TetIndList[4 * i + 1]] + dtdtdlambda * (1.0 / mass) * dC1;
				TempPosition[TetIndList[4 * i + 2]] = TempPosition[TetIndList[4 * i + 2]] + dtdtdlambda * (1.0 / mass) * dC2;
				TempPosition[TetIndList[4 * i + 3]] = TempPosition[TetIndList[4 * i + 3]] + dtdtdlambda * (1.0 / mass) * dC3;

				Lamdalist[i] += dtdtdlambda / (dt * dt);
			}
		}
	}

	void FemProjectJC(std::vector<fvec3>& TempPosition)
	{

		std::vector<fvec3> dx;
		dx.resize(4 * 6 * (N - 1) * (N - 1) * (N - 1));

#pragma omp parallel for shared(dx)
		for (uint32_t i = 0; i < 6 * (N - 1) * (N - 1) * (N - 1); i++) {
			fvec3 X0 = RestPositionList[TetIndList[4 * i + 0]];
			fvec3 X1 = RestPositionList[TetIndList[4 * i + 1]];
			fvec3 X2 = RestPositionList[TetIndList[4 * i + 2]];
			fvec3 X3 = RestPositionList[TetIndList[4 * i + 3]];

			fvec3 x0 = TempPosition[TetIndList[4 * i + 0]];
			fvec3 x1 = TempPosition[TetIndList[4 * i + 1]];
			fvec3 x2 = TempPosition[TetIndList[4 * i + 2]];
			fvec3 x3 = TempPosition[TetIndList[4 * i + 3]];

			fmat3 A = AList[i];
			fmat3 F = mat3(x1 - x0, x2 - x0, x3 - x0) * A;
			fmat3 E = 0.5 * (F.transpose() * F - fmat3::indentity());

			float V = VList[i];

			float W;
			fmat3 B;

			if (MaterialInd == 0) {
				W = 5.0 * V * (mu * E.sqlength() + 0.5 * lambda * E.trace() * E.trace());
				B = 5.0 * V * (2 * mu * F * E + lambda * E.trace() * F);
			} else {
				float J	   = F.det();
				float logJ = std::log(J);
				if (J < 0.0)
					logJ = 0.0;
				W = V * (0.5 * mu * (F.sqlength() - 3) - mu * logJ + 0.5 * lambda * logJ * logJ);
				B = V * (mu * F - mu * (F.inverse()).transpose() + lambda * logJ * (F.inverse()).transpose());
			}

			//std::cout << W << std::endl;

			if (W > 0.0) {
				float C = std::sqrt(2.0 * W);

				fmat3 BAt = B * A.transpose();

				fvec3 dC1 = (1.0 / C) * fvec3(BAt.m[0], BAt.m[3], BAt.m[6]);
				fvec3 dC2 = (1.0 / C) * fvec3(BAt.m[1], BAt.m[4], BAt.m[7]);
				fvec3 dC3 = (1.0 / C) * fvec3(BAt.m[2], BAt.m[5], BAt.m[8]);
				fvec3 dC0 = -(dC1 + dC2 + dC3);

				float dtdtdlambda = (-C - Lamdalist[i]) / ((dC0.sqlength() + dC1.sqlength() + dC2.sqlength() + dC3.sqlength()) / mass + 1.0 / (dt * dt));
				dtdtdlambda *= 1.2;

				dx[4 * i + 0] = dtdtdlambda * (1.0 / mass) * dC0;
				dx[4 * i + 1] = dtdtdlambda * (1.0 / mass) * dC1;
				dx[4 * i + 2] = dtdtdlambda * (1.0 / mass) * dC2;
				dx[4 * i + 3] = dtdtdlambda * (1.0 / mass) * dC3;

				Lamdalist[i] += dtdtdlambda / (dt * dt);
			} else {
				dx[4 * i + 0] = fvec3(0.0);
				dx[4 * i + 1] = fvec3(0.0);
				dx[4 * i + 2] = fvec3(0.0);
				dx[4 * i + 3] = fvec3(0.0);
			}
		}

		for (uint32_t i = 0; i < 4 * 6 * (N - 1) * (N - 1) * (N - 1); i++) {
			TempPosition[TetIndList[i]] = TempPosition[TetIndList[i]] + dx[i];
		}
	}

	void FixedProjection(std::vector<fvec3>& TempPosition)
	{

		//for (uint32_t i = N * N * (N - 1); i < N * N * N; i++) {
		//	TempPosition[i] = RestPositionList[i];
		//}

#pragma omp parallel
		for (uint32_t i = 0; i < N * N * N; i++) {
			if (TempPosition[i].y < -12.0)
				TempPosition[i].y = -12.0;
		}

		if (rightwall) {

#pragma omp parallel
			for (uint32_t y = 0; y < N; y++) {
				for (uint32_t z = 0; z < N; z++) {
					uint32_t Ind	    = N * N * y + N * z + N - 1;
					TempPosition[Ind]   = fvec3(rightrot, 0.0, 0.0).rotation() * (RestPositionList[Ind] - bias) + bias;
					TempPosition[Ind].x = bias.x + rightedge;
				}
			}

#pragma omp parallel
			for (uint32_t y = 0; y < N; y++) {
				for (uint32_t z = 0; z < N; z++) {
					uint32_t Ind	  = N * N * y + N * z;
					TempPosition[Ind] = RestPositionList[Ind];
				}
			}
		}
	}
};

namespace Physics {

int32_t solver = 0;

void timestep(CubeMesh& CM)
{

	CM.ClearLamda();

	uint32_t NodeSize = CM.N * CM.N * CM.N;
	std::vector<fvec3> tempp(NodeSize);

#pragma omp parallel
	for (uint32_t i = 0; i < NodeSize; i++) {
		fvec3 velocity = CM.VelocityList[i] + dt * fvec3(0.0, -9.8, 0.0);
		tempp[i]       = CM.PositionList[i] + dt * velocity;
	}

	if (solver == 0)
		for (uint32_t x = 0; x < 5; x++) {
			CM.FemProjectGS(tempp);
			CM.FixedProjection(tempp);
		}
	else
		for (uint32_t x = 0; x < 5; x++) {
			CM.FemProjectJC(tempp);
			CM.FixedProjection(tempp);
		}

#pragma omp parallel
	for (uint32_t i = 0; i < NodeSize; i++) {
		CM.VelocityList[i] = (tempp[i] - CM.PositionList[i]) / dt;
		CM.VelocityList[i] = 0.99 * CM.VelocityList[i];
		CM.PositionList[i] = tempp[i];
	}
}
};

int main(int argc, char const* argv[])
{
	using namespace std;

	Visualizer::Init(1024, 1024);

	IMGUI_CHECKVERSION();
	ImGui::CreateContext();
	ImGuiIO& io    = ImGui::GetIO();
	io.IniFilename = NULL;
	(void)io;

	ImGui::StyleColorsDark();

	ImGui_ImplOpenGL3_Init("#version 460 core");
	ImGui_ImplGlfw_InitForOpenGL(Visualizer::GetWindowPtr(), true);

	Renderer3D::Init();

	omp_set_num_threads(8);

	uint32_t cagelist[24] = { 0, 1, 1, 2, 2, 3, 3, 0, 4, 5, 5, 6,
		6, 7, 7, 4, 0, 4, 3, 7, 1, 5, 2, 6 };
	linevertarray cage(8, nullptr, 24, cagelist);
	cage.setposition(0, -15.0f, -15.0f, -15.0f);
	cage.setposition(1, 15.0f, -15.0f, -15.0f);
	cage.setposition(2, 15.0f, 15.0f, -15.0f);
	cage.setposition(3, -15.0f, 15.0f, -15.0f);
	cage.setposition(4, -15.0f, -15.0f, 15.0f);
	cage.setposition(5, 15.0f, -15.0f, 15.0f);
	cage.setposition(6, 15.0f, 15.0f, 15.0f);
	cage.setposition(7, -15.0f, 15.0f, 15.0f);
	cage.settype(0);

	uint32_t floorlist[6] = { 0, 1, 2, 0, 2, 3 };
	trianglevertarray floor(6, nullptr, 6, floorlist);
	floor.setposition(0, -12.0f, -12.0f, -12.0f);
	floor.setposition(1, -12.0f, -12.0f, 12.0f);
	floor.setposition(2, 12.0f, -12.0f, 12.0f);
	floor.setposition(3, 12.0f, -12.0f, -12.0f);
	for (uint32_t i = 0; i < 4; i++) {
		floor.setnormal(i, 0.0f, 1.0f, 0.0f);
	}
	floor.setuv(0, 0.0f, 1.0f);
	floor.setuv(1, 0.0f, 0.0f);
	floor.setuv(2, 1.0f, 0.0f);
	floor.setuv(3, 1.0f, 1.0f);
	floor.settype(2);

	uint8_t* images0 = new uint8_t[1024 * 4 * 1024];
	for (uint32_t i = 0; i < 1024; i++) {
		for (uint32_t j = 0; j < 1024; j++) {
			if ((((i / 32) % 2) + ((j / 32) % 2)) == 1) {
				images0[1024 * 4 * i + 4 * j + 0] = 255;
				images0[1024 * 4 * i + 4 * j + 1] = 255;
				images0[1024 * 4 * i + 4 * j + 2] = 255;
				images0[1024 * 4 * i + 4 * j + 3] = 255;
			} else {
				images0[1024 * 4 * i + 4 * j + 0] = 50;
				images0[1024 * 4 * i + 4 * j + 1] = 50;
				images0[1024 * 4 * i + 4 * j + 2] = 50;
				images0[1024 * 4 * i + 4 * j + 3] = 255;
			}
		}
	}
	texture simasima0(1024, 1024, images0);

	uint8_t* images1 = new uint8_t[1024 * 4 * 1024];
	for (uint32_t i = 0; i < 1024; i++) {
		for (uint32_t j = 0; j < 1024; j++) {
			if ((((i / 64) % 2) + ((j / 64) % 2)) == 1) {
				images1[1024 * 4 * i + 4 * j + 0] = 255;
				images1[1024 * 4 * i + 4 * j + 1] = 255;
				images1[1024 * 4 * i + 4 * j + 2] = 255;
				images1[1024 * 4 * i + 4 * j + 3] = 255;
			} else {
				images1[1024 * 4 * i + 4 * j + 0] = 50;
				images1[1024 * 4 * i + 4 * j + 1] = 50;
				images1[1024 * 4 * i + 4 * j + 2] = 50;
				images1[1024 * 4 * i + 4 * j + 3] = 255;
			}
		}
	}
	texture simasima1(1024, 1024, images1);

	fvec3 camerap(0.0, 15.0, 20.0);

	//myphysics.rbodyvec.reserve(10);

	CubeMesh CM0(10, 5.0, fvec3(-8.0, -5.0, 0.0));

	//init

	std::vector<Renderer3D::drawobject> shadowlist;
	std::vector<Renderer3D::drawobject> edgelist;
	std::vector<Renderer3D::drawobject> renderlist;

	shadowlist.emplace_back(Renderer3D::drawobject { floor, nullptr, nullptr, nullptr });
	shadowlist.emplace_back(Renderer3D::drawobject { CM0.tva, nullptr, nullptr, nullptr });

	renderlist.emplace_back(Renderer3D::drawobject { floor, &simasima0, nullptr, nullptr });
	renderlist.emplace_back(Renderer3D::drawobject { cage, nullptr, nullptr, nullptr });
	renderlist.emplace_back(Renderer3D::drawobject { CM0.tva, nullptr, nullptr, nullptr });
	renderlist.emplace_back(Renderer3D::drawobject { CM0.lva, nullptr, nullptr, nullptr });

	//rendering loop

	double ctime = 0.0;
	double vtime = 0.0;

	while (!Visualizer::Is_Closed()) {
		Renderer3D::Clear();
		//imgui reset
		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();
		ImGui::SetNextWindowSize(ImVec2(300, 400), ImGuiCond_FirstUseEver);

		static bool is_stop = true;

		static float mousemovex = 0.5 * 3.14;
		static float mousemovey = 0.5 * 3.14;
		static float cameraL	= 30.0f;

		float camerasinp;
		float cameracosp;
		float camerasint;
		float cameracost;

		static bool nextframe = false;

		//physics
		ctime = Visualizer::GetTime();
		if (!is_stop || nextframe) {
			Physics::timestep(CM0);
			vtime += dt;
			nextframe = false;

			CM0.UpdataVertarray();
			CM0.Setdata();
		}

		//imgui
		{

			ImGui::Begin("SoftCube");

			ImGui::Text("FPS : %.1f", ImGui::GetIO().Framerate);

			ImGui::Checkbox("stop", &is_stop);
			if (ImGui::Button("nextframe")) {
				nextframe = true;
			}

			ImVec2 mousemove = ImGui::GetMouseDragDelta(1);

			mousemovex += mousemove.x / (1024 * 5);
			mousemovey += mousemove.y / (1024 * 5);
			ImGui::Text(" theta = %.1f", mousemovey);
			ImGui::Text(" phi = %.1f", mousemovex);

			camerasinp = std::sin(mousemovex);
			cameracosp = std::cos(mousemovex);
			camerasint = std::sin(mousemovey);
			cameracost = std::cos(mousemovey);

			float dlength = io.MouseWheel;
			cameraL += dlength * 1.0f;
			if (cameraL < 10.0)
				cameraL = 10.0;
			if (cameraL > 80.0)
				cameraL = 80.0;

			camerap.x = cameraL * camerasint * cameracosp;
			camerap.z = cameraL * camerasint * camerasinp;
			camerap.y = cameraL * cameracost;

			ImGui::Text(" x = %.1f", camerap.x);
			ImGui::Text(" y = %.1f", camerap.y);
			ImGui::Text(" z = %.1f", camerap.z);

			ImGui::Text("camera length = %.1f", cameraL);

			ImGui::Combo("Solver", &(Physics::solver), "Gauss-Seidel\0Jacobi\0\0");

			ImGui::Combo("Material", &(CM0.MaterialInd), "Saint Venant-Kirchhoff\0Neo-Hookean\0\0");
			ImGui::SliderFloat("mu", &mu, 0.0f, 500.0f, "%.4f", ImGuiSliderFlags_Logarithmic);
			ImGui::SliderFloat("lambda", &lambda, 0.0f, 500.0f, "%.4f", ImGuiSliderFlags_Logarithmic);

			ImGui::Checkbox("right edge", &rightwall);
			ImGui::SliderFloat("rightedge", &rightedge, 0.0f, 20.0f);
			ImGui::SliderFloat("rightrot", &rightrot, 0.001f, 2.0 * 3.1415f);

			ImGui::Text("realtime = %.1f", ctime);
			ImGui::Text("virtualtime = %.1f", vtime);

			if (ImGui::Button("reset")) {

				vtime = 0.0;

				for (uint32_t i = 0; i < (CM0.N * CM0.N * CM0.N); i++) {
					CM0.PositionList[i] = CM0.RestPositionList[i];
					CM0.VelocityList[i] = fvec3(0.0);
				}
			}

			ImGui::End();
		}

		//renderer set

		Renderer3D::setcposi(camerap);
		Renderer3D::updateUniformobj();

		//rendering

		Renderer3D::Draw(shadowlist, edgelist, renderlist);

		ImGui::Render();
		ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

		//swapbuff
		Visualizer::SwapBuffer();
		//wait event
		Visualizer::PollEvent();
	}

	return 0;
}
