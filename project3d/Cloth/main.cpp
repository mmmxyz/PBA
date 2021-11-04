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

float mu      = 40;
float lambda  = 80;
float bendCof = 1.0;
float rho     = 0.01;

class ClothMesh {
    public:
	std::vector<fvec3> PositionList;
	std::vector<fvec2> RestPositionList;
	std::vector<fvec3> VelocityList;
	std::vector<uint32_t> TriIndList;
	std::vector<uint32_t> InnerEdgeIndList;
	std::vector<fvec4> InnerEdgeCList;

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

	std::vector<float> ElasticLamdalist;
	std::vector<float> BendLamdalist;

	std::vector<fmat2> AList;
	std::vector<float> VList;

	float mass;

	const uint32_t N;
	ClothMesh(uint32_t N, float length, const fvec3& bias)
	    : N(N)
	    , lva()
	    , tva()
	{

		mass /= N * N;

		PositionList.reserve(N * N);
		RestPositionList.reserve(N * N);
		VelocityList.reserve(N * N);
		TriIndList.reserve(3 * 2 * (N - 1) * (N - 1));
		InnerEdgeIndList.reserve(4 * ((N - 1) * (N - 1) + 2 * (N - 2) * (N - 2)));
		InnerEdgeIndList.reserve((N - 1) * (N - 1) + 2 * (N - 2) * (N - 2));

		AList.reserve(2 * (N - 1) * (N - 1));
		VList.reserve(2 * (N - 1) * (N - 1));

		ElasticLamdalist.resize(3 * 2 * (N - 1) * (N - 1));
		BendLamdalist.resize((N - 1) * (N - 1) + 2 * (N - 2) * (N - 2));

		mass = length * length * rho / (N * N);

		for (uint32_t y = 0; y < N; y++) {
			float vy = -0.5 * length + length * (y / float(N - 1));
			for (uint32_t x = 0; x < N; x++) {
				float vx = -0.5 * length + length * (x / float(N - 1));
				PositionList.emplace_back(fvec3(vx, vy, -vy + 0.5) + bias);
				RestPositionList.emplace_back(fvec2(vx, vy) + fvec2(bias.x, bias.y));
				VelocityList.emplace_back(fvec3(0.0));
			}
		}

		for (uint32_t y = 0; y < N - 1; y++) {
			for (uint32_t x = 0; x < N - 1; x++) {
				uint32_t Ind0 = N * y + x;
				uint32_t Ind1 = Ind0 + 1;
				uint32_t Ind2 = Ind0 + N;
				uint32_t Ind3 = Ind0 + N + 1;

				TriIndList.emplace_back(Ind0);
				TriIndList.emplace_back(Ind1);
				TriIndList.emplace_back(Ind2);

				TriIndList.emplace_back(Ind2);
				TriIndList.emplace_back(Ind1);
				TriIndList.emplace_back(Ind3);
			}
		}
		for (uint32_t y = 0; y < N - 1; y++) {
			for (uint32_t x = 0; x < N - 1; x++) {
				uint32_t Ind0 = N * y + x;
				uint32_t Ind1 = Ind0 + 1;
				uint32_t Ind2 = Ind0 + N;
				uint32_t Ind3 = Ind0 + N + 1;

				InnerEdgeIndList.emplace_back(Ind1);
				InnerEdgeIndList.emplace_back(Ind2);
				InnerEdgeIndList.emplace_back(Ind0);
				InnerEdgeIndList.emplace_back(Ind3);
			}
		}
		for (uint32_t y = 1; y < N - 1; y++) {
			for (uint32_t x = 1; x < N - 1; x++) {
				uint32_t Ind0 = N * y + x;
				uint32_t Ind1 = Ind0 + N;
				uint32_t Ind2 = Ind0 + N - 1;
				uint32_t Ind3 = Ind0 + 1;

				InnerEdgeIndList.emplace_back(Ind0);
				InnerEdgeIndList.emplace_back(Ind1);
				InnerEdgeIndList.emplace_back(Ind2);
				InnerEdgeIndList.emplace_back(Ind3);
			}
		}
		for (uint32_t y = 1; y < N - 1; y++) {
			for (uint32_t x = 1; x < N - 1; x++) {
				uint32_t Ind0 = N * y + x;
				uint32_t Ind1 = Ind0 + 1;
				uint32_t Ind2 = Ind0 + N;
				uint32_t Ind3 = Ind0 - N + 1;

				InnerEdgeIndList.emplace_back(Ind0);
				InnerEdgeIndList.emplace_back(Ind1);
				InnerEdgeIndList.emplace_back(Ind2);
				InnerEdgeIndList.emplace_back(Ind3);
			}
		}

		for (uint32_t i = 0; i < (N - 1) * (N - 1) + 2 * (N - 2) * (N - 2); i++) {
			fvec2 X0 = RestPositionList[InnerEdgeIndList[4 * i + 0]];
			fvec2 X1 = RestPositionList[InnerEdgeIndList[4 * i + 1]];
			fvec2 X2 = RestPositionList[InnerEdgeIndList[4 * i + 2]];
			fvec2 X3 = RestPositionList[InnerEdgeIndList[4 * i + 3]];

			float angle210 = std::acos(((X2 - X1).normalize()).dot((X0 - X1).normalize()));
			float angle201 = std::acos(((X2 - X0).normalize()).dot((X1 - X0).normalize()));
			float angle310 = std::acos(((X3 - X1).normalize()).dot((X0 - X1).normalize()));
			float angle301 = std::acos(((X3 - X0).normalize()).dot((X1 - X0).normalize()));

			float cot210 = 1.0 / std::tan(angle210);
			float cot201 = 1.0 / std::tan(angle201);
			float cot310 = 1.0 / std::tan(angle310);
			float cot301 = 1.0 / std::tan(angle301);

			InnerEdgeCList.emplace_back(fvec4(cot210 + cot310, cot301 + cot201, -cot210 - cot201, -cot310 - cot301));
		}

		for (uint32_t i = 0; i < 2 * (N - 1) * (N - 1); i++) {
			fvec2 X0 = RestPositionList[TriIndList[3 * i + 0]];
			fvec2 X1 = RestPositionList[TriIndList[3 * i + 1]];
			fvec2 X2 = RestPositionList[TriIndList[3 * i + 2]];

			AList.emplace_back(fmat2(X1 - X0, X2 - X0).inverse());
			VList.emplace_back((X1 - X0).cross(X2 - X0) / 2.0);
		}

		tvsize = 2 * N * N;
		tvdata = new vertex[tvsize];
		tisize = 2 * 2 * 3 * (N - 1) * (N - 1);
		tilist = new uint32_t[tisize];

		for (uint32_t y = 0; y < N - 1; y++) {
			for (uint32_t x = 0; x < N - 1; x++) {
				uint32_t Ind0			  = N * y + x;
				tilist[6 * ((N - 1) * y + x) + 0] = Ind0;
				tilist[6 * ((N - 1) * y + x) + 1] = Ind0 + 1;
				tilist[6 * ((N - 1) * y + x) + 2] = Ind0 + N;

				tilist[6 * ((N - 1) * y + x) + 3] = Ind0 + N;
				tilist[6 * ((N - 1) * y + x) + 4] = Ind0 + 1;
				tilist[6 * ((N - 1) * y + x) + 5] = Ind0 + N + 1;
			}
		}

		for (uint32_t y = 0; y < N - 1; y++) {
			for (uint32_t x = 0; x < N - 1; x++) {
				uint32_t Ind0						  = N * y + x;
				tilist[6 * (N - 1) * (N - 1) + 6 * ((N - 1) * y + x) + 0] = N * N + Ind0;
				tilist[6 * (N - 1) * (N - 1) + 6 * ((N - 1) * y + x) + 1] = N * N + Ind0 + N;
				tilist[6 * (N - 1) * (N - 1) + 6 * ((N - 1) * y + x) + 2] = N * N + Ind0 + 1;

				tilist[6 * (N - 1) * (N - 1) + 6 * ((N - 1) * y + x) + 3] = N * N + Ind0 + N;
				tilist[6 * (N - 1) * (N - 1) + 6 * ((N - 1) * y + x) + 4] = N * N + Ind0 + N + 1;
				tilist[6 * (N - 1) * (N - 1) + 6 * ((N - 1) * y + x) + 5] = N * N + Ind0 + 1;
			}
		}

		for (uint32_t i = 0; i < 2 * N * N; i++) {
			tvdata[i].color[0] = 0.8;
			tvdata[i].color[1] = 0.8;
			tvdata[i].color[2] = 0.8;
			tvdata[i].color[3] = 1.0;
			tvdata[i].type	   = 1;
		}

		lvsize = N * N;
		lvdata = new vertex[lvsize];
		lisize = (6 * (N - 1) + 2) * (N - 1) + 2 * (N - 1);
		lilist = new uint32_t[lisize];

		for (uint32_t y = 0; y < N - 1; y++) {
			for (uint32_t x = 0; x < N - 1; x++) {
				uint32_t Ind0				  = N * y + x;
				lilist[(6 * (N - 1) + 2) * y + 6 * x + 0] = Ind0;
				lilist[(6 * (N - 1) + 2) * y + 6 * x + 1] = Ind0 + 1;

				lilist[(6 * (N - 1) + 2) * y + 6 * x + 2] = Ind0 + 1;
				lilist[(6 * (N - 1) + 2) * y + 6 * x + 3] = Ind0 + N;

				lilist[(6 * (N - 1) + 2) * y + 6 * x + 4] = Ind0 + N;
				lilist[(6 * (N - 1) + 2) * y + 6 * x + 5] = Ind0;
			}
			lilist[(6 * (N - 1) + 2) * y + 6 * (N - 1) + 0] = N * y;
			lilist[(6 * (N - 1) + 2) * y + 6 * (N - 1) + 1] = N * y + N;
		}
		for (uint32_t x = 0; x < N - 1; x++) {
			lilist[(6 * (N - 1) + 2) * (N - 1) + 2 * x + 0] = N * (N - 1) + x;
			lilist[(6 * (N - 1) + 2) * (N - 1) + 2 * x + 1] = N * (N - 1) + x + 1;
		}

		for (uint32_t i = 0; i < N * N; i++) {
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

	void
	UpdataVertarray()
	{

		for (uint32_t y = 0; y < N; y++) {
			for (uint32_t x = 0; x < N; x++) {
				fvec3 v				      = PositionList[N * y + x];
				tvdata[N * y + x].position[0]	      = v.x;
				tvdata[N * y + x].position[1]	      = v.y;
				tvdata[N * y + x].position[2]	      = v.z;
				tvdata[N * N + N * y + x].position[0] = v.x;
				tvdata[N * N + N * y + x].position[1] = v.y;
				tvdata[N * N + N * y + x].position[2] = v.z;
				lvdata[N * y + x].position[0]	      = v.x;
				lvdata[N * y + x].position[1]	      = v.y;
				lvdata[N * y + x].position[2]	      = v.z;
			}
		}

		//normal

		std::vector<fvec3> NormalSet(2 * N * N);
		for (auto& x : NormalSet)
			x = fvec3(0.0);

		for (uint32_t i = 0; i < 4 * (N - 1) * (N - 1); i++) {
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

		for (uint32_t i = 0; i < 2 * N * N; i++) {
			tvdata[i].normal[0] = NormalSet[i].x;
			tvdata[i].normal[1] = NormalSet[i].y;
			tvdata[i].normal[2] = NormalSet[i].z;
			tvdata[i].position[0] += 0.005 * NormalSet[i].x;
			tvdata[i].position[1] += 0.005 * NormalSet[i].y;
			tvdata[i].position[2] += 0.005 * NormalSet[i].z;
		}
		for (uint32_t i = 0; i < N * N; i++) {
			lvdata[i].position[0] += 0.005 * NormalSet[i].x;
			lvdata[i].position[1] += 0.005 * NormalSet[i].y;
			lvdata[i].position[2] += 0.005 * NormalSet[i].z;
		}
	}

	void Setdata()
	{
		tva.setdata(tvdata);
		lva.setdata(lvdata);
	}

	void ClearLamda()
	{
		for (auto& x : ElasticLamdalist) {
			x = 0.0;
		}
		for (auto& x : BendLamdalist) {
			x = 0.0;
		}
	}

	void FemElasticProjectGS(std::vector<fvec3>& TempPosition)
	{
		for (uint32_t i = 0; i < 2 * (N - 1) * (N - 1); i++) {
			fvec2 X0 = RestPositionList[TriIndList[3 * i + 0]];
			fvec2 X1 = RestPositionList[TriIndList[3 * i + 1]];
			fvec2 X2 = RestPositionList[TriIndList[3 * i + 2]];

			fvec3 x0 = TempPosition[TriIndList[3 * i + 0]];
			fvec3 x1 = TempPosition[TriIndList[3 * i + 1]];
			fvec3 x2 = TempPosition[TriIndList[3 * i + 2]];

			fmat2 A	 = AList[i];
			fmat32 F = mat32(x1 - x0, x2 - x0) * A;

			float C00 = F.m[2 * 0 + 0] * F.m[2 * 0 + 0] + F.m[2 * 1 + 0] * F.m[2 * 1 + 0] + F.m[2 * 2 + 0] * F.m[2 * 2 + 0];
			float C01 = F.m[2 * 0 + 0] * F.m[2 * 0 + 1] + F.m[2 * 1 + 0] * F.m[2 * 1 + 1] + F.m[2 * 2 + 0] * F.m[2 * 2 + 1];
			float C11 = F.m[2 * 0 + 1] * F.m[2 * 0 + 1] + F.m[2 * 1 + 1] * F.m[2 * 1 + 1] + F.m[2 * 2 + 1] * F.m[2 * 2 + 1];

			fmat2 E;
			E.m[0] = 0.5 * (C00 - 1);
			E.m[1] = 0.5 * C01;
			E.m[2] = 0.5 * C01;
			E.m[3] = 0.5 * (C11 - 1);

			float V = VList[i];
			V *= 0.01f;

			float W	 = V * (mu * E.sqlength() + 0.5 * lambda * E.trace() * E.trace());
			fmat32 B = V * (2 * mu * F * E + lambda * E.trace() * F);

			//std::cout << W << std::endl;
			//std::cout << V << std::endl;

			if (W > 0.0) {
				float C = std::sqrt(2.0 * W);

				fmat32 BAt = B * A.transpose();

				fvec3 dC1 = (1.0 / C) * fvec3(BAt.m[0], BAt.m[2], BAt.m[4]);
				fvec3 dC2 = (1.0 / C) * fvec3(BAt.m[1], BAt.m[3], BAt.m[5]);
				fvec3 dC0 = -(dC1 + dC2);

				//std::cout << dC1 << std::endl;

				float dtdtdlambda = (-C - ElasticLamdalist[i]) / ((dC0.sqlength() + dC1.sqlength() + dC2.sqlength()) / mass + 1.0 / (dt * dt));

				TempPosition[TriIndList[3 * i + 0]] = TempPosition[TriIndList[3 * i + 0]] + dtdtdlambda * (1.0 / mass) * dC0;
				TempPosition[TriIndList[3 * i + 1]] = TempPosition[TriIndList[3 * i + 1]] + dtdtdlambda * (1.0 / mass) * dC1;
				TempPosition[TriIndList[3 * i + 2]] = TempPosition[TriIndList[3 * i + 2]] + dtdtdlambda * (1.0 / mass) * dC2;

				ElasticLamdalist[i] += dtdtdlambda / (dt * dt);
			}
		}
	}

	void FemElasticProjectJC(std::vector<fvec3>& TempPosition)
	{
		std::vector<fvec3> dx;
		dx.resize(3 * 2 * (N - 1) * (N - 1));

#pragma omp parallel
		for (uint32_t i = 0; i < 2 * (N - 1) * (N - 1); i++) {
			fvec2 X0 = RestPositionList[TriIndList[3 * i + 0]];
			fvec2 X1 = RestPositionList[TriIndList[3 * i + 1]];
			fvec2 X2 = RestPositionList[TriIndList[3 * i + 2]];

			fvec3 x0 = TempPosition[TriIndList[3 * i + 0]];
			fvec3 x1 = TempPosition[TriIndList[3 * i + 1]];
			fvec3 x2 = TempPosition[TriIndList[3 * i + 2]];

			fmat2 A	 = AList[i];
			fmat32 F = mat32(x1 - x0, x2 - x0) * A;

			float C00 = F.m[2 * 0 + 0] * F.m[2 * 0 + 0] + F.m[2 * 1 + 0] * F.m[2 * 1 + 0] + F.m[2 * 2 + 0] * F.m[2 * 2 + 0];
			float C01 = F.m[2 * 0 + 0] * F.m[2 * 0 + 1] + F.m[2 * 1 + 0] * F.m[2 * 1 + 1] + F.m[2 * 2 + 0] * F.m[2 * 2 + 1];
			float C11 = F.m[2 * 0 + 1] * F.m[2 * 0 + 1] + F.m[2 * 1 + 1] * F.m[2 * 1 + 1] + F.m[2 * 2 + 1] * F.m[2 * 2 + 1];

			fmat2 E;
			E.m[0] = 0.5 * (C00 - 1);
			E.m[1] = 0.5 * C01;
			E.m[2] = 0.5 * C01;
			E.m[3] = 0.5 * (C11 - 1);

			float V = VList[i];
			V *= 0.01f;

			float W	 = V * (mu * E.sqlength() + 0.5 * lambda * E.trace() * E.trace());
			fmat32 B = V * (2 * mu * F * E + lambda * E.trace() * F);

			//std::cout << W << std::endl;
			//std::cout << V << std::endl;

			if (W > 0.0) {
				float C = std::sqrt(2.0 * W);

				fmat32 BAt = B * A.transpose();

				fvec3 dC1 = (1.0 / C) * fvec3(BAt.m[0], BAt.m[2], BAt.m[4]);
				fvec3 dC2 = (1.0 / C) * fvec3(BAt.m[1], BAt.m[3], BAt.m[5]);
				fvec3 dC0 = -(dC1 + dC2);

				//std::cout << dC1 << std::endl;

				float dtdtdlambda = (-C - ElasticLamdalist[i]) / ((dC0.sqlength() + dC1.sqlength() + dC2.sqlength()) / mass + 1.0 / (dt * dt));

				dx[3 * i + 0] = dtdtdlambda * (1.0 / mass) * dC0;
				dx[3 * i + 1] = dtdtdlambda * (1.0 / mass) * dC1;
				dx[3 * i + 2] = dtdtdlambda * (1.0 / mass) * dC2;

				ElasticLamdalist[i] += dtdtdlambda / (dt * dt);
			} else {
				dx[3 * i + 0] = fvec3(0.0);
				dx[3 * i + 1] = fvec3(0.0);
				dx[3 * i + 2] = fvec3(0.0);
			}
		}

#pragma omp parallel
		for (uint32_t i = 0; i < 3 * 2 * (N - 1) * (N - 1); i++) {
			TempPosition[TriIndList[i]] = TempPosition[TriIndList[i]] + dx[i];
		}
	}

	void FemBendProjectionGS(std::vector<fvec3>& TempPosition)
	{

		for (uint32_t i = 0; i < (N - 1) * (N - 1) + 2 * (N - 2) * (N - 2); i++) {
			fvec3 x0 = TempPosition[InnerEdgeIndList[4 * i + 0]];
			fvec3 x1 = TempPosition[InnerEdgeIndList[4 * i + 1]];
			fvec3 x2 = TempPosition[InnerEdgeIndList[4 * i + 2]];
			fvec3 x3 = TempPosition[InnerEdgeIndList[4 * i + 3]];

			fvec4 Cot = InnerEdgeCList[i];

			fmat4 X	   = fmat4(fvec4(x0), fvec4(x1), fvec4(x2), fvec4(x3));
			fvec4 XCot = X * Cot;

			XCot = bendCof * XCot;

			float Q = XCot.x * XCot.x + XCot.y * XCot.y + XCot.z * XCot.z;

			std::cout << Q << std::endl;

			if (Q > 0.0) {
				float C = std::sqrt(Q);

				fvec3 dC0 = fvec3(
				    XCot.x * Cot.x,
				    XCot.y * Cot.x,
				    XCot.z * Cot.x);
				fvec3 dC1 = fvec3(
				    XCot.x * Cot.y,
				    XCot.y * Cot.y,
				    XCot.z * Cot.y);
				fvec3 dC2 = fvec3(
				    XCot.x * Cot.z,
				    XCot.y * Cot.z,
				    XCot.z * Cot.z);
				fvec3 dC3 = fvec3(
				    XCot.x * Cot.w,
				    XCot.y * Cot.w,
				    XCot.z * Cot.w);

				float dtdtdlambda = (-C - BendLamdalist[i]) / ((dC0.sqlength() + dC1.sqlength() + dC2.sqlength() + dC3.sqlength()) / mass + 1.0 / (dt * dt));

				TempPosition[InnerEdgeIndList[4 * i + 0]] = TempPosition[InnerEdgeIndList[4 * i + 0]] + dtdtdlambda * (1.0 / mass) * dC0;
				TempPosition[InnerEdgeIndList[4 * i + 1]] = TempPosition[InnerEdgeIndList[4 * i + 1]] + dtdtdlambda * (1.0 / mass) * dC1;
				TempPosition[InnerEdgeIndList[4 * i + 2]] = TempPosition[InnerEdgeIndList[4 * i + 2]] + dtdtdlambda * (1.0 / mass) * dC2;
				TempPosition[InnerEdgeIndList[4 * i + 3]] = TempPosition[InnerEdgeIndList[4 * i + 3]] + dtdtdlambda * (1.0 / mass) * dC3;

				BendLamdalist[i] += dtdtdlambda / (dt * dt);
			}
		}
	}

	void FemBendProjectionJC(std::vector<fvec3>& TempPosition)
	{

		for (uint32_t i = 0; i < (N - 1) * (N - 1) + 2 * (N - 2) * (N - 2); i++) {
			fvec3 x0 = TempPosition[InnerEdgeIndList[4 * i + 0]];
			fvec3 x1 = TempPosition[InnerEdgeIndList[4 * i + 1]];
			fvec3 x2 = TempPosition[InnerEdgeIndList[4 * i + 2]];
			fvec3 x3 = TempPosition[InnerEdgeIndList[4 * i + 3]];

			fvec4 Cot = InnerEdgeCList[i];

			fmat4 X	   = fmat4(fvec4(x0), fvec4(x1), fvec4(x2), fvec4(x3));
			fvec4 XCot = X * Cot;

			XCot = bendCof * XCot;

			float Q = XCot.x * XCot.x + XCot.y * XCot.y + XCot.z * XCot.z;

			std::cout << Q << std::endl;

			if (Q > 0.0) {
				float C = std::sqrt(Q);

				fvec3 dC0 = fvec3(
				    XCot.x * Cot.x,
				    XCot.y * Cot.x,
				    XCot.z * Cot.x);
				fvec3 dC1 = fvec3(
				    XCot.x * Cot.y,
				    XCot.y * Cot.y,
				    XCot.z * Cot.y);
				fvec3 dC2 = fvec3(
				    XCot.x * Cot.z,
				    XCot.y * Cot.z,
				    XCot.z * Cot.z);
				fvec3 dC3 = fvec3(
				    XCot.x * Cot.w,
				    XCot.y * Cot.w,
				    XCot.z * Cot.w);

				float dtdtdlambda = (-C - BendLamdalist[i]) / ((dC0.sqlength() + dC1.sqlength() + dC2.sqlength() + dC3.sqlength()) / mass + 1.0 / (dt * dt));

				TempPosition[InnerEdgeIndList[4 * i + 0]] = TempPosition[InnerEdgeIndList[4 * i + 0]] + dtdtdlambda * (1.0 / mass) * dC0;
				TempPosition[InnerEdgeIndList[4 * i + 1]] = TempPosition[InnerEdgeIndList[4 * i + 1]] + dtdtdlambda * (1.0 / mass) * dC1;
				TempPosition[InnerEdgeIndList[4 * i + 2]] = TempPosition[InnerEdgeIndList[4 * i + 2]] + dtdtdlambda * (1.0 / mass) * dC2;
				TempPosition[InnerEdgeIndList[4 * i + 3]] = TempPosition[InnerEdgeIndList[4 * i + 3]] + dtdtdlambda * (1.0 / mass) * dC3;

				BendLamdalist[i] += dtdtdlambda / (dt * dt);
			}
		}
	}

	void FixedProjection(std::vector<fvec3>& TempPosition)
	{

		//for (uint32_t i = N * (N - 1); i < N * N; i++) {
		//	TempPosition[i].x = RestPositionList[i].x;
		//	TempPosition[i].y = RestPositionList[i].y;
		//	TempPosition[i].z = 0.0f;
		//}

		TempPosition[N * (N - 1)].x = RestPositionList[N * (N - 1)].x;
		TempPosition[N * (N - 1)].y = RestPositionList[N * (N - 1)].y;
		TempPosition[N * (N - 1)].z = 0.0f;
		TempPosition[N * N - 1].x   = RestPositionList[N * N - 1].x;
		TempPosition[N * N - 1].y   = RestPositionList[N * N - 1].y;
		TempPosition[N * N - 1].z   = 0.0f;

		//for (uint32_t i = 0; i < N * N; i++) {
		//	if (TempPosition[i].y < -12.0)
		//		TempPosition[i].y = -12.0;
		//}
	}
};

namespace Physics {

int32_t solver = 0;

void timestep(ClothMesh& CM)
{

	CM.ClearLamda();

	uint32_t NodeSize = CM.N * CM.N;
	std::vector<fvec3> tempp(NodeSize);

#pragma omp parallel
	for (uint32_t i = 0; i < NodeSize; i++) {
		fvec3 velocity = CM.VelocityList[i] + dt * fvec3(0.0, -9.8, 0.0);
		tempp[i]       = CM.PositionList[i] + dt * velocity;
	}

#pragma omp parallel
	if (solver == 0)
		for (uint32_t x = 0; x < 5; x++) {
			CM.FixedProjection(tempp);
			CM.FemElasticProjectGS(tempp);
			CM.FemBendProjection(tempp);
		}
	else
		for (uint32_t x = 0; x < 5; x++) {
			CM.FixedProjection(tempp);
			CM.FemElasticProjectJC(tempp);
			CM.FemBendProjection(tempp);
		}

#pragma omp parallel
	for (uint32_t i = 0; i < NodeSize; i++) {
		CM.VelocityList[i] = (tempp[i] - CM.PositionList[i]) / dt;
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

	ClothMesh CM0(30, 8.0, fvec3(0.0, 0.0, 0.0));

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

			ImGui::Begin("Cloth");

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

			ImGui::SliderFloat("mu", &mu, 0.0f, 500.0f, "%.4f", ImGuiSliderFlags_Logarithmic);
			ImGui::SliderFloat("lambda", &lambda, 0.0f, 500.0f, "%.4f", ImGuiSliderFlags_Logarithmic);
			ImGui::SliderFloat("Bending", &bendCof, 0.0f, 500.0f, "%.4f", ImGuiSliderFlags_Logarithmic);

			ImGui::Text("realtime = %.1f", ctime);
			ImGui::Text("virtualtime = %.1f", vtime);

			if (ImGui::Button("reset")) {

				vtime = 0.0;

				for (uint32_t i = 0; i < (CM0.N * CM0.N); i++) {
					CM0.PositionList[i].x = CM0.RestPositionList[i].x;
					CM0.PositionList[i].y = CM0.RestPositionList[i].y;
					CM0.PositionList[i].z = CM0.RestPositionList[i].y + 0.5;
					CM0.VelocityList[i]   = fvec3(0.0);
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
