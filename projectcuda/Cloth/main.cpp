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

#include "utils/meshgenerator/meshgenerator.hpp"

#include "utils/fem/fem.hpp"

#include "kernel.hpp"

#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>

constexpr uint32_t FPS = 60;
constexpr float dt     = 1.0 / FPS;
//glfwSwapInterval(1); なので60FPS以下になる。
//これはモニタやGPUに依るかもしれない。

float mu      = 90;
float lambda  = 50;
float bendCof = 0.00100;

float SpringCof = 3.0;

float rho = 0.005;

float edgedist = 0.3;

class ClothMesh {
    public:
	std::vector<fvec3> PositionList;
	std::vector<fvec2> RestPositionList;
	std::vector<fvec3> VelocityList;
	std::vector<fvec3> TempPositionList;

	std::vector<fvec3> dx;

	uint32_t vertsize;
	uint32_t* tilist;
	uint32_t* tilist2;
	uint32_t tisize;
	uint32_t* edgelist;
	uint32_t edgesize;
	uint32_t* InnerEdgeList;
	fvec4* InnerEdgeCList;
	uint32_t InnerEdgesize;

	linevertarray lva;

	trianglevertarray tva;
	trianglevertarray tva2;

	std::vector<float> ElasticLamdalist;
	std::vector<float> MassSpringLamdalist;
	std::vector<float> BendLamdalist;

	std::vector<fmat2> AList;
	std::vector<float> VList;

	std::vector<uint32_t> fix0;
	std::vector<uint32_t> fix1;

	float mass;

	const uint32_t N;

	ClothMesh(uint32_t N, float length, const fvec3& bias)
	    : N(N)
	    , lva()
	    , tva()
	    , tva2()
	{

		fvec3* vertdata;
		fvec2* Restvertdata;

		ClothFemMesh(N, length, &vertdata, &Restvertdata, vertsize, &tilist, tisize, &edgelist, edgesize, &InnerEdgeList, InnerEdgesize, &InnerEdgeCList);

		PositionList.reserve(vertsize);
		RestPositionList.reserve(vertsize);
		TempPositionList.reserve(vertsize);
		VelocityList.reserve(vertsize);

		dx.reserve(std::max(std::max(edgesize, tisize), InnerEdgesize));

		AList.reserve(tisize / 3);
		VList.reserve(tisize / 3);

		ElasticLamdalist.resize(tisize / 3);
		BendLamdalist.resize(InnerEdgesize / 4);
		MassSpringLamdalist.resize(edgesize / 2);

		mass = length * length * rho / (N * N);

		for (uint32_t i = 0; i < vertsize; i++) {
			PositionList[i]	    = vertdata[i];
			TempPositionList[i] = vertdata[i];
			RestPositionList[i] = Restvertdata[i];
			VelocityList[i]	    = fvec3(0.0);
		}
		for (uint32_t i = 0; i < vertsize; i++) {

			if ((PositionList[i] - PositionList[N * (N - 1)]).length() < 0.5)
				fix0.push_back(i);

			if ((PositionList[i] - PositionList[N * N - 1]).length() < 0.5)
				fix1.push_back(i);
		}

		tilist2 = new uint32_t[tisize];
		for (uint32_t i = 0; i < tisize / 3; i++) {
			tilist2[3 * i + 0] = tilist[3 * i + 0];
			tilist2[3 * i + 1] = tilist[3 * i + 2];
			tilist2[3 * i + 2] = tilist[3 * i + 1];
		}

		for (uint32_t i = 0; i < tisize / 3; i++) {
			fvec2 X0 = Restvertdata[tilist[3 * i + 0]];
			fvec2 X1 = Restvertdata[tilist[3 * i + 1]];
			fvec2 X2 = Restvertdata[tilist[3 * i + 2]];

			AList.emplace_back(fmat2(X1 - X0, X2 - X0).inverse());
			VList.emplace_back((X1 - X0).cross(X2 - X0) / 2.0);
		}

		tva.resetvertarray(vertsize, tisize, tilist);
		tva2.resetvertarray(vertsize, tisize, tilist2);

		tva.settype(1);
		tva.setcolor(0.8, 0.8, 0.8, 1.0);
		tva2.settype(1);
		tva2.setcolor(0.8, 0.8, 0.8, 1.0);

		lva.resetvertarray(vertsize, edgesize, edgelist);
		lva.settype(0);
		lva.setcolor(0.0, 0.0, 0.0, 1.0);

		this->UpdataVertarray();
		this->Setdata();

		delete[] vertdata;
		delete[] Restvertdata;
	}

	void
	UpdataVertarray()
	{

#pragma omp parallel for
		for (uint32_t i = 0; i < vertsize; i++) {
			fvec3 v		    = PositionList[i];
			tva[i].position[0]  = v.x;
			tva[i].position[1]  = v.y;
			tva[i].position[2]  = v.z;
			tva2[i].position[0] = v.x;
			tva2[i].position[1] = v.y;
			tva2[i].position[2] = v.z;
			lva[i].position[0]  = v.x;
			lva[i].position[1]  = v.y;
			lva[i].position[2]  = v.z;
		}

		//normal

		std::vector<fvec3> NormalSet(vertsize);
#pragma omp parallel for
		for (auto& x : NormalSet)
			x = fvec3(0.0);

		for (uint32_t i = 0; i < tisize / 3; i++) {
			fvec3 v0 = fvec3(tva[tilist[3 * i + 0]].position);
			fvec3 v1 = fvec3(tva[tilist[3 * i + 1]].position);
			fvec3 v2 = fvec3(tva[tilist[3 * i + 2]].position);

			fvec3 normal = (v1 - v0).cross(v2 - v0);
			if (normal.sqlength() > 0.000001)
				normal = normal.normalize();

			NormalSet[tilist[3 * i + 0]] = NormalSet[tilist[3 * i + 0]] + normal;
			NormalSet[tilist[3 * i + 1]] = NormalSet[tilist[3 * i + 1]] + normal;
			NormalSet[tilist[3 * i + 2]] = NormalSet[tilist[3 * i + 2]] + normal;
		}

#pragma omp parallel for
		for (uint32_t i = 0; i < vertsize; i++) {
			tva[i].normal[0] = NormalSet[i].x;
			tva[i].normal[1] = NormalSet[i].y;
			tva[i].normal[2] = NormalSet[i].z;
			tva[i].position[0] += 0.001 * NormalSet[i].x;
			tva[i].position[1] += 0.001 * NormalSet[i].y;
			tva[i].position[2] += 0.001 * NormalSet[i].z;

			tva2[i].normal[0] = -NormalSet[i].x;
			tva2[i].normal[1] = -NormalSet[i].y;
			tva2[i].normal[2] = -NormalSet[i].z;
			tva2[i].position[0] -= 0.001 * NormalSet[i].x;
			tva2[i].position[1] -= 0.001 * NormalSet[i].y;
			tva2[i].position[2] -= 0.001 * NormalSet[i].z;
		}

#pragma omp parallel for
		for (uint32_t i = 0; i < vertsize; i++) {
			lva[i].position[0] += 0.005 * NormalSet[i].x;
			lva[i].position[1] += 0.005 * NormalSet[i].y;
			lva[i].position[2] += 0.005 * NormalSet[i].z;
		}
	}

	void Setdata()
	{
		tva.vboupdate();
		tva2.vboupdate();
		lva.vboupdate();
	}

	void ClearLamda()
	{
#pragma omp parallel for
		for (auto& x : ElasticLamdalist) {
			x = 0.0;
		}
#pragma omp parallel for
		for (auto& x : BendLamdalist) {
			x = 0.0;
		}
#pragma omp parallel for
		for (auto& x : MassSpringLamdalist) {
			x = 0.0;
		}
	}

	void MassSpringSolverGS(std::vector<fvec3>& TempPosition)
	{

		for (uint32_t i = 0; i < edgesize / 2; i++) {
			fvec2 X0 = RestPositionList[edgelist[2 * i + 0]];
			fvec2 X1 = RestPositionList[edgelist[2 * i + 1]];

			fvec3 x0 = TempPosition[edgelist[2 * i + 0]];
			fvec3 x1 = TempPosition[edgelist[2 * i + 1]];

			float C = SpringCof * ((x1 - x0).sqlength() - (X1 - X0).sqlength());

			if (std::abs(C) > 0.001) {

				fvec3 dC0 = SpringCof * (x0 - x1);
				fvec3 dC1 = SpringCof * (x1 - x0);

				float dtdtlambda = (-C - MassSpringLamdalist[i]) / ((dC0.sqlength() + dC1.sqlength()) / mass + 1.0 / (dt * dt));

				TempPosition[edgelist[2 * i + 0]] = TempPosition[edgelist[2 * i + 0]] + dtdtlambda * (1.0 / mass) * dC0;
				TempPosition[edgelist[2 * i + 1]] = TempPosition[edgelist[2 * i + 1]] + dtdtlambda * (1.0 / mass) * dC1;

				MassSpringLamdalist[i] += dtdtlambda / (dt * dt);
			}
		}
	}

	void MassSpringSolverJC(std::vector<fvec3>& TempPosition)
	{

#pragma omp parallel for
		for (uint32_t i = 0; i < edgesize / 2; i++) {
			fvec2 X0 = RestPositionList[edgelist[2 * i + 0]];
			fvec2 X1 = RestPositionList[edgelist[2 * i + 1]];

			fvec3 x0 = TempPosition[edgelist[2 * i + 0]];
			fvec3 x1 = TempPosition[edgelist[2 * i + 1]];

			float C = SpringCof * ((x1 - x0).sqlength() - (X1 - X0).sqlength());

			if (std::abs(C) > 0.001) {

				fvec3 dC0 = SpringCof * (x0 - x1);
				fvec3 dC1 = SpringCof * (x1 - x0);

				float dtdtlambda = (-C - MassSpringLamdalist[i]) / ((dC0.sqlength() + dC1.sqlength()) / mass + 1.0 / (dt * dt));
				dtdtlambda *= 0.4;

				dx[2 * i + 0] = dtdtlambda * (1.0 / mass) * dC0;
				dx[2 * i + 1] = dtdtlambda * (1.0 / mass) * dC1;

				MassSpringLamdalist[i] += dtdtlambda / (dt * dt);
			} else {
				dx[2 * i + 0] = fvec3(0.0);
				dx[2 * i + 1] = fvec3(0.0);
			}
		}

		for (uint32_t i = 0; i < edgesize; i++) {
			TempPosition[edgelist[i]] = TempPosition[edgelist[i]] + dx[i];
		}
	}

	void FemElasticProjectGS(std::vector<fvec3>& TempPosition)
	{
		for (uint32_t i = 0; i < tisize / 3; i++) {
			fvec2 X0 = RestPositionList[tilist[3 * i + 0]];
			fvec2 X1 = RestPositionList[tilist[3 * i + 1]];
			fvec2 X2 = RestPositionList[tilist[3 * i + 2]];

			fvec3 x0 = TempPosition[tilist[3 * i + 0]];
			fvec3 x1 = TempPosition[tilist[3 * i + 1]];
			fvec3 x2 = TempPosition[tilist[3 * i + 2]];

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

			if (W > 0.0) {
				float C = std::sqrt(2.0 * W);

				fmat32 BAt = B * A.transpose();

				fvec3 dC1 = (1.0 / C) * fvec3(BAt.m[0], BAt.m[2], BAt.m[4]);
				fvec3 dC2 = (1.0 / C) * fvec3(BAt.m[1], BAt.m[3], BAt.m[5]);
				fvec3 dC0 = -(dC1 + dC2);

				float dtdtdlambda = (-C - ElasticLamdalist[i]) / ((dC0.sqlength() + dC1.sqlength() + dC2.sqlength()) / mass + 1.0 / (dt * dt));
				dtdtdlambda *= 1.0;

				TempPosition[tilist[3 * i + 0]] = TempPosition[tilist[3 * i + 0]] + dtdtdlambda * (1.0 / mass) * dC0;
				TempPosition[tilist[3 * i + 1]] = TempPosition[tilist[3 * i + 1]] + dtdtdlambda * (1.0 / mass) * dC1;
				TempPosition[tilist[3 * i + 2]] = TempPosition[tilist[3 * i + 2]] + dtdtdlambda * (1.0 / mass) * dC2;

				ElasticLamdalist[i] += dtdtdlambda / (dt * dt);
			}
		}
	}

	void FemElasticProjectJC(std::vector<fvec3>& TempPosition, const uint32_t X)
	{

#pragma omp parallel for
		for (uint32_t i = 0; i < tisize / 3; i++) {

			fvec2 X0 = RestPositionList[tilist[3 * i + 0]];
			fvec2 X1 = RestPositionList[tilist[3 * i + 1]];
			fvec2 X2 = RestPositionList[tilist[3 * i + 2]];

			fvec3 x0 = TempPosition[tilist[3 * i + 0]];
			fvec3 x1 = TempPosition[tilist[3 * i + 1]];
			fvec3 x2 = TempPosition[tilist[3 * i + 2]];

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

			if (W > 0.0) {
				float C = std::sqrt(2.0 * W);

				fmat32 BAt = B * A.transpose();

				fvec3 dC1 = (1.0 / C) * fvec3(BAt.m[0], BAt.m[2], BAt.m[4]);
				fvec3 dC2 = (1.0 / C) * fvec3(BAt.m[1], BAt.m[3], BAt.m[5]);
				fvec3 dC0 = -(dC1 + dC2);

				float dtdtdlambda = (-C - ElasticLamdalist[i]) / ((dC0.sqlength() + dC1.sqlength() + dC2.sqlength()) / mass + 1.0 / (dt * dt));
				dtdtdlambda *= (X) / 40.0;

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

		for (uint32_t i = 0; i < tisize; i++) {
			TempPosition[tilist[i]] = TempPosition[tilist[i]] + dx[i];
		}
	}

	void FemBendProjectionGS(std::vector<fvec3>& TempPosition)
	{

		for (uint32_t i = 0; i < InnerEdgesize / 4; i++) {
			fvec3 x0 = TempPosition[InnerEdgeList[4 * i + 0]];
			fvec3 x1 = TempPosition[InnerEdgeList[4 * i + 1]];
			fvec3 x2 = TempPosition[InnerEdgeList[4 * i + 2]];
			fvec3 x3 = TempPosition[InnerEdgeList[4 * i + 3]];

			fvec4 Cot = InnerEdgeCList[i];

			fmat4 X	   = fmat4(fvec4(x0), fvec4(x1), fvec4(x2), fvec4(x3));
			fvec4 XCot = X * Cot;

			float Q = bendCof * XCot.sqlength();

			if (Q > 0.0) {
				float C = std::sqrt(2.0 * Q);

				fvec3 dC0 = (1.0 / C) * bendCof * fvec3(XCot.x * Cot.x, XCot.y * Cot.x, XCot.z * Cot.x);
				fvec3 dC1 = (1.0 / C) * bendCof * fvec3(XCot.x * Cot.y, XCot.y * Cot.y, XCot.z * Cot.y);
				fvec3 dC2 = (1.0 / C) * bendCof * fvec3(XCot.x * Cot.z, XCot.y * Cot.z, XCot.z * Cot.z);
				fvec3 dC3 = (1.0 / C) * bendCof * fvec3(XCot.x * Cot.w, XCot.y * Cot.w, XCot.z * Cot.w);

				float dtdtdlambda = (-C - BendLamdalist[i]) / ((dC0.sqlength() + dC1.sqlength() + dC2.sqlength() + dC3.sqlength()) / mass + 1.0 / (dt * dt));

				TempPosition[InnerEdgeList[4 * i + 0]] = TempPosition[InnerEdgeList[4 * i + 0]] + dtdtdlambda * (1.0 / mass) * dC0;
				TempPosition[InnerEdgeList[4 * i + 1]] = TempPosition[InnerEdgeList[4 * i + 1]] + dtdtdlambda * (1.0 / mass) * dC1;
				TempPosition[InnerEdgeList[4 * i + 2]] = TempPosition[InnerEdgeList[4 * i + 2]] + dtdtdlambda * (1.0 / mass) * dC2;
				TempPosition[InnerEdgeList[4 * i + 3]] = TempPosition[InnerEdgeList[4 * i + 3]] + dtdtdlambda * (1.0 / mass) * dC3;

				BendLamdalist[i] += dtdtdlambda / (dt * dt);
			}
		}
	}

	void FemBendProjectionJC(std::vector<fvec3>& TempPosition)
	{

#pragma omp parallel for
		for (uint32_t i = 0; i < InnerEdgesize / 4; i++) {
			fvec3 x0 = TempPosition[InnerEdgeList[4 * i + 0]];
			fvec3 x1 = TempPosition[InnerEdgeList[4 * i + 1]];
			fvec3 x2 = TempPosition[InnerEdgeList[4 * i + 2]];
			fvec3 x3 = TempPosition[InnerEdgeList[4 * i + 3]];

			fvec4 Cot = InnerEdgeCList[i];

			fmat4 X	   = fmat4(fvec4(x0), fvec4(x1), fvec4(x2), fvec4(x3));
			fvec4 XCot = X * Cot;

			float Q = bendCof * XCot.sqlength();

			if (Q > 0.0) {
				float C = std::sqrt(2.0 * Q);

				fvec3 dC0 = (1.0 / C) * bendCof * fvec3(XCot.x * Cot.x, XCot.y * Cot.x, XCot.z * Cot.x);
				fvec3 dC1 = (1.0 / C) * bendCof * fvec3(XCot.x * Cot.y, XCot.y * Cot.y, XCot.z * Cot.y);
				fvec3 dC2 = (1.0 / C) * bendCof * fvec3(XCot.x * Cot.z, XCot.y * Cot.z, XCot.z * Cot.z);
				fvec3 dC3 = (1.0 / C) * bendCof * fvec3(XCot.x * Cot.w, XCot.y * Cot.w, XCot.z * Cot.w);

				float dtdtdlambda = (-C - BendLamdalist[i]) / ((dC0.sqlength() + dC1.sqlength() + dC2.sqlength() + dC3.sqlength()) / mass + 1.0 / (dt * dt));
				dtdtdlambda *= 0.5;

				dx[4 * i + 0] = dtdtdlambda * (1.0 / mass) * dC0;
				dx[4 * i + 1] = dtdtdlambda * (1.0 / mass) * dC1;
				dx[4 * i + 2] = dtdtdlambda * (1.0 / mass) * dC2;
				dx[4 * i + 3] = dtdtdlambda * (1.0 / mass) * dC3;

				BendLamdalist[i] += dtdtdlambda / (dt * dt);
			} else {
				dx[4 * i + 0] = fvec3(0.0);
				dx[4 * i + 1] = fvec3(0.0);
				dx[4 * i + 2] = fvec3(0.0);
				dx[4 * i + 3] = fvec3(0.0);
			}
		}

		for (uint32_t i = 0; i < InnerEdgesize; i++) {
			TempPosition[InnerEdgeList[i]] = TempPosition[InnerEdgeList[i]] + dx[i];
		}
	}

	void FixedProjection(std::vector<fvec3>& TempPosition)
	{

		for (auto ind : fix0) {
			TempPosition[ind].x = RestPositionList[ind].x + 0.5 * edgedist;
			TempPosition[ind].y = RestPositionList[ind].y;
			TempPosition[ind].z = 0.0f;
		}
		for (auto ind : fix1) {
			TempPosition[ind].x = RestPositionList[ind].x - 0.5 * edgedist;
			TempPosition[ind].y = RestPositionList[ind].y;
			TempPosition[ind].z = 0.0f;
		}
	}
};

namespace Physics {

int32_t solver = 4;

void timestep(ClothMesh& CM)
{

	CM.ClearLamda();

	uint32_t NodeSize = CM.vertsize;
	std::vector<fvec3> tempp(NodeSize);

#pragma omp parallel for
	for (uint32_t i = 0; i < NodeSize; i++) {
		fvec3 velocity = CM.VelocityList[i] + dt * fvec3(0.0, -9.8, 0.0);
		tempp[i]       = CM.PositionList[i] + dt * velocity;
	}

	if (solver == 0) {
		for (uint32_t x = 0; x < 10; x++) {
			CM.FemElasticProjectGS(tempp);
			CM.FemBendProjectionGS(tempp);
			CM.FixedProjection(tempp);
		}
	} else if (solver == 1) {
		for (uint32_t x = 0; x < 10; x++) {
			CM.FemElasticProjectJC(tempp, x + 8);
			CM.FemBendProjectionJC(tempp);
			CM.FixedProjection(tempp);
		}
	} else if (solver == 2) {
		for (uint32_t x = 0; x < 10; x++) {
			CM.MassSpringSolverGS(tempp);
			CM.FemBendProjectionGS(tempp);
			CM.FixedProjection(tempp);
		}
	} else if (solver == 3) {
		for (uint32_t x = 0; x < 10; x++) {
			CM.MassSpringSolverJC(tempp);
			CM.FemBendProjectionJC(tempp);
			CM.FixedProjection(tempp);
		}
	} else if (solver == 4) {
		ClearLambdaGPU();
		for (uint32_t x = 0; x < 20; x++) {
			CM.FixedProjection(tempp);
			FemElasticProjectGPU(tempp.data(), lambda, mu);
			if (x % 5 == 0)
				FemBendProjectGPU(tempp.data(), bendCof);
		}
	}

#pragma omp parallel for
	for (uint32_t i = 0; i < NodeSize; i++) {
		CM.VelocityList[i] = (tempp[i] - CM.PositionList[i]) / dt;
		CM.VelocityList[i] = (1.0 - 0.7 * dt) * CM.VelocityList[i];
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

	omp_set_num_threads(20);

	uint32_t cagelist[24] = { 0, 1, 1, 2, 2, 3, 3, 0, 4, 5, 5, 6,
		6, 7, 7, 4, 0, 4, 3, 7, 1, 5, 2, 6 };
	linevertarray cage(8, 24, cagelist);
	cage.setposition(0, -15.0f, -15.0f, -15.0f);
	cage.setposition(1, 15.0f, -15.0f, -15.0f);
	cage.setposition(2, 15.0f, 15.0f, -15.0f);
	cage.setposition(3, -15.0f, 15.0f, -15.0f);
	cage.setposition(4, -15.0f, -15.0f, 15.0f);
	cage.setposition(5, 15.0f, -15.0f, 15.0f);
	cage.setposition(6, 15.0f, 15.0f, 15.0f);
	cage.setposition(7, -15.0f, 15.0f, 15.0f);
	cage.settype(0);
	cage.vboupdate();

	uint32_t floorlist[6] = { 0, 1, 2, 0, 2, 3 };
	trianglevertarray floor(6, 6, floorlist);
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
	floor.vboupdate();

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

	Renderer3D::setclookat(fvec3(0.0, 2.0, 0.0));
	Renderer3D::setLightint(600.0);

	ClothMesh CM0(51, 8.0, fvec3(0.0, 3.0, 0.0));

	//CUDA Init
	MeshInfo mInfo(CM0.vertsize, CM0.RestPositionList.data(), CM0.tilist, CM0.tisize, CM0.edgelist, CM0.edgesize, CM0.InnerEdgeList, CM0.InnerEdgeCList, CM0.InnerEdgesize, CM0.AList.data(), CM0.VList.data(), CM0.mass, dt);
	Init(mInfo);

	//Render object

	std::vector<Renderer3D::drawobject> shadowlist;
	std::vector<Renderer3D::drawobject> edgelist;
	std::vector<Renderer3D::drawobject> renderlist;

	shadowlist.emplace_back(floor);
	shadowlist.emplace_back(CM0.tva);
	shadowlist.emplace_back(CM0.tva2);

	renderlist.emplace_back(floor, simasima0);
	renderlist.emplace_back(cage);
	renderlist.emplace_back(CM0.tva);
	renderlist.emplace_back(CM0.tva2);
	renderlist.emplace_back(CM0.lva);

	//rendering loop

	double ctime = 0.0;
	double vtime = 0.0;

	double prevtime = 0.0;

	while (!Visualizer::Is_Closed()) {
		Renderer3D::Clear();
		//imgui reset
		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();

		static bool is_stop = true;

		static int32_t renderstatus = 0;

		static float mousemovex = 0.5 * 3.14;
		static float mousemovey = 0.5 * 3.14;
		static float cameraL	= 30.0f;

		float camerasinp;
		float cameracosp;
		float camerasint;
		float cameracost;

		static bool nextframe = false;

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
			if (cameraL < 1.0)
				cameraL = 1.0;
			if (cameraL > 80.0)
				cameraL = 80.0;

			camerap.x = cameraL * camerasint * cameracosp;
			camerap.z = cameraL * camerasint * camerasinp;
			camerap.y = cameraL * cameracost;

			ImGui::Text(" x = %.1f", camerap.x);
			ImGui::Text(" y = %.1f", camerap.y);
			ImGui::Text(" z = %.1f", camerap.z);

			ImGui::Text("camera length = %.1f", cameraL);

			if (ImGui::Combo("Render Status", &(renderstatus), "SurfaceEdge\0Surface\0Edge\0triangle\0edge\0\0")) {
				if (renderstatus == 0) {
					renderlist[2].renderswitch = true;
					renderlist[3].renderswitch = true;
					renderlist[4].renderswitch = true;
					CM0.lva.setcolor(0.0, 0.0, 0.0, 1.0);
					CM0.Setdata();
				} else if (renderstatus == 1) {
					renderlist[2].renderswitch = true;
					renderlist[3].renderswitch = true;
					renderlist[4].renderswitch = false;
				} else if (renderstatus == 2) {
					renderlist[2].renderswitch = false;
					renderlist[3].renderswitch = false;
					renderlist[4].renderswitch = true;
					CM0.lva.setcolor(1.0, 1.0, 1.0, 1.0);
					CM0.Setdata();
				} else if (renderstatus == 3) {
					renderlist[2].renderswitch = false;
					renderlist[3].renderswitch = false;
					renderlist[4].renderswitch = false;
				} else if (renderstatus == 4) {
					renderlist[2].renderswitch = false;
					renderlist[3].renderswitch = false;
					renderlist[4].renderswitch = false;
				}
			}

			ImGui::Combo("Solver", &(Physics::solver), "Continuum Gauss-Seidel\0Continuum Jacobi\0MassSpring Gauss-Seidel\0MassSpring Jacobi\0GPU Continuum\0\0");

			ImGui::SliderFloat("mu", &mu, 0.0f, 500.0f, "%.4f", ImGuiSliderFlags_Logarithmic);
			ImGui::SliderFloat("lambda", &lambda, 0.0f, 500.0f, "%.4f", ImGuiSliderFlags_Logarithmic);
			ImGui::SliderFloat("Bending", &bendCof, 0.0f, 1.0f, "%.4f", ImGuiSliderFlags_Logarithmic);
			ImGui::SliderFloat("Spring", &SpringCof, 0.0f, 10.0f, "%.4f", ImGuiSliderFlags_Logarithmic);

			ImGui::SliderFloat("length", &edgedist, -6.0f, 6.0f);

			ImGui::Text("virtualtime = %.1f", vtime);

			if (ImGui::Button("reset")) {

				vtime = 0.0;

				for (uint32_t i = 0; i < CM0.vertsize; i++) {
					CM0.PositionList[i].x = CM0.RestPositionList[i].x;
					CM0.PositionList[i].y = 4.00000;
					CM0.PositionList[i].z = 1.00010 * CM0.RestPositionList[i].y - 4.0;
					CM0.VelocityList[i]   = fvec3(0.0);
				}
			}

			ImGui::End();
		}

		//physics
		ctime = Visualizer::GetTime();
		if ((!is_stop || nextframe) && ctime - prevtime > 0.8 * dt) {

			prevtime = ctime;

			Physics::timestep(CM0);
			vtime += dt;
			nextframe = false;

			CM0.UpdataVertarray();
			CM0.Setdata();
		}

		//renderer set

		if (renderstatus == 3) {
			for (uint32_t i = 0; i < CM0.tisize / 3; i++) {

				fvec3 x0 = CM0.PositionList[CM0.tilist[3 * i + 0]];
				fvec3 x1 = CM0.PositionList[CM0.tilist[3 * i + 1]];
				fvec3 x2 = CM0.PositionList[CM0.tilist[3 * i + 2]];

				fvec3 cm = (x0 + x1 + x2) / 3.0;

				x0 = (x0 - cm) * 0.5 + cm;
				x1 = (x1 - cm) * 0.5 + cm;
				x2 = (x2 - cm) * 0.5 + cm;

				Renderer3D::DrawTriangle(x0, x1, x2);
				Renderer3D::DrawTriangle(x0, x2, x1);
			}
		}
		if (renderstatus == 4) {
			for (uint32_t i = 0; i < CM0.edgesize / 2; i++) {

				fvec3 x0 = CM0.PositionList[CM0.edgelist[2 * i + 0]];
				fvec3 x1 = CM0.PositionList[CM0.edgelist[2 * i + 1]];

				fvec3 cm = (x0 + x1) / 2.0;

				x0 = (x0 - cm) * 0.5 + cm;
				x1 = (x1 - cm) * 0.5 + cm;

				Renderer3D::DrawLine(x0, x1);
			}
		}

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
