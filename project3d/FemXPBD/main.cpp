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
#include "utils/mathfunc/polardecompose.hpp"

#include "utils/fileloader/TetGenLoader.hpp"

#include "utils/meshgenerator/meshgenerator.hpp"

#include "utils/fem/fem.hpp"

#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>

constexpr uint32_t FPS = 60;
constexpr float dt     = 1.0 / FPS;
//glfwSwapInterval(1); なので60FPS以下になる。
//これはモニタやGPUに依るかもしれない。

float mu     = 80;
float lambda = 80;
float rho    = 0.1;

class DeformableMesh {
    public:
	std::vector<fvec3> PositionList;
	std::vector<fvec3> RestPositionList;
	std::vector<fvec3> VelocityList;
	std::vector<fvec3> TempPositionList;

	uint32_t vertsize;
	uint32_t* elementlist;
	uint32_t elementsize; //要素数x4
	uint32_t* tilist;     //三角形数x3
	uint32_t tisize;
	uint32_t* lilist; //辺数x2
	uint32_t lisize;

	linevertarray lva;

	trianglevertarray tva;

	std::vector<float> Lamdalist;

	std::vector<fmat3> AList;
	std::vector<float> VList;
	std::vector<fquaternion> qList;

	int32_t MaterialInd = 0;

	float mass;
	//mass of nodepoint

	DeformableMesh()
	    : lva()
	    , tva()
	{

		fvec3* vertdata;

		LoadTetGentoTetrahedraMesh("../../../resource/Bunny", &vertdata, vertsize, &elementlist, elementsize, &tilist, tisize, &lilist, lisize, fvec3(0.0, -8.0, 0.0), 6.0);
		//LoadTetGentoTetrahedraMesh("../../../resource/Dragon", &vertdata, vertsize, &elementlist, elementsize, &tilist, tisize, &lilist, lisize, fvec3(0.0, -6.0, 0.0), 18.0);
		//CubeTetrahedra(6, 10.0, &vertdata, vertsize, &elementlist, elementsize, &tilist, tisize, &lilist, lisize, fvec3(0.0, 0.0, 0.0));

		tva.resetvertarray(vertsize, tisize, tilist);
		lva.resetvertarray(vertsize, lisize, lilist);
		tva.settype(1);
		lva.setcolor(0.1, 0.1, 0.1, 1.0);

		RestPositionList.resize(vertsize);
		PositionList.resize(vertsize);
		VelocityList.resize(vertsize);
		TempPositionList.resize(vertsize);

		for (uint32_t i = 0; i < vertsize; i++) {
			RestPositionList[i] = vertdata[i];
			PositionList[i]	    = vertdata[i];
			VelocityList[i]	    = fvec3(0.0);
		}

		AList.resize(elementsize / 4);
		VList.resize(elementsize / 4);
		qList.resize(elementsize / 4);
		Lamdalist.resize(elementsize / 4);

		for (uint32_t i = 0; i < elementsize / 4; i++) {
			fvec3 X0 = RestPositionList[elementlist[4 * i + 0]];
			fvec3 X1 = RestPositionList[elementlist[4 * i + 1]];
			fvec3 X2 = RestPositionList[elementlist[4 * i + 2]];
			fvec3 X3 = RestPositionList[elementlist[4 * i + 3]];

			AList[i] = fmat3(X1 - X0, X2 - X0, X3 - X0).inverse();
			VList[i] = fvec3::STP(X1 - X0, X2 - X0, X3 - X0) / 6.0;
			qList[i] = fquaternion(0.0, 0.0, 0.0, 1.0);
		}

		mass = 0.00001;

		this->UpdataVertarray();

		this->Setdata();
	}

	void UpdataVertarray()
	{

		//normal

		std::vector<fvec3> NormalSet(vertsize);

#pragma omp parallel for
		for (auto& x : NormalSet)
			x = fvec3(0.0);

		for (uint32_t i = 0; i < tisize / 3; i++) {
			fvec3 v0 = PositionList[tilist[3 * i + 0]];
			fvec3 v1 = PositionList[tilist[3 * i + 1]];
			fvec3 v2 = PositionList[tilist[3 * i + 2]];

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

			tva[i].position[0] = PositionList[i].x;
			tva[i].position[1] = PositionList[i].y;
			tva[i].position[2] = PositionList[i].z;

			lva[i].position[0] = PositionList[i].x + 0.001 * NormalSet[i].x;
			lva[i].position[1] = PositionList[i].y + 0.001 * NormalSet[i].y;
			lva[i].position[2] = PositionList[i].z + 0.001 * NormalSet[i].z;
		}
	}

	void Setdata()
	{
		tva.vboupdate();
		lva.vboupdate();
	}

	void ClearLamda()
	{
#pragma omp parallel for
		for (auto& x : Lamdalist) {
			x = 0.0;
		}
	}

	void FemProjectGS()
	{
		for (uint32_t i = 0; i < elementsize / 4; i++) {
			const fvec3& X0 = RestPositionList[elementlist[4 * i + 0]];
			const fvec3& X1 = RestPositionList[elementlist[4 * i + 1]];
			const fvec3& X2 = RestPositionList[elementlist[4 * i + 2]];
			const fvec3& X3 = RestPositionList[elementlist[4 * i + 3]];

			const fvec3& x0 = TempPositionList[elementlist[4 * i + 0]];
			const fvec3& x1 = TempPositionList[elementlist[4 * i + 1]];
			const fvec3& x2 = TempPositionList[elementlist[4 * i + 2]];
			const fvec3& x3 = TempPositionList[elementlist[4 * i + 3]];

			const fmat3& A = AList[i];
			const fmat3 F  = mat3(x1 - x0, x2 - x0, x3 - x0) * A;
			const fmat3 E  = 0.5 * (F.transpose() * F - fmat3::indentity());

			const float& V = VList[i];

			float W;
			fvec3 dx0, dx1, dx2, dx3;
			bool ValidEnergy;

			if (MaterialInd == 0) {
				ValidEnergy = FemElasticDxStVenant(F, E, A, V, lambda, mu, W, dx0, dx1, dx2, dx3);
			} else if (MaterialInd == 1) {
				ValidEnergy = FemElasticDxNeoHookean(F, E, A, V, lambda, mu, W, dx0, dx1, dx2, dx3);
			} else if (MaterialInd == 2) {
				ValidEnergy = FemElasticDxCoRotational(F, E, A, qList[i], V, lambda, mu, W, dx0, dx1, dx2, dx3);
			}

			if (ValidEnergy && W > 0.0001) {
				float C = std::sqrt(2.0 * W);

				//fmat3 BAt = B * A.transpose();
				//fvec3 dC1 = (1.0 / C) * fvec3(BAt.m[0], BAt.m[3], BAt.m[6]);
				//fvec3 dC2 = (1.0 / C) * fvec3(BAt.m[1], BAt.m[4], BAt.m[7]);
				//fvec3 dC3 = (1.0 / C) * fvec3(BAt.m[2], BAt.m[5], BAt.m[8]);
				//fvec3 dC0 = -(dC1 + dC2 + dC3);

				fvec3 dC1 = (1.0 / C) * dx1;
				fvec3 dC2 = (1.0 / C) * dx2;
				fvec3 dC3 = (1.0 / C) * dx3;
				fvec3 dC0 = (1.0 / C) * dx0;

				float dtdtdlambda = (-C - Lamdalist[i]) / ((dC0.sqlength() + dC1.sqlength() + dC2.sqlength() + dC3.sqlength()) / mass + 1.0 / (dt * dt));

				TempPositionList[elementlist[4 * i + 0]] = TempPositionList[elementlist[4 * i + 0]] + dtdtdlambda * (1.0 / mass) * dC0;
				TempPositionList[elementlist[4 * i + 1]] = TempPositionList[elementlist[4 * i + 1]] + dtdtdlambda * (1.0 / mass) * dC1;
				TempPositionList[elementlist[4 * i + 2]] = TempPositionList[elementlist[4 * i + 2]] + dtdtdlambda * (1.0 / mass) * dC2;
				TempPositionList[elementlist[4 * i + 3]] = TempPositionList[elementlist[4 * i + 3]] + dtdtdlambda * (1.0 / mass) * dC3;

				Lamdalist[i] += dtdtdlambda / (dt * dt);
			}

			/*
			if (F.det() < 0.00001 && lambda > 0.00001 && mu > 0.00001 && MaterialInd != 2) {
				fquaternion q = ExtractRotation(F, 1, qList[i]);
				qList[i]      = q;
				fmat3 R	      = q.qtomat();

				fmat3 S	  = R.transpose() * F;
				float trS = S.trace();

				W	= 0.5 * (F.sqlength() - 2.0 * trS + 3);
				fmat3 B = (F - R);

				float C = std::sqrt(2.0 * W);

				fmat3 BAt = B * A.transpose();

				fvec3 dC1 = fvec3(BAt.m[0], BAt.m[3], BAt.m[6]);
				fvec3 dC2 = fvec3(BAt.m[1], BAt.m[4], BAt.m[7]);
				fvec3 dC3 = fvec3(BAt.m[2], BAt.m[5], BAt.m[8]);
				fvec3 dC0 = -(dC1 + dC2 + dC3);

				float dtdtdlambda = -C / ((dC0.sqlength() + dC1.sqlength() + dC2.sqlength() + dC3.sqlength()));

				TempPositionList[elementlist[4 * i + 0]] = TempPositionList[elementlist[4 * i + 0]] + dtdtdlambda * dC0;
				TempPositionList[elementlist[4 * i + 1]] = TempPositionList[elementlist[4 * i + 1]] + dtdtdlambda * dC1;
				TempPositionList[elementlist[4 * i + 2]] = TempPositionList[elementlist[4 * i + 2]] + dtdtdlambda * dC2;
				TempPositionList[elementlist[4 * i + 3]] = TempPositionList[elementlist[4 * i + 3]] + dtdtdlambda * dC3;
			}
			*/
		}
	}

	void FemProjectJC()
	{

		std::vector<fvec3> dx(elementsize);

#pragma omp parallel for
		for (uint32_t i = 0; i < elementsize / 4; i++) {

			const fvec3& X0 = RestPositionList[elementlist[4 * i + 0]];
			const fvec3& X1 = RestPositionList[elementlist[4 * i + 1]];
			const fvec3& X2 = RestPositionList[elementlist[4 * i + 2]];
			const fvec3& X3 = RestPositionList[elementlist[4 * i + 3]];

			const fvec3& x0 = TempPositionList[elementlist[4 * i + 0]];
			const fvec3& x1 = TempPositionList[elementlist[4 * i + 1]];
			const fvec3& x2 = TempPositionList[elementlist[4 * i + 2]];
			const fvec3& x3 = TempPositionList[elementlist[4 * i + 3]];

			const fmat3& A = AList[i];
			const fmat3 F  = mat3(x1 - x0, x2 - x0, x3 - x0) * A;
			const fmat3 E  = 0.5 * (F.transpose() * F - fmat3::indentity());

			const float& V = VList[i];

			float W;
			fvec3 dx0, dx1, dx2, dx3;
			bool ValidEnergy;

			/*
			if (MaterialInd == 0) {
				W = 5.0 * V * (mu * E.sqlength() + 0.5 * lambda * E.trace() * E.trace());
				B = 5.0 * V * (2 * mu * F * E + lambda * E.trace() * F);
			} else if (MaterialInd == 1) {
				float J	   = F.det();
				float logJ = std::log(J);
				if (J < 0.00001)
					logJ = 0.0;
				W = V * (0.5 * mu * (F.sqlength() - 3) - mu * logJ + 0.5 * lambda * logJ * logJ);

				if (std::abs(F.det()) < 0.00001)
					W = 0.0;
				B = V * (mu * F - mu * (F.inverse()).transpose() + lambda * logJ * (F.inverse()).transpose());
			} else if (MaterialInd == 2) {
				fquaternion q = ExtractRotation(F, 3, qList[i]);
				qList[i]      = q;
				fmat3 R	      = q.qtomat();
				fmat3 S	      = R.transpose() * F;
				float trS     = S.trace();
				W	      = V * (0.5 * mu * (F.sqlength() - 2.0 * trS + 3) + 0.5 * lambda * (trS * trS - 6.0 * trS + 9.0));
				B	      = V * (mu * (F - R) + lambda * (trS - 3) * R);
			}
			*/
			if (MaterialInd == 0) {
				ValidEnergy = FemElasticDxStVenant(F, E, A, V, lambda, mu, W, dx0, dx1, dx2, dx3);
			} else if (MaterialInd == 1) {
				ValidEnergy = FemElasticDxNeoHookean(F, E, A, V, lambda, mu, W, dx0, dx1, dx2, dx3);

				//float J	   = F.det();
				//float logJ = std::log(J);
				//if (J < 0.00001)
				//	logJ = 0.0;
				//W = V * (0.5 * mu * (F.sqlength() - 3) - mu * logJ + 0.5 * lambda * logJ * logJ);

				//if (std::abs(F.det()) < 0.00001)
				//	W = 0.0;
				//B = V * (mu * F - mu * (F.inverse()).transpose() + lambda * logJ * (F.inverse()).transpose());
			} else if (MaterialInd == 2) {
				//fquaternion q = ExtractRotation(F, 3, qList[i]);
				//qList[i]      = q;
				ValidEnergy = FemElasticDxCoRotational(F, E, A, qList[i], V, lambda, mu, W, dx0, dx1, dx2, dx3);

				//fmat3 R	  = q.qtomat();
				//fmat3 S	  = R.transpose() * F;
				//float trS = S.trace();
				//W	  = V * (0.5 * mu * (F.sqlength() - 2.0 * trS + 3) + 0.5 * lambda * (trS * trS - 6.0 * trS + 9.0));
				//B	  = V * (mu * (F - R) + lambda * (trS - 3) * R);
			}

			//std::cout << W << std::endl;

			if (ValidEnergy && W > 0.0001) {
				float C = std::sqrt(2.0 * W);

				//fmat3 BAt = B * A.transpose();
				//fvec3 dC1 = (1.0 / C) * fvec3(BAt.m[0], BAt.m[3], BAt.m[6]);
				//fvec3 dC2 = (1.0 / C) * fvec3(BAt.m[1], BAt.m[4], BAt.m[7]);
				//fvec3 dC3 = (1.0 / C) * fvec3(BAt.m[2], BAt.m[5], BAt.m[8]);
				//fvec3 dC0 = -(dC1 + dC2 + dC3);

				fvec3 dC1 = (1.0 / C) * dx1;
				fvec3 dC2 = (1.0 / C) * dx2;
				fvec3 dC3 = (1.0 / C) * dx3;
				fvec3 dC0 = (1.0 / C) * dx0;

				float dtdtdlambda = (-C - Lamdalist[i]) / ((dC0.sqlength() + dC1.sqlength() + dC2.sqlength() + dC3.sqlength()) / mass + 1.0 / (dt * dt));
				dtdtdlambda *= 0.15;

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

			/*
			if (F.det() < 0.00001 && lambda > 0.00001 && mu > 0.00001) {
				fquaternion q = ExtractRotation(F, 1, qList[i]);
				qList[i]      = q;
				fmat3 R	      = q.qtomat();

				fmat3 S	  = R.transpose() * F;
				float trS = S.trace();

				W	= 0.5 * (F.sqlength() - 2.0 * trS + 3);
				fmat3 B = (F - R);

				float C = std::sqrt(2.0 * W);

				fmat3 BAt = B * A.transpose();

				fvec3 dC1 = fvec3(BAt.m[0], BAt.m[3], BAt.m[6]);
				fvec3 dC2 = fvec3(BAt.m[1], BAt.m[4], BAt.m[7]);
				fvec3 dC3 = fvec3(BAt.m[2], BAt.m[5], BAt.m[8]);
				fvec3 dC0 = -(dC1 + dC2 + dC3);

				float dtdtdlambda = -C / ((dC0.sqlength() + dC1.sqlength() + dC2.sqlength() + dC3.sqlength()));
				dtdtdlambda *= 0.35;

				dx[4 * i + 0] = dtdtdlambda * dC0;
				dx[4 * i + 1] = dtdtdlambda * dC1;
				dx[4 * i + 2] = dtdtdlambda * dC2;
				dx[4 * i + 3] = dtdtdlambda * dC3;
			}
			*/
		}

#pragma omp parallel for
		for (uint32_t i = 0; i < elementsize; i++) {
			TempPositionList[elementlist[i]] = TempPositionList[elementlist[i]] + dx[i];
		}
	}

	void FemFixInversionGS()
	{
		for (uint32_t i = 0; i < elementsize / 4; i++) {
			const fvec3& X0 = RestPositionList[elementlist[4 * i + 0]];
			const fvec3& X1 = RestPositionList[elementlist[4 * i + 1]];
			const fvec3& X2 = RestPositionList[elementlist[4 * i + 2]];
			const fvec3& X3 = RestPositionList[elementlist[4 * i + 3]];

			const fvec3& x0 = TempPositionList[elementlist[4 * i + 0]];
			const fvec3& x1 = TempPositionList[elementlist[4 * i + 1]];
			const fvec3& x2 = TempPositionList[elementlist[4 * i + 2]];
			const fvec3& x3 = TempPositionList[elementlist[4 * i + 3]];

			const fmat3& A = AList[i];
			const fmat3 F  = mat3(x1 - x0, x2 - x0, x3 - x0) * A;
			const fmat3 E  = 0.5 * (F.transpose() * F - fmat3::indentity());

			if (F.det() < 0.00001 && lambda > 0.00001 && mu > 0.00001 && MaterialInd != 2) {
				fquaternion q = ExtractRotation(F, 3, qList[i]);
				qList[i]      = q;
				fmat3 R	      = q.qtomat();

				fmat3 S	  = R.transpose() * F;
				float trS = S.trace();

				float W = 0.5 * (F.sqlength() - 2.0 * trS + 3);
				fmat3 B = (F - R);

				float C = std::sqrt(2.0 * W);

				fmat3 BAt = B * A.transpose();

				fvec3 dC1 = fvec3(BAt.m[0], BAt.m[3], BAt.m[6]);
				fvec3 dC2 = fvec3(BAt.m[1], BAt.m[4], BAt.m[7]);
				fvec3 dC3 = fvec3(BAt.m[2], BAt.m[5], BAt.m[8]);
				fvec3 dC0 = -(dC1 + dC2 + dC3);

				float dtdtdlambda = -C / ((dC0.sqlength() + dC1.sqlength() + dC2.sqlength() + dC3.sqlength()));

				TempPositionList[elementlist[4 * i + 0]] = TempPositionList[elementlist[4 * i + 0]] + dtdtdlambda * dC0;
				TempPositionList[elementlist[4 * i + 1]] = TempPositionList[elementlist[4 * i + 1]] + dtdtdlambda * dC1;
				TempPositionList[elementlist[4 * i + 2]] = TempPositionList[elementlist[4 * i + 2]] + dtdtdlambda * dC2;
				TempPositionList[elementlist[4 * i + 3]] = TempPositionList[elementlist[4 * i + 3]] + dtdtdlambda * dC3;
			}
		}
	}

	void FemFixInversionJC()
	{
		std::vector<fvec3> dx(elementsize);

#pragma omp parallel for
		for (uint32_t i = 0; i < elementsize / 4; i++) {

			const fvec3& X0 = RestPositionList[elementlist[4 * i + 0]];
			const fvec3& X1 = RestPositionList[elementlist[4 * i + 1]];
			const fvec3& X2 = RestPositionList[elementlist[4 * i + 2]];
			const fvec3& X3 = RestPositionList[elementlist[4 * i + 3]];

			const fvec3& x0 = TempPositionList[elementlist[4 * i + 0]];
			const fvec3& x1 = TempPositionList[elementlist[4 * i + 1]];
			const fvec3& x2 = TempPositionList[elementlist[4 * i + 2]];
			const fvec3& x3 = TempPositionList[elementlist[4 * i + 3]];

			const fmat3& A = AList[i];
			const fmat3 F  = mat3(x1 - x0, x2 - x0, x3 - x0) * A;
			const fmat3 E  = 0.5 * (F.transpose() * F - fmat3::indentity());

			if (F.det() < 0.00001 && lambda > 0.00001 && mu > 0.00001) {
				fquaternion q = ExtractRotation(F, 3, qList[i]);
				qList[i]      = q;
				fmat3 R	      = q.qtomat();

				fmat3 S	  = R.transpose() * F;
				float trS = S.trace();

				float W = 0.5 * (F.sqlength() - 2.0 * trS + 3);
				fmat3 B = (F - R);

				float C = std::sqrt(2.0 * W);

				fmat3 BAt = B * A.transpose();

				fvec3 dC1 = fvec3(BAt.m[0], BAt.m[3], BAt.m[6]);
				fvec3 dC2 = fvec3(BAt.m[1], BAt.m[4], BAt.m[7]);
				fvec3 dC3 = fvec3(BAt.m[2], BAt.m[5], BAt.m[8]);
				fvec3 dC0 = -(dC1 + dC2 + dC3);

				float dtdtdlambda = -C / ((dC0.sqlength() + dC1.sqlength() + dC2.sqlength() + dC3.sqlength()));
				dtdtdlambda *= 0.35;

				dx[4 * i + 0] = dtdtdlambda * dC0;
				dx[4 * i + 1] = dtdtdlambda * dC1;
				dx[4 * i + 2] = dtdtdlambda * dC2;
				dx[4 * i + 3] = dtdtdlambda * dC3;
			}
		}

#pragma omp parallel for
		for (uint32_t i = 0; i < elementsize; i++) {
			TempPositionList[elementlist[i]] = TempPositionList[elementlist[i]] + dx[i];
		}
	}

	void FixedProjection()
	{

		//for (uint32_t i = N * N * (N - 1); i < N * N * N; i++) {
		//	TempPosition[i] = RestPositionList[i];
		//}

#pragma omp parallel for
		for (uint32_t i = 0; i < vertsize; i++) {
			fvec3& position = TempPositionList[i];
			if (position.y < -12.0)
				position.y = -11.99;
			if (position.y > 12.0)
				position.y = 11.99;

			if (position.x < -12.0)
				position.x = -11.99;
			if (position.x > 12.0)
				position.x = 11.99;

			if (position.z < -12.0)
				position.z = -11.99;
			if (position.z > 12.0)
				position.z = 11.99;
		}

		/*
		if (rightwall) {

#pragma omp parallel for
			for (uint32_t y = 0; y < N; y++) {
				for (uint32_t z = 0; z < N; z++) {
					uint32_t Ind	    = N * N * y + N * z + N - 1;
					TempPosition[Ind]   = fvec3(rightrot, 0.0, 0.0).rotation() * (RestPositionList[Ind] - bias) + bias;
					TempPosition[Ind].x = bias.x + rightedge;
				}
			}

#pragma omp parallel for
			for (uint32_t y = 0; y < N; y++) {
				for (uint32_t z = 0; z < N; z++) {
					uint32_t Ind	  = N * N * y + N * z;
					TempPosition[Ind] = RestPositionList[Ind];
				}
			}
		}
	*/
	}
};

namespace Physics {

int32_t solver = 0;

void timestep(DeformableMesh& CM)
{

	CM.ClearLamda();

	uint32_t NodeSize = CM.vertsize;

#pragma omp parallel for
	for (uint32_t i = 0; i < NodeSize; i++) {
		fvec3 velocity	       = CM.VelocityList[i] + dt * fvec3(0.0, -9.8, 0.0);
		CM.TempPositionList[i] = CM.PositionList[i] + dt * velocity;
	}

	if (solver == 0)
		for (uint32_t x = 0; x < 5; x++) {
			CM.FemProjectGS();
			CM.FemFixInversionGS();
			CM.FixedProjection();
		}
	else if (solver == 1)
		for (uint32_t x = 0; x < 15; x++) {
			CM.FemProjectJC();
			CM.FemFixInversionJC();
			CM.FixedProjection();
		}

#pragma omp parallel for
	for (uint32_t i = 0; i < NodeSize; i++) {
		CM.VelocityList[i] = (CM.TempPositionList[i] - CM.PositionList[i]) / dt;
		CM.VelocityList[i] = (1.0 - 0.7 * dt) * CM.VelocityList[i];
		CM.PositionList[i] = CM.TempPositionList[i];
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

	fvec3 camerap(0.0, 15.0, 20.0);

	//myphysics.rbodyvec.reserve(10);

	DeformableMesh CM0;

	//init

	std::vector<Renderer3D::drawobject> shadowlist;
	std::vector<Renderer3D::drawobject> edgelist;
	std::vector<Renderer3D::drawobject> renderlist;

	shadowlist.emplace_back(floor);
	shadowlist.emplace_back(CM0.tva);

	renderlist.emplace_back(floor, simasima0);
	renderlist.emplace_back(cage);
	renderlist.emplace_back(CM0.tva);
	renderlist.emplace_back(CM0.lva);

	Renderer3D::setLightint(500.0);

	//random

	std::random_device seedgen;
	std::mt19937 engine(seedgen());
	std::uniform_real_distribution<> udist(-50.0, 50.0);

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
		ImGui::SetNextWindowSize(ImVec2(300, 450), ImGuiCond_FirstUseEver);

		static bool is_stop = true;

		static float mousemovex = 0.5 * 3.14;
		static float mousemovey = 0.5 * 3.14;
		static float cameraL	= 30.0f;

		float camerasinp;
		float cameracosp;
		float camerasint;
		float cameracost;

		static bool nextframe = false;

		static int32_t renderstatus = 0;

		static float realdt = 0.1;

		//imgui
		{

			ImGui::Begin("Deformable Body");

			ImGui::Text("Rendering FPS : %.1f", ImGui::GetIO().Framerate);
			ImGui::Text("Simulation FPS : %.1f", 1.0 / realdt);

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

			if (ImGui::Combo("Render Status", &(renderstatus), "Surface\0Volume\0\0")) {
				if (renderstatus == 0) {
					renderlist[2].renderswitch = true;
					renderlist[3].renderswitch = true;
				} else if (renderstatus == 1) {
					renderlist[2].renderswitch = false;
					renderlist[3].renderswitch = false;
				}
			}

			ImGui::Combo("Solver", &(Physics::solver), "Gauss-Seidel\0Jacobi\0\0");

			ImGui::Combo("Material", &(CM0.MaterialInd), "Saint Venant-Kirchhoff\0Neo-Hookean\0Co-Rotational\0\0");
			ImGui::SliderFloat("mu", &mu, 0.0f, 150.0f, "%.4f", ImGuiSliderFlags_Logarithmic);
			ImGui::SliderFloat("lambda", &lambda, 0.0f, 150.0f, "%.4f", ImGuiSliderFlags_Logarithmic);

			ImGui::Text("realtime = %.1f", ctime);
			ImGui::Text("virtualtime = %.1f", vtime);

			if (ImGui::Button("reset")) {
				vtime = 0.0;

				for (uint32_t i = 0; i < CM0.vertsize; i++) {
					CM0.PositionList[i] = CM0.RestPositionList[i];
					CM0.VelocityList[i] = fvec3(0.0);

					for (auto& q : CM0.qList)
						q = fquaternion(0.0, 0.0, 0.0, 1.0);
				}
				CM0.UpdataVertarray();
				CM0.Setdata();
			}

			if (ImGui::Button("zero squash")) {
				vtime = 0.0;

				for (uint32_t i = 0; i < CM0.vertsize; i++) {
					CM0.PositionList[i] = fvec3(0.0);

					CM0.VelocityList[i] = fvec3(0.0);

					for (auto& q : CM0.qList)
						q = fquaternion(0.0, 0.0, 0.0, 1.0);
				}
				CM0.UpdataVertarray();
				CM0.Setdata();
			}

			if (ImGui::Button("random")) {
				vtime = 0.0;

				for (uint32_t i = 0; i < CM0.vertsize; i++) {
					CM0.PositionList[i].x = std::clamp(udist(engine), -12.0, 12.0);
					CM0.PositionList[i].y = std::clamp(udist(engine), -12.0, 12.0);
					CM0.PositionList[i].z = std::clamp(udist(engine), -12.0, 12.0);

					CM0.VelocityList[i] = fvec3(0.0);

					for (auto& q : CM0.qList)
						q = fquaternion(0.0, 0.0, 0.0, 1.0);
				}
				CM0.UpdataVertarray();
				CM0.Setdata();
			}

			ImGui::End();
		}

		//physics
		ctime = Visualizer::GetTime();

		if ((!is_stop || nextframe) && ctime - prevtime > 0.8 * dt) {

			realdt	 = ctime - prevtime;
			prevtime = ctime;

			Physics::timestep(CM0);
			vtime += dt;
			nextframe = false;

			CM0.UpdataVertarray();
			CM0.Setdata();
		}

		if (renderstatus == 1) {
			for (uint32_t i = 0; i < CM0.elementsize / 4; i++) {
				fvec3 x0 = CM0.PositionList[CM0.elementlist[4 * i + 0]];
				fvec3 x1 = CM0.PositionList[CM0.elementlist[4 * i + 1]];
				fvec3 x2 = CM0.PositionList[CM0.elementlist[4 * i + 2]];
				fvec3 x3 = CM0.PositionList[CM0.elementlist[4 * i + 3]];

				fvec3 cm = (x0 + x1 + x2 + x3) / 4.0;

				x0 = (x0 - cm) * 0.5 + cm;
				x1 = (x1 - cm) * 0.5 + cm;
				x2 = (x2 - cm) * 0.5 + cm;
				x3 = (x3 - cm) * 0.5 + cm;

				Renderer3D::DrawTetrahedron(x0, x1, x2, x3);
			}
		}

		//for (uint32_t i = 0; i < CM0.vertsize; i++) {
		//	std::cout << CM0.PositionList[i] << std::endl;
		//}

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
