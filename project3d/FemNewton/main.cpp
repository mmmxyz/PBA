#include <iostream>
#include <cmath>
#include <cstdint>
#include <vector>
#include <deque>
#include <algorithm>
#include <random>
#include <omp.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

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
	std::vector<fvec3> extPositionList;

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
	std::vector<bool> ValidationList;

	int32_t MaterialInd = 1;

	float mass;
	//mass of nodepoint

	Eigen::SparseMatrix<float> GlobalHessian;
	Eigen::MatrixXf GlobalB;
	Eigen::MatrixXf DecentDirecton;

	Eigen::BiCGSTAB<Eigen::SparseMatrix<float>> Solver;

	struct HessianGradient {
		Eigen::Matrix3f Hessians[16];
		Eigen::Vector3f Gradient[4];
	};
	std::vector<HessianGradient> ElementHGList;

	DeformableMesh()
	    : lva()
	    , tva()
	{

		fvec3* vertdata;

		LoadTetGentoTetrahedraMesh("../../../resource/LightBunny", &vertdata, vertsize, &elementlist, elementsize, &tilist, tisize, &lilist, lisize, fvec3(0.0, -8.0, 0.0), 6.0);
		//LoadTetGentoTetrahedraMesh("../../../resource/Bunny", &vertdata, vertsize, &elementlist, elementsize, &tilist, tisize, &lilist, lisize, fvec3(0.0, -8.0, 0.0), 6.0);
		//LoadTetGentoTetrahedraMesh("../../../resource/Dragon", &vertdata, vertsize, &elementlist, elementsize, &tilist, tisize, &lilist, lisize, fvec3(0.0, -8.0, 0.0), 18.0);
		//CubeTetrahedra(5, 10.0, &vertdata, vertsize, &elementlist, elementsize, &tilist, tisize, &lilist, lisize, fvec3(0.0, 0.0, 0.0));

		tva.resetvertarray(vertsize, tisize, tilist);
		lva.resetvertarray(vertsize, lisize, lilist);
		tva.settype(1);
		lva.setcolor(0.1, 0.1, 0.1, 1.0);

		RestPositionList.resize(vertsize);
		PositionList.resize(vertsize);
		VelocityList.resize(vertsize);
		TempPositionList.resize(vertsize);
		extPositionList.resize(vertsize);

		GlobalHessian.resize(3 * vertsize, 3 * vertsize);
		GlobalB.resize(3 * vertsize, 1);
		DecentDirecton.resize(3 * vertsize, 1);

		//init Hessian
		std::vector<Eigen::Triplet<float>> tripletList;
		for (uint32_t x = 0; x < elementsize / 4; x++) {

			for (uint32_t k = 0; k < 4; k++)
				for (uint32_t l = 0; l < 4; l++)
					for (uint32_t i = 0; i < 3; i++)
						for (uint32_t j = 0; j < 3; j++)
							tripletList.emplace_back(3 * elementlist[4 * x + k] + i, 3 * elementlist[4 * x + l] + j, 0.0);
		}
		GlobalHessian.setFromTriplets(tripletList.begin(), tripletList.end());

		for (uint32_t i = 0; i < vertsize; i++) {
			RestPositionList[i] = vertdata[i];
			PositionList[i]	    = vertdata[i];
			VelocityList[i]	    = fvec3(0.0);
		}

		AList.resize(elementsize / 4);
		VList.resize(elementsize / 4);
		qList.resize(elementsize / 4);
		Lamdalist.resize(elementsize / 4);
		ValidationList.resize(elementsize / 4);

		ElementHGList.resize(elementsize / 4);

		for (uint32_t i = 0; i < elementsize / 4; i++) {
			fvec3 X0 = RestPositionList[elementlist[4 * i + 0]];
			fvec3 X1 = RestPositionList[elementlist[4 * i + 1]];
			fvec3 X2 = RestPositionList[elementlist[4 * i + 2]];
			fvec3 X3 = RestPositionList[elementlist[4 * i + 3]];

			AList[i]	  = fmat3(X1 - X0, X2 - X0, X3 - X0).inverse();
			VList[i]	  = fvec3::STP(X1 - X0, X2 - X0, X3 - X0) / 6.0;
			qList[i]	  = fquaternion(0.0, 0.0, 0.0, 1.0);
			ValidationList[i] = true;

			ElementHGList[i].Hessians[3 * 0 + 0] = Eigen::MatrixXf::Zero(3, 3);
			ElementHGList[i].Hessians[3 * 0 + 1] = Eigen::MatrixXf::Zero(3, 3);
			ElementHGList[i].Hessians[3 * 0 + 2] = Eigen::MatrixXf::Zero(3, 3);
			ElementHGList[i].Hessians[3 * 1 + 0] = Eigen::MatrixXf::Zero(3, 3);
			ElementHGList[i].Hessians[3 * 1 + 1] = Eigen::MatrixXf::Zero(3, 3);
			ElementHGList[i].Hessians[3 * 1 + 2] = Eigen::MatrixXf::Zero(3, 3);
			ElementHGList[i].Hessians[3 * 2 + 0] = Eigen::MatrixXf::Zero(3, 3);
			ElementHGList[i].Hessians[3 * 2 + 1] = Eigen::MatrixXf::Zero(3, 3);
			ElementHGList[i].Hessians[3 * 2 + 2] = Eigen::MatrixXf::Zero(3, 3);
			ElementHGList[i].Hessians[3 * 3 + 0] = Eigen::MatrixXf::Zero(3, 3);
			ElementHGList[i].Hessians[3 * 3 + 1] = Eigen::MatrixXf::Zero(3, 3);
			ElementHGList[i].Hessians[3 * 3 + 2] = Eigen::MatrixXf::Zero(3, 3);
			ElementHGList[i].Gradient[0]	     = Eigen::MatrixXf::Zero(3, 1);
			ElementHGList[i].Gradient[1]	     = Eigen::MatrixXf::Zero(3, 1);
			ElementHGList[i].Gradient[2]	     = Eigen::MatrixXf::Zero(3, 1);
			ElementHGList[i].Gradient[3]	     = Eigen::MatrixXf::Zero(3, 1);
		}

		//mass = 0.000001;
		mass = 0.05;
		//mass = 0.5;

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

	void ClearHessian()
	{

#pragma omp parallel for
		for (uint32_t i = 0; i < GlobalHessian.outerSize(); i++)
			for (Eigen::SparseMatrix<float>::InnerIterator itr(GlobalHessian, i); itr; ++itr)
				itr.valueRef() = 0.0;

		//for (uint32_t i = 0; i < elementsize / 4; i++)
		//	ValidationList[i] = true;
	}

	//bool EvaluateLocalEnergy(const uint32_t& index, const Eigen::Matrix<float, 12, 1>& Dx, const float& alpha, float& Energy)
	//{

	//	Energy = 0.0;

	//	const uint32_t& Ind0 = elementlist[4 * index + 0];
	//	const uint32_t& Ind1 = elementlist[4 * index + 1];
	//	const uint32_t& Ind2 = elementlist[4 * index + 2];
	//	const uint32_t& Ind3 = elementlist[4 * index + 3];

	//	const fvec3& X0 = RestPositionList[Ind0];
	//	const fvec3& X1 = RestPositionList[Ind1];
	//	const fvec3& X2 = RestPositionList[Ind2];
	//	const fvec3& X3 = RestPositionList[Ind3];

	//	const fvec3 x0 = TempPositionList[Ind0] + alpha * fvec3(Dx(3 * 0 + 0), Dx(3 * 0 + 1), Dx(3 * 0 + 2));
	//	const fvec3 x1 = TempPositionList[Ind1] + alpha * fvec3(Dx(3 * 1 + 0), Dx(3 * 1 + 1), Dx(3 * 1 + 2));
	//	const fvec3 x2 = TempPositionList[Ind2] + alpha * fvec3(Dx(3 * 2 + 0), Dx(3 * 2 + 1), Dx(3 * 2 + 2));
	//	const fvec3 x3 = TempPositionList[Ind3] + alpha * fvec3(Dx(3 * 3 + 0), Dx(3 * 3 + 1), Dx(3 * 3 + 2));

	//	const fvec3& s0 = extPositionList[Ind0];
	//	const fvec3& s1 = extPositionList[Ind1];
	//	const fvec3& s2 = extPositionList[Ind2];
	//	const fvec3& s3 = extPositionList[Ind3];

	//	const fmat3& A = AList[index];
	//	const fmat3 F  = mat3(x1 - x0, x2 - x0, x3 - x0) * A;
	//	const fmat3 E  = 0.5 * (F.transpose() * F - fmat3::indentity());

	//	const float& V = VList[index];

	//	float W;
	//	fvec3 dx0, dx1, dx2, dx3;

	//	fmat3 AFinv, AAt;
	//	float logJ;
	//	bool ValidEnergy;
	//	ValidEnergy = FemElasticHessianNeoHookean(F, E, A, V, lambda, mu, W, logJ, dx0, dx1, dx2, dx3, AFinv, AAt);

	//	if (!ValidEnergy)
	//		return false;

	//	Energy += (mass / (2.0 * dt * dt)) * ((x0 - s0).sqlength());
	//	Energy += (mass / (2.0 * dt * dt)) * ((x1 - s1).sqlength());
	//	Energy += (mass / (2.0 * dt * dt)) * ((x2 - s2).sqlength());
	//	Energy += (mass / (2.0 * dt * dt)) * ((x3 - s3).sqlength());

	//	Energy += W;

	//	return true;
	//}

	//void FemProjectLocalNewton()
	//{
	//	for (uint32_t i = 0; i < elementsize / 4; i++) {
	//		const uint32_t Ind0 = elementlist[4 * i + 0];
	//		const uint32_t Ind1 = elementlist[4 * i + 1];
	//		const uint32_t Ind2 = elementlist[4 * i + 2];
	//		const uint32_t Ind3 = elementlist[4 * i + 3];

	//		const fvec3& X0 = RestPositionList[Ind0];
	//		const fvec3& X1 = RestPositionList[Ind1];
	//		const fvec3& X2 = RestPositionList[Ind2];
	//		const fvec3& X3 = RestPositionList[Ind3];

	//		const fvec3& x0 = TempPositionList[Ind0];
	//		const fvec3& x1 = TempPositionList[Ind1];
	//		const fvec3& x2 = TempPositionList[Ind2];
	//		const fvec3& x3 = TempPositionList[Ind3];

	//		const fvec3& s0 = extPositionList[Ind0];
	//		const fvec3& s1 = extPositionList[Ind1];
	//		const fvec3& s2 = extPositionList[Ind2];
	//		const fvec3& s3 = extPositionList[Ind3];

	//		const fmat3& A = AList[i];
	//		const fmat3 F  = mat3(x1 - x0, x2 - x0, x3 - x0) * A;
	//		const fmat3 E  = 0.5 * (F.transpose() * F - fmat3::indentity());

	//		const float& V = VList[i];

	//		float W;
	//		fvec3 dx0, dx1, dx2, dx3;

	//		//if (MaterialInd == 0) {
	//		//	FemElasticDxStVenant(F, E, A, V, lambda, mu, W, dx0, dx1, dx2, dx3);
	//		//} else if (MaterialInd == 1) {
	//		fmat3 AFinv, AAt;
	//		float logJ;
	//		bool ValidEnergy;
	//		ValidEnergy = FemElasticHessianNeoHookean(F, E, A, V, lambda, mu, W, logJ, dx0, dx1, dx2, dx3, AFinv, AAt);

	//		if (!ValidEnergy)
	//			continue;

	//		const float rowAAt1sum = AAt.m[0] + AAt.m[1] + AAt.m[2];
	//		const float rowAAt2sum = AAt.m[3] + AAt.m[4] + AAt.m[5];
	//		const float rowAAt3sum = AAt.m[6] + AAt.m[7] + AAt.m[8];

	//		const Eigen::Vector3f vectorizedAFinv0(
	//		    -AFinv.m[0] - AFinv.m[3] - AFinv.m[6],
	//		    -AFinv.m[1] - AFinv.m[4] - AFinv.m[7],
	//		    -AFinv.m[2] - AFinv.m[5] - AFinv.m[8]);
	//		const Eigen::Vector3f vectorizedAFinv1(AFinv.m[0], AFinv.m[1], AFinv.m[2]);
	//		const Eigen::Vector3f vectorizedAFinv2(AFinv.m[3], AFinv.m[4], AFinv.m[5]);
	//		const Eigen::Vector3f vectorizedAFinv3(AFinv.m[6], AFinv.m[7], AFinv.m[8]);

	//		Eigen::Matrix<float, 12, 12> Hessian = (mass / (dt * dt)) * Eigen::Matrix<float, 12, 12>::Identity();

	//		//Hessian += (lambda + mu - lambda * logJ) * vectorizedAFinv * vectorizedAFinv.transpose();

	//		/*
	//		Hessian.block<3, 3>(3 * 0, 3 * 0) += V * (lambda + mu - lambda * logJ) * (vectorizedAFinv0 * vectorizedAFinv0.transpose());

	//		Hessian.block<3, 3>(3 * 0, 3 * 1) += V * (lambda + mu - lambda * logJ) * (vectorizedAFinv1 * vectorizedAFinv0.transpose());
	//		Hessian.block<3, 3>(3 * 1, 3 * 0) += V * (lambda + mu - lambda * logJ) * (vectorizedAFinv0 * vectorizedAFinv1.transpose());

	//		Hessian.block<3, 3>(3 * 0, 3 * 2) += V * (lambda + mu - lambda * logJ) * (vectorizedAFinv2 * vectorizedAFinv0.transpose());
	//		Hessian.block<3, 3>(3 * 2, 3 * 0) += V * (lambda + mu - lambda * logJ) * (vectorizedAFinv0 * vectorizedAFinv2.transpose());

	//		Hessian.block<3, 3>(3 * 0, 3 * 3) += V * (lambda + mu - lambda * logJ) * (vectorizedAFinv3 * vectorizedAFinv0.transpose());
	//		Hessian.block<3, 3>(3 * 3, 3 * 0) += V * (lambda + mu - lambda * logJ) * (vectorizedAFinv0 * vectorizedAFinv3.transpose());

	//		Hessian.block<3, 3>(3 * 1, 3 * 1) += V * (lambda + mu - lambda * logJ) * (vectorizedAFinv1 * vectorizedAFinv1.transpose());

	//		Hessian.block<3, 3>(3 * 1, 3 * 2) += V * (lambda + mu - lambda * logJ) * (vectorizedAFinv2 * vectorizedAFinv1.transpose());
	//		Hessian.block<3, 3>(3 * 2, 3 * 1) += V * (lambda + mu - lambda * logJ) * (vectorizedAFinv1 * vectorizedAFinv2.transpose());

	//		Hessian.block<3, 3>(3 * 1, 3 * 3) += V * (lambda + mu - lambda * logJ) * (vectorizedAFinv3 * vectorizedAFinv1.transpose());
	//		Hessian.block<3, 3>(3 * 3, 3 * 1) += V * (lambda + mu - lambda * logJ) * (vectorizedAFinv1 * vectorizedAFinv3.transpose());

	//		Hessian.block<3, 3>(3 * 2, 3 * 2) += V * (lambda + mu - lambda * logJ) * (vectorizedAFinv2 * vectorizedAFinv2.transpose());

	//		Hessian.block<3, 3>(3 * 2, 3 * 3) += V * (lambda + mu - lambda * logJ) * (vectorizedAFinv3 * vectorizedAFinv2.transpose());
	//		Hessian.block<3, 3>(3 * 3, 3 * 2) += V * (lambda + mu - lambda * logJ) * (vectorizedAFinv2 * vectorizedAFinv3.transpose());

	//		Hessian.block<3, 3>(3 * 3, 3 * 3) += V * (lambda + mu - lambda * logJ) * (vectorizedAFinv3 * vectorizedAFinv3.transpose());
	//		*/

	//		Hessian.block<3, 3>(3 * 0, 3 * 0) += V * (lambda + mu - lambda * logJ) * (vectorizedAFinv0 * vectorizedAFinv0.transpose());

	//		Hessian.block<3, 3>(3 * 0, 3 * 1) += V * (lambda * (vectorizedAFinv0 * vectorizedAFinv1.transpose()) + (mu - lambda * logJ) * (vectorizedAFinv1 * vectorizedAFinv0.transpose()));
	//		Hessian.block<3, 3>(3 * 1, 3 * 0) += V * (lambda * (vectorizedAFinv1 * vectorizedAFinv0.transpose()) + (mu - lambda * logJ) * (vectorizedAFinv0 * vectorizedAFinv1.transpose()));

	//		Hessian.block<3, 3>(3 * 0, 3 * 2) += V * (lambda * (vectorizedAFinv0 * vectorizedAFinv2.transpose()) + (mu - lambda * logJ) * (vectorizedAFinv2 * vectorizedAFinv0.transpose()));
	//		Hessian.block<3, 3>(3 * 2, 3 * 0) += V * (lambda * (vectorizedAFinv2 * vectorizedAFinv0.transpose()) + (mu - lambda * logJ) * (vectorizedAFinv0 * vectorizedAFinv2.transpose()));

	//		Hessian.block<3, 3>(3 * 0, 3 * 3) += V * (lambda * (vectorizedAFinv0 * vectorizedAFinv3.transpose()) + (mu - lambda * logJ) * (vectorizedAFinv3 * vectorizedAFinv0.transpose()));
	//		Hessian.block<3, 3>(3 * 3, 3 * 0) += V * (lambda * (vectorizedAFinv3 * vectorizedAFinv0.transpose()) + (mu - lambda * logJ) * (vectorizedAFinv0 * vectorizedAFinv3.transpose()));

	//		Hessian.block<3, 3>(3 * 1, 3 * 1) += V * (lambda + mu - lambda * logJ) * (vectorizedAFinv1 * vectorizedAFinv1.transpose());

	//		Hessian.block<3, 3>(3 * 1, 3 * 2) += V * (lambda * (vectorizedAFinv1 * vectorizedAFinv2.transpose()) + (mu - lambda * logJ) * (vectorizedAFinv2 * vectorizedAFinv1.transpose()));
	//		Hessian.block<3, 3>(3 * 2, 3 * 1) += V * (lambda * (vectorizedAFinv2 * vectorizedAFinv1.transpose()) + (mu - lambda * logJ) * (vectorizedAFinv1 * vectorizedAFinv2.transpose()));

	//		Hessian.block<3, 3>(3 * 1, 3 * 3) += V * (lambda * (vectorizedAFinv1 * vectorizedAFinv3.transpose()) + (mu - lambda * logJ) * (vectorizedAFinv3 * vectorizedAFinv1.transpose()));
	//		Hessian.block<3, 3>(3 * 3, 3 * 1) += V * (lambda * (vectorizedAFinv3 * vectorizedAFinv1.transpose()) + (mu - lambda * logJ) * (vectorizedAFinv1 * vectorizedAFinv3.transpose()));

	//		Hessian.block<3, 3>(3 * 2, 3 * 2) += V * (lambda + mu - lambda * logJ) * (vectorizedAFinv2 * vectorizedAFinv2.transpose());

	//		Hessian.block<3, 3>(3 * 2, 3 * 3) += V * (lambda * (vectorizedAFinv2 * vectorizedAFinv3.transpose()) + (mu - lambda * logJ) * (vectorizedAFinv3 * vectorizedAFinv2.transpose()));
	//		Hessian.block<3, 3>(3 * 3, 3 * 2) += V * (lambda * (vectorizedAFinv3 * vectorizedAFinv2.transpose()) + (mu - lambda * logJ) * (vectorizedAFinv2 * vectorizedAFinv3.transpose()));

	//		Hessian.block<3, 3>(3 * 3, 3 * 3) += V * (lambda + mu - lambda * logJ) * (vectorizedAFinv3 * vectorizedAFinv3.transpose());

	//		//////

	//		Hessian.block<3, 3>(3 * 0, 3 * 0) += V * Eigen::Matrix3f::Identity() * (rowAAt1sum + rowAAt2sum + rowAAt3sum) * mu;

	//		Hessian.block<3, 3>(3 * 0, 3 * 1) -= V * Eigen::Matrix3f::Identity() * rowAAt1sum * mu;
	//		Hessian.block<3, 3>(3 * 1, 3 * 0) -= V * Eigen::Matrix3f::Identity() * rowAAt1sum * mu;

	//		Hessian.block<3, 3>(3 * 0, 3 * 2) -= V * Eigen::Matrix3f::Identity() * rowAAt2sum * mu;
	//		Hessian.block<3, 3>(3 * 2, 3 * 0) -= V * Eigen::Matrix3f::Identity() * rowAAt2sum * mu;

	//		Hessian.block<3, 3>(3 * 0, 3 * 3) -= V * Eigen::Matrix3f::Identity() * rowAAt3sum * mu;
	//		Hessian.block<3, 3>(3 * 3, 3 * 0) -= V * Eigen::Matrix3f::Identity() * rowAAt3sum * mu;

	//		Hessian.block<3, 3>(3 * 1, 3 * 1) += V * Eigen::Matrix3f::Identity() * AAt.m[0] * mu;

	//		Hessian.block<3, 3>(3 * 1, 3 * 2) += V * Eigen::Matrix3f::Identity() * AAt.m[1] * mu;
	//		Hessian.block<3, 3>(3 * 2, 3 * 1) += V * Eigen::Matrix3f::Identity() * AAt.m[3] * mu;

	//		Hessian.block<3, 3>(3 * 1, 3 * 3) += V * Eigen::Matrix3f::Identity() * AAt.m[2] * mu;
	//		Hessian.block<3, 3>(3 * 3, 3 * 1) += V * Eigen::Matrix3f::Identity() * AAt.m[6] * mu;

	//		Hessian.block<3, 3>(3 * 2, 3 * 2) += V * Eigen::Matrix3f::Identity() * AAt.m[4] * mu;

	//		Hessian.block<3, 3>(3 * 2, 3 * 3) += V * Eigen::Matrix3f::Identity() * AAt.m[5] * mu;
	//		Hessian.block<3, 3>(3 * 3, 3 * 2) += V * Eigen::Matrix3f::Identity() * AAt.m[7] * mu;

	//		Hessian.block<3, 3>(3 * 3, 3 * 3) += V * Eigen::Matrix3f::Identity() * AAt.m[8] * mu;

	//		//Hessian *= V;

	//		/*
	//		if (i == 0) {
	//			std::cout << Hessian << std::endl;

	//			fvec3 rowAFinv1 = fvec3(AFinv.m[0], AFinv.m[1], AFinv.m[2]);
	//			fvec3 rowAFinv2 = fvec3(AFinv.m[3], AFinv.m[4], AFinv.m[5]);
	//			fvec3 rowAFinv3 = fvec3(AFinv.m[6], AFinv.m[7], AFinv.m[8]);

	//			std::cout << V * (lambda + mu - lambda * logJ) * ((rowAFinv1 + rowAFinv2 + rowAFinv3).tensorproduct(rowAFinv1 + rowAFinv2 + rowAFinv3)) << std::endl;

	//			std::cout << V * (lambda + mu - lambda * logJ) * (rowAFinv1.tensorproduct(rowAFinv1)) << std::endl;
	//			std::cout << V * (lambda + mu - lambda * logJ) * (rowAFinv2.tensorproduct(rowAFinv2)) << std::endl;
	//			std::cout << V * (lambda + mu - lambda * logJ) * (rowAFinv3.tensorproduct(rowAFinv3)) << std::endl;
	//		}
	//		*/

	//		//fvec3 rowAFinv1 = fvec3(AFinv.m[0], AFinv.m[1], AFinv.m[2]);
	//		//fvec3 rowAFinv2 = fvec3(AFinv.m[3], AFinv.m[4], AFinv.m[5]);
	//		//fvec3 rowAFinv3 = fvec3(AFinv.m[6], AFinv.m[7], AFinv.m[8]);

	//		//H1 = (lambda + mu - lambda * logJ) * (rowAFinv1.tensorproduct(rowAFinv1)) + mu * AAt.m[0] * fmat3::indentity();
	//		//H2 = (lambda + mu - lambda * logJ) * (rowAFinv2.tensorproduct(rowAFinv2)) + mu * AAt.m[1] * fmat3::indentity();
	//		//H3 = (lambda + mu - lambda * logJ) * (rowAFinv3.tensorproduct(rowAFinv3)) + mu * AAt.m[2] * fmat3::indentity();
	//		//H0 = (lambda + mu - lambda * logJ) * ((rowAFinv1 + rowAFinv2 + rowAFinv3).tensorproduct(rowAFinv1 + rowAFinv2 + rowAFinv3)) + mu * AAtsum * fmat3::indentity();

	//		//} else if (MaterialInd == 2) {
	//		//	FemElasticDxCoRotational(F, E, A, qList[i], V, lambda, mu, W, dx0, dx1, dx2, dx3);
	//		//}

	//		if (ValidEnergy && W > 0.001) {

	//			Eigen::Matrix<float, 12, 1> Dx;
	//			Dx(0)  = -(mass / (dt * dt)) * (x0.x - s0.x) - dx0.x;
	//			Dx(1)  = -(mass / (dt * dt)) * (x0.y - s0.y) - dx0.y;
	//			Dx(2)  = -(mass / (dt * dt)) * (x0.z - s0.z) - dx0.z;
	//			Dx(3)  = -(mass / (dt * dt)) * (x1.x - s1.x) - dx1.x;
	//			Dx(4)  = -(mass / (dt * dt)) * (x1.y - s1.y) - dx1.y;
	//			Dx(5)  = -(mass / (dt * dt)) * (x1.z - s1.z) - dx1.z;
	//			Dx(6)  = -(mass / (dt * dt)) * (x2.x - s2.x) - dx2.x;
	//			Dx(7)  = -(mass / (dt * dt)) * (x2.y - s2.y) - dx2.y;
	//			Dx(8)  = -(mass / (dt * dt)) * (x2.z - s2.z) - dx2.z;
	//			Dx(9)  = -(mass / (dt * dt)) * (x3.x - s3.x) - dx3.x;
	//			Dx(10) = -(mass / (dt * dt)) * (x3.y - s3.y) - dx3.y;
	//			Dx(11) = -(mass / (dt * dt)) * (x3.z - s3.z) - dx3.z;

	//			Eigen::Matrix<float, 12, 1> HinvDx = Hessian.partialPivLu().solve(Dx);

	//			//std::cout << "!!!!!!!!!!!!!!!!!!!!!" << std::endl;
	//			//std::cout << Hessian << std::endl;
	//			//std::cout << Dx << std::endl;
	//			//std::cout << HinvDx << std::endl;

	//			//std::cout << Dx(0) + Dx(3) + Dx(6) + Dx(9) << std::endl;
	//			//std::cout << Dx(1) + Dx(4) + Dx(7) + Dx(10) << std::endl;
	//			//std::cout << Dx(2) + Dx(5) + Dx(8) + Dx(11) << std::endl;

	//			//std::cout << HinvDx(0) + HinvDx(3) + HinvDx(6) + HinvDx(9) << std::endl;
	//			//std::cout << HinvDx(1) + HinvDx(4) + HinvDx(7) + HinvDx(10) << std::endl;
	//			//std::cout << HinvDx(2) + HinvDx(5) + HinvDx(8) + HinvDx(11) << std::endl;

	//			float CurrentEnergy;
	//			EvaluateLocalEnergy(i, HinvDx, 0.0, CurrentEnergy);

	//			float DxoHinvDx = Dx.dot(HinvDx);

	//			for (float alpha = 1.00; alpha > 0.0001; alpha *= 0.1) {

	//				float Energy;

	//				if (EvaluateLocalEnergy(i, HinvDx, alpha, Energy) && Energy < (CurrentEnergy - 0.3 * alpha * DxoHinvDx)) {

	//					//std::cout << TempPositionList[Ind0] << std::endl;

	//					TempPositionList[Ind0] = TempPositionList[Ind0] + alpha * fvec3(HinvDx(0), HinvDx(1), HinvDx(2));
	//					TempPositionList[Ind1] = TempPositionList[Ind1] + alpha * fvec3(HinvDx(3), HinvDx(4), HinvDx(5));
	//					TempPositionList[Ind2] = TempPositionList[Ind2] + alpha * fvec3(HinvDx(6), HinvDx(7), HinvDx(8));
	//					TempPositionList[Ind3] = TempPositionList[Ind3] + alpha * fvec3(HinvDx(9), HinvDx(10), HinvDx(11));

	//					//std::cout << "------------------------" << std::endl;
	//					//std::cout << alpha * fvec3(HinvDx(0), HinvDx(1), HinvDx(2)) << " " << dx0 << std::endl;
	//					//std::cout << alpha * fvec3(HinvDx(3), HinvDx(4), HinvDx(5)) << " " << dx1 << std::endl;
	//					//std::cout << alpha * fvec3(HinvDx(6), HinvDx(7), HinvDx(8)) << " " << dx2 << std::endl;
	//					//std::cout << alpha * fvec3(HinvDx(9), HinvDx(10), HinvDx(11)) << " " << dx3 << std::endl;

	//					//std::cout << HinvDx(0) + HinvDx(3) + HinvDx(6) + HinvDx(9) << std::endl;
	//					//std::cout << HinvDx(1) + HinvDx(4) + HinvDx(7) + HinvDx(10) << std::endl;
	//					//std::cout << HinvDx(2) + HinvDx(5) + HinvDx(8) + HinvDx(11) << std::endl;

	//					////std::cout << TempPositionList[Ind0] << std::endl;
	//					//std::cout << std::endl;

	//					//std::cout << alpha << " " << CurrentEnergy - Energy << std::endl;
	//					//std::cout << alpha * HinvDx << std::endl;
	//					//std::cout << HinvDx.dot(Dx) / (HinvDx.norm() * Dx.norm()) << std::endl;
	//					//std::cout << std::endl;

	//					break;
	//				}
	//			}
	//		}
	//	}
	//}

	void BlockInsertion(uint32_t i, uint32_t j, const Eigen::Matrix<float, 3, 3>& Mat)
	{
		for (uint32_t k = 0; k < 3; k++)
			for (uint32_t l = 0; l < 3; l++)
				GlobalHessian.coeffRef(3 * i + k, 3 * j + l) += Mat(k, l);
	}

	void BlockInsertion(uint32_t i, uint32_t j, const float& value)
	{
		GlobalHessian.coeffRef(3 * i + 0, 3 * j + 0) += value;
		GlobalHessian.coeffRef(3 * i + 1, 3 * j + 1) += value;
		GlobalHessian.coeffRef(3 * i + 2, 3 * j + 2) += value;
	}

	bool EvaluateEnergy(const float& alpha, float& Energy)
	{
		Energy = 0.0;

		for (uint32_t i = 0; i < vertsize; i++)
			Energy += (mass / (2.0 * dt * dt)) * ((alpha * fvec3(DecentDirecton(3 * i + 0), DecentDirecton(3 * i + 1), DecentDirecton(3 * i + 2)) + TempPositionList[i] - extPositionList[i]).sqlength());

		for (uint32_t i = 0; i < elementsize / 4; i++) {

			if (!ValidationList[i])
				continue;

			const uint32_t Ind0 = elementlist[4 * i + 0];
			const uint32_t Ind1 = elementlist[4 * i + 1];
			const uint32_t Ind2 = elementlist[4 * i + 2];
			const uint32_t Ind3 = elementlist[4 * i + 3];

			const fvec3& X0 = RestPositionList[Ind0];
			const fvec3& X1 = RestPositionList[Ind1];
			const fvec3& X2 = RestPositionList[Ind2];
			const fvec3& X3 = RestPositionList[Ind3];

			const fvec3 x0 = TempPositionList[Ind0] + alpha * fvec3(DecentDirecton(3 * Ind0 + 0), DecentDirecton(3 * Ind0 + 1), DecentDirecton(3 * Ind0 + 2));
			const fvec3 x1 = TempPositionList[Ind1] + alpha * fvec3(DecentDirecton(3 * Ind1 + 0), DecentDirecton(3 * Ind1 + 1), DecentDirecton(3 * Ind1 + 2));
			const fvec3 x2 = TempPositionList[Ind2] + alpha * fvec3(DecentDirecton(3 * Ind2 + 0), DecentDirecton(3 * Ind2 + 1), DecentDirecton(3 * Ind2 + 2));
			const fvec3 x3 = TempPositionList[Ind3] + alpha * fvec3(DecentDirecton(3 * Ind3 + 0), DecentDirecton(3 * Ind3 + 1), DecentDirecton(3 * Ind3 + 2));

			const fmat3& A = AList[i];
			const fmat3 F  = mat3(x1 - x0, x2 - x0, x3 - x0) * A;
			const fmat3 E  = 0.5 * (F.transpose() * F - fmat3::indentity());

			const float& V = VList[i];

			float W;
			bool ValidEnergy;
			ValidEnergy = FemElasticEnergyNeoHookean(F, E, A, V, lambda, mu, W);

			if (ValidEnergy)
				Energy += W;
			else
				return false;
		}

		return true;
	}

	void EvaluateElementHG()
	{

		ClearHessian();

#pragma omp parallel for
		for (uint32_t i = 0; i < vertsize; i++) {
			GlobalHessian.coeffRef(3 * i + 0, 3 * i + 0) = (mass / (dt * dt));
			GlobalHessian.coeffRef(3 * i + 1, 3 * i + 1) = (mass / (dt * dt));
			GlobalHessian.coeffRef(3 * i + 2, 3 * i + 2) = (mass / (dt * dt));

			GlobalB(3 * i + 0) = (mass / (dt * dt)) * (TempPositionList[i].x - extPositionList[i].x);
			GlobalB(3 * i + 1) = (mass / (dt * dt)) * (TempPositionList[i].y - extPositionList[i].y);
			GlobalB(3 * i + 2) = (mass / (dt * dt)) * (TempPositionList[i].z - extPositionList[i].z);
		}

#pragma omp parallel for
		for (uint32_t i = 0; i < elementsize / 4; i++) {

			const uint32_t& Ind0 = elementlist[4 * i + 0];
			const uint32_t& Ind1 = elementlist[4 * i + 1];
			const uint32_t& Ind2 = elementlist[4 * i + 2];
			const uint32_t& Ind3 = elementlist[4 * i + 3];

			const fvec3& X0 = RestPositionList[Ind0];
			const fvec3& X1 = RestPositionList[Ind1];
			const fvec3& X2 = RestPositionList[Ind2];
			const fvec3& X3 = RestPositionList[Ind3];

			const fvec3& x0 = TempPositionList[Ind0];
			const fvec3& x1 = TempPositionList[Ind1];
			const fvec3& x2 = TempPositionList[Ind2];
			const fvec3& x3 = TempPositionList[Ind3];

			const fvec3& s0 = extPositionList[Ind0];
			const fvec3& s1 = extPositionList[Ind1];
			const fvec3& s2 = extPositionList[Ind2];
			const fvec3& s3 = extPositionList[Ind3];

			const fmat3& A = AList[i];
			const fmat3 F  = mat3(x1 - x0, x2 - x0, x3 - x0) * A;
			const fmat3 E  = 0.5 * (F.transpose() * F - fmat3::indentity());

			const float& V = VList[i];

			float W;
			fvec3 dx0, dx1, dx2, dx3;

			fmat3 AFinv, AAt;
			float logJ;
			bool ValidEnergy;
			ValidEnergy = FemElasticHessianNeoHookean(F, E, A, V, lambda, mu, W, logJ, dx0, dx1, dx2, dx3, AFinv, AAt);

			if (!ValidEnergy) {
				ValidationList[i] = false;
				continue;
			}

			ValidationList[i] = true;

			const float rowAAt1sum = AAt.m[0] + AAt.m[1] + AAt.m[2];
			const float rowAAt2sum = AAt.m[3] + AAt.m[4] + AAt.m[5];
			const float rowAAt3sum = AAt.m[6] + AAt.m[7] + AAt.m[8];

			const Eigen::Vector3f vectorizedAFinv0(
			    -AFinv.m[0] - AFinv.m[3] - AFinv.m[6],
			    -AFinv.m[1] - AFinv.m[4] - AFinv.m[7],
			    -AFinv.m[2] - AFinv.m[5] - AFinv.m[8]);
			const Eigen::Vector3f vectorizedAFinv1(AFinv.m[0], AFinv.m[1], AFinv.m[2]);
			const Eigen::Vector3f vectorizedAFinv2(AFinv.m[3], AFinv.m[4], AFinv.m[5]);
			const Eigen::Vector3f vectorizedAFinv3(AFinv.m[6], AFinv.m[7], AFinv.m[8]);

			ElementHGList[i].Hessians[4 * 0 + 0] = (mass / (dt * dt)) * Eigen::Matrix<float, 3, 3>::Identity();
			ElementHGList[i].Hessians[4 * 1 + 1] = (mass / (dt * dt)) * Eigen::Matrix<float, 3, 3>::Identity();
			ElementHGList[i].Hessians[4 * 2 + 2] = (mass / (dt * dt)) * Eigen::Matrix<float, 3, 3>::Identity();
			ElementHGList[i].Hessians[4 * 3 + 3] = (mass / (dt * dt)) * Eigen::Matrix<float, 3, 3>::Identity();

			ElementHGList[i].Hessians[4 * 0 + 0] += V * (lambda + mu - lambda * logJ) * (vectorizedAFinv0 * vectorizedAFinv0.transpose());

			ElementHGList[i].Hessians[4 * 0 + 1] = V * (lambda * (vectorizedAFinv0 * vectorizedAFinv1.transpose()) + (mu - lambda * logJ) * (vectorizedAFinv1 * vectorizedAFinv0.transpose()));
			ElementHGList[i].Hessians[4 * 1 + 0] = V * (lambda * (vectorizedAFinv1 * vectorizedAFinv0.transpose()) + (mu - lambda * logJ) * (vectorizedAFinv0 * vectorizedAFinv1.transpose()));

			ElementHGList[i].Hessians[4 * 0 + 2] = V * (lambda * (vectorizedAFinv0 * vectorizedAFinv2.transpose()) + (mu - lambda * logJ) * (vectorizedAFinv2 * vectorizedAFinv0.transpose()));
			ElementHGList[i].Hessians[4 * 2 + 0] = V * (lambda * (vectorizedAFinv2 * vectorizedAFinv0.transpose()) + (mu - lambda * logJ) * (vectorizedAFinv0 * vectorizedAFinv2.transpose()));

			ElementHGList[i].Hessians[4 * 0 + 3] = V * (lambda * (vectorizedAFinv0 * vectorizedAFinv3.transpose()) + (mu - lambda * logJ) * (vectorizedAFinv3 * vectorizedAFinv0.transpose()));
			ElementHGList[i].Hessians[4 * 3 + 0] = V * (lambda * (vectorizedAFinv3 * vectorizedAFinv0.transpose()) + (mu - lambda * logJ) * (vectorizedAFinv0 * vectorizedAFinv3.transpose()));

			ElementHGList[i].Hessians[4 * 1 + 1] += V * (lambda + mu - lambda * logJ) * (vectorizedAFinv1 * vectorizedAFinv1.transpose());

			ElementHGList[i].Hessians[4 * 1 + 2] = V * (lambda * (vectorizedAFinv1 * vectorizedAFinv2.transpose()) + (mu - lambda * logJ) * (vectorizedAFinv2 * vectorizedAFinv1.transpose()));
			ElementHGList[i].Hessians[4 * 2 + 1] = V * (lambda * (vectorizedAFinv2 * vectorizedAFinv1.transpose()) + (mu - lambda * logJ) * (vectorizedAFinv1 * vectorizedAFinv2.transpose()));

			ElementHGList[i].Hessians[4 * 1 + 3] = V * (lambda * (vectorizedAFinv1 * vectorizedAFinv3.transpose()) + (mu - lambda * logJ) * (vectorizedAFinv3 * vectorizedAFinv1.transpose()));
			ElementHGList[i].Hessians[4 * 3 + 1] = V * (lambda * (vectorizedAFinv3 * vectorizedAFinv1.transpose()) + (mu - lambda * logJ) * (vectorizedAFinv1 * vectorizedAFinv3.transpose()));

			ElementHGList[i].Hessians[4 * 2 + 2] += V * (lambda + mu - lambda * logJ) * (vectorizedAFinv2 * vectorizedAFinv2.transpose());

			ElementHGList[i].Hessians[4 * 2 + 3] = V * (lambda * (vectorizedAFinv2 * vectorizedAFinv3.transpose()) + (mu - lambda * logJ) * (vectorizedAFinv3 * vectorizedAFinv2.transpose()));
			ElementHGList[i].Hessians[4 * 3 + 2] = V * (lambda * (vectorizedAFinv3 * vectorizedAFinv2.transpose()) + (mu - lambda * logJ) * (vectorizedAFinv2 * vectorizedAFinv3.transpose()));

			ElementHGList[i].Hessians[4 * 3 + 3] += V * (lambda + mu - lambda * logJ) * (vectorizedAFinv3 * vectorizedAFinv3.transpose());

			//////

			ElementHGList[i].Hessians[4 * 0 + 0] += V * Eigen::Matrix3f::Identity() * (rowAAt1sum + rowAAt2sum + rowAAt3sum) * mu;

			ElementHGList[i].Hessians[4 * 0 + 1] -= V * Eigen::Matrix3f::Identity() * rowAAt1sum * mu;
			ElementHGList[i].Hessians[4 * 1 + 0] -= V * Eigen::Matrix3f::Identity() * rowAAt1sum * mu;

			ElementHGList[i].Hessians[4 * 0 + 2] -= V * Eigen::Matrix3f::Identity() * rowAAt2sum * mu;
			ElementHGList[i].Hessians[4 * 2 + 0] -= V * Eigen::Matrix3f::Identity() * rowAAt2sum * mu;

			ElementHGList[i].Hessians[4 * 0 + 3] -= V * Eigen::Matrix3f::Identity() * rowAAt3sum * mu;
			ElementHGList[i].Hessians[4 * 3 + 0] -= V * Eigen::Matrix3f::Identity() * rowAAt3sum * mu;

			ElementHGList[i].Hessians[4 * 1 + 1] += V * Eigen::Matrix3f::Identity() * AAt.m[0] * mu;

			ElementHGList[i].Hessians[4 * 1 + 2] += V * Eigen::Matrix3f::Identity() * AAt.m[1] * mu;
			ElementHGList[i].Hessians[4 * 2 + 1] += V * Eigen::Matrix3f::Identity() * AAt.m[3] * mu;

			ElementHGList[i].Hessians[4 * 1 + 3] += V * Eigen::Matrix3f::Identity() * AAt.m[2] * mu;
			ElementHGList[i].Hessians[4 * 3 + 1] += V * Eigen::Matrix3f::Identity() * AAt.m[6] * mu;

			ElementHGList[i].Hessians[4 * 2 + 2] += V * Eigen::Matrix3f::Identity() * AAt.m[4] * mu;

			ElementHGList[i].Hessians[4 * 2 + 3] += V * Eigen::Matrix3f::Identity() * AAt.m[5] * mu;
			ElementHGList[i].Hessians[4 * 3 + 2] += V * Eigen::Matrix3f::Identity() * AAt.m[7] * mu;

			ElementHGList[i].Hessians[4 * 3 + 3] += V * Eigen::Matrix3f::Identity() * AAt.m[8] * mu;

			ElementHGList[i].Gradient[0](0) = -dx0.x;
			ElementHGList[i].Gradient[0](1) = -dx0.y;
			ElementHGList[i].Gradient[0](2) = -dx0.z;

			ElementHGList[i].Gradient[1](0) = -dx1.x;
			ElementHGList[i].Gradient[1](1) = -dx1.y;
			ElementHGList[i].Gradient[1](2) = -dx1.z;

			ElementHGList[i].Gradient[2](0) = -dx2.x;
			ElementHGList[i].Gradient[2](1) = -dx2.y;
			ElementHGList[i].Gradient[2](2) = -dx2.z;

			ElementHGList[i].Gradient[3](0) = -dx3.x;
			ElementHGList[i].Gradient[3](1) = -dx3.y;
			ElementHGList[i].Gradient[3](2) = -dx3.z;
		}
	}

	void MergeHessian()
	{
		for (uint32_t i = 0; i < elementsize / 4; i++) {

			const uint32_t& Ind0 = elementlist[4 * i + 0];
			const uint32_t& Ind1 = elementlist[4 * i + 1];
			const uint32_t& Ind2 = elementlist[4 * i + 2];
			const uint32_t& Ind3 = elementlist[4 * i + 3];

			BlockInsertion(Ind0, Ind0, ElementHGList[i].Hessians[4 * 0 + 0]);

			BlockInsertion(Ind1, Ind0, ElementHGList[i].Hessians[4 * 1 + 0]);
			BlockInsertion(Ind0, Ind1, ElementHGList[i].Hessians[4 * 0 + 1]);

			BlockInsertion(Ind2, Ind0, ElementHGList[i].Hessians[4 * 2 + 0]);
			BlockInsertion(Ind0, Ind2, ElementHGList[i].Hessians[4 * 0 + 2]);

			BlockInsertion(Ind3, Ind0, ElementHGList[i].Hessians[4 * 3 + 0]);
			BlockInsertion(Ind0, Ind3, ElementHGList[i].Hessians[4 * 0 + 3]);

			BlockInsertion(Ind1, Ind1, ElementHGList[i].Hessians[4 * 1 + 1]);

			BlockInsertion(Ind2, Ind1, ElementHGList[i].Hessians[4 * 2 + 1]);
			BlockInsertion(Ind1, Ind2, ElementHGList[i].Hessians[4 * 1 + 2]);

			BlockInsertion(Ind3, Ind1, ElementHGList[i].Hessians[4 * 3 + 1]);
			BlockInsertion(Ind1, Ind3, ElementHGList[i].Hessians[4 * 1 + 3]);

			BlockInsertion(Ind2, Ind2, ElementHGList[i].Hessians[4 * 2 + 2]);

			BlockInsertion(Ind2, Ind3, ElementHGList[i].Hessians[4 * 2 + 3]);
			BlockInsertion(Ind3, Ind2, ElementHGList[i].Hessians[4 * 3 + 2]);

			BlockInsertion(Ind3, Ind3, ElementHGList[i].Hessians[4 * 3 + 3]);

			//////
		}

		Solver.compute(GlobalHessian);
	}

	void MergeGradient()
	{
		for (uint32_t i = 0; i < elementsize / 4; i++) {

			const uint32_t& Ind0 = elementlist[4 * i + 0];
			const uint32_t& Ind1 = elementlist[4 * i + 1];
			const uint32_t& Ind2 = elementlist[4 * i + 2];
			const uint32_t& Ind3 = elementlist[4 * i + 3];

			//////

			GlobalB(3 * Ind0 + 0) += ElementHGList[i].Gradient[0](0);
			GlobalB(3 * Ind0 + 1) += ElementHGList[i].Gradient[0](1);
			GlobalB(3 * Ind0 + 2) += ElementHGList[i].Gradient[0](2);

			GlobalB(3 * Ind1 + 0) += ElementHGList[i].Gradient[1](0);
			GlobalB(3 * Ind1 + 1) += ElementHGList[i].Gradient[1](1);
			GlobalB(3 * Ind1 + 2) += ElementHGList[i].Gradient[1](2);

			GlobalB(3 * Ind2 + 0) += ElementHGList[i].Gradient[2](0);
			GlobalB(3 * Ind2 + 1) += ElementHGList[i].Gradient[2](1);
			GlobalB(3 * Ind2 + 2) += ElementHGList[i].Gradient[2](2);

			GlobalB(3 * Ind3 + 0) += ElementHGList[i].Gradient[3](0);
			GlobalB(3 * Ind3 + 1) += ElementHGList[i].Gradient[3](1);
			GlobalB(3 * Ind3 + 2) += ElementHGList[i].Gradient[3](2);
		}
	}

	void EvaluateElementHessian()
	{

		ClearHessian();
#pragma omp parallel for
		for (uint32_t i = 0; i < vertsize; i++) {
			GlobalHessian.coeffRef(3 * i + 0, 3 * i + 0) = (mass / (dt * dt));
			GlobalHessian.coeffRef(3 * i + 1, 3 * i + 1) = (mass / (dt * dt));
			GlobalHessian.coeffRef(3 * i + 2, 3 * i + 2) = (mass / (dt * dt));
		}

#pragma omp parallel for
		for (uint32_t i = 0; i < elementsize / 4; i++) {

			const uint32_t& Ind0 = elementlist[4 * i + 0];
			const uint32_t& Ind1 = elementlist[4 * i + 1];
			const uint32_t& Ind2 = elementlist[4 * i + 2];
			const uint32_t& Ind3 = elementlist[4 * i + 3];

			const fvec3& X0 = RestPositionList[Ind0];
			const fvec3& X1 = RestPositionList[Ind1];
			const fvec3& X2 = RestPositionList[Ind2];
			const fvec3& X3 = RestPositionList[Ind3];

			const fvec3& x0 = TempPositionList[Ind0];
			const fvec3& x1 = TempPositionList[Ind1];
			const fvec3& x2 = TempPositionList[Ind2];
			const fvec3& x3 = TempPositionList[Ind3];

			const fvec3& s0 = extPositionList[Ind0];
			const fvec3& s1 = extPositionList[Ind1];
			const fvec3& s2 = extPositionList[Ind2];
			const fvec3& s3 = extPositionList[Ind3];

			const fmat3& A = AList[i];
			const fmat3 F  = mat3(x1 - x0, x2 - x0, x3 - x0) * A;
			const fmat3 E  = 0.5 * (F.transpose() * F - fmat3::indentity());

			const float& V = VList[i];

			float W;
			fvec3 dx0, dx1, dx2, dx3;

			fmat3 AFinv, AAt;
			float logJ;
			bool ValidEnergy;
			ValidEnergy = FemElasticHessianNeoHookean(F, E, A, V, lambda, mu, W, logJ, dx0, dx1, dx2, dx3, AFinv, AAt);

			if (!ValidEnergy) {
				ValidationList[i] = false;
				continue;
			}

			ValidationList[i] = true;

			const float rowAAt1sum = AAt.m[0] + AAt.m[1] + AAt.m[2];
			const float rowAAt2sum = AAt.m[3] + AAt.m[4] + AAt.m[5];
			const float rowAAt3sum = AAt.m[6] + AAt.m[7] + AAt.m[8];

			const Eigen::Vector3f vectorizedAFinv0(
			    -AFinv.m[0] - AFinv.m[3] - AFinv.m[6],
			    -AFinv.m[1] - AFinv.m[4] - AFinv.m[7],
			    -AFinv.m[2] - AFinv.m[5] - AFinv.m[8]);
			const Eigen::Vector3f vectorizedAFinv1(AFinv.m[0], AFinv.m[1], AFinv.m[2]);
			const Eigen::Vector3f vectorizedAFinv2(AFinv.m[3], AFinv.m[4], AFinv.m[5]);
			const Eigen::Vector3f vectorizedAFinv3(AFinv.m[6], AFinv.m[7], AFinv.m[8]);

			ElementHGList[i].Hessians[4 * 0 + 0] = (mass / (dt * dt)) * Eigen::Matrix<float, 3, 3>::Identity();
			ElementHGList[i].Hessians[4 * 1 + 1] = (mass / (dt * dt)) * Eigen::Matrix<float, 3, 3>::Identity();
			ElementHGList[i].Hessians[4 * 2 + 2] = (mass / (dt * dt)) * Eigen::Matrix<float, 3, 3>::Identity();
			ElementHGList[i].Hessians[4 * 3 + 3] = (mass / (dt * dt)) * Eigen::Matrix<float, 3, 3>::Identity();

			ElementHGList[i].Hessians[4 * 0 + 0] += V * (lambda + mu - lambda * logJ) * (vectorizedAFinv0 * vectorizedAFinv0.transpose());

			ElementHGList[i].Hessians[4 * 0 + 1] = V * (lambda * (vectorizedAFinv0 * vectorizedAFinv1.transpose()) + (mu - lambda * logJ) * (vectorizedAFinv1 * vectorizedAFinv0.transpose()));
			ElementHGList[i].Hessians[4 * 1 + 0] = V * (lambda * (vectorizedAFinv1 * vectorizedAFinv0.transpose()) + (mu - lambda * logJ) * (vectorizedAFinv0 * vectorizedAFinv1.transpose()));

			ElementHGList[i].Hessians[4 * 0 + 2] = V * (lambda * (vectorizedAFinv0 * vectorizedAFinv2.transpose()) + (mu - lambda * logJ) * (vectorizedAFinv2 * vectorizedAFinv0.transpose()));
			ElementHGList[i].Hessians[4 * 2 + 0] = V * (lambda * (vectorizedAFinv2 * vectorizedAFinv0.transpose()) + (mu - lambda * logJ) * (vectorizedAFinv0 * vectorizedAFinv2.transpose()));

			ElementHGList[i].Hessians[4 * 0 + 3] = V * (lambda * (vectorizedAFinv0 * vectorizedAFinv3.transpose()) + (mu - lambda * logJ) * (vectorizedAFinv3 * vectorizedAFinv0.transpose()));
			ElementHGList[i].Hessians[4 * 3 + 0] = V * (lambda * (vectorizedAFinv3 * vectorizedAFinv0.transpose()) + (mu - lambda * logJ) * (vectorizedAFinv0 * vectorizedAFinv3.transpose()));

			ElementHGList[i].Hessians[4 * 1 + 1] += V * (lambda + mu - lambda * logJ) * (vectorizedAFinv1 * vectorizedAFinv1.transpose());

			ElementHGList[i].Hessians[4 * 1 + 2] = V * (lambda * (vectorizedAFinv1 * vectorizedAFinv2.transpose()) + (mu - lambda * logJ) * (vectorizedAFinv2 * vectorizedAFinv1.transpose()));
			ElementHGList[i].Hessians[4 * 2 + 1] = V * (lambda * (vectorizedAFinv2 * vectorizedAFinv1.transpose()) + (mu - lambda * logJ) * (vectorizedAFinv1 * vectorizedAFinv2.transpose()));

			ElementHGList[i].Hessians[4 * 1 + 3] = V * (lambda * (vectorizedAFinv1 * vectorizedAFinv3.transpose()) + (mu - lambda * logJ) * (vectorizedAFinv3 * vectorizedAFinv1.transpose()));
			ElementHGList[i].Hessians[4 * 3 + 1] = V * (lambda * (vectorizedAFinv3 * vectorizedAFinv1.transpose()) + (mu - lambda * logJ) * (vectorizedAFinv1 * vectorizedAFinv3.transpose()));

			ElementHGList[i].Hessians[4 * 2 + 2] += V * (lambda + mu - lambda * logJ) * (vectorizedAFinv2 * vectorizedAFinv2.transpose());

			ElementHGList[i].Hessians[4 * 2 + 3] = V * (lambda * (vectorizedAFinv2 * vectorizedAFinv3.transpose()) + (mu - lambda * logJ) * (vectorizedAFinv3 * vectorizedAFinv2.transpose()));
			ElementHGList[i].Hessians[4 * 3 + 2] = V * (lambda * (vectorizedAFinv3 * vectorizedAFinv2.transpose()) + (mu - lambda * logJ) * (vectorizedAFinv2 * vectorizedAFinv3.transpose()));

			ElementHGList[i].Hessians[4 * 3 + 3] += V * (lambda + mu - lambda * logJ) * (vectorizedAFinv3 * vectorizedAFinv3.transpose());

			//////

			ElementHGList[i].Hessians[4 * 0 + 0] += V * Eigen::Matrix3f::Identity() * (rowAAt1sum + rowAAt2sum + rowAAt3sum) * mu;

			ElementHGList[i].Hessians[4 * 0 + 1] -= V * Eigen::Matrix3f::Identity() * rowAAt1sum * mu;
			ElementHGList[i].Hessians[4 * 1 + 0] -= V * Eigen::Matrix3f::Identity() * rowAAt1sum * mu;

			ElementHGList[i].Hessians[4 * 0 + 2] -= V * Eigen::Matrix3f::Identity() * rowAAt2sum * mu;
			ElementHGList[i].Hessians[4 * 2 + 0] -= V * Eigen::Matrix3f::Identity() * rowAAt2sum * mu;

			ElementHGList[i].Hessians[4 * 0 + 3] -= V * Eigen::Matrix3f::Identity() * rowAAt3sum * mu;
			ElementHGList[i].Hessians[4 * 3 + 0] -= V * Eigen::Matrix3f::Identity() * rowAAt3sum * mu;

			ElementHGList[i].Hessians[4 * 1 + 1] += V * Eigen::Matrix3f::Identity() * AAt.m[0] * mu;

			ElementHGList[i].Hessians[4 * 1 + 2] += V * Eigen::Matrix3f::Identity() * AAt.m[1] * mu;
			ElementHGList[i].Hessians[4 * 2 + 1] += V * Eigen::Matrix3f::Identity() * AAt.m[3] * mu;

			ElementHGList[i].Hessians[4 * 1 + 3] += V * Eigen::Matrix3f::Identity() * AAt.m[2] * mu;
			ElementHGList[i].Hessians[4 * 3 + 1] += V * Eigen::Matrix3f::Identity() * AAt.m[6] * mu;

			ElementHGList[i].Hessians[4 * 2 + 2] += V * Eigen::Matrix3f::Identity() * AAt.m[4] * mu;

			ElementHGList[i].Hessians[4 * 2 + 3] += V * Eigen::Matrix3f::Identity() * AAt.m[5] * mu;
			ElementHGList[i].Hessians[4 * 3 + 2] += V * Eigen::Matrix3f::Identity() * AAt.m[7] * mu;

			ElementHGList[i].Hessians[4 * 3 + 3] += V * Eigen::Matrix3f::Identity() * AAt.m[8] * mu;
		}
	}

	void EvaluateElementGradient()
	{

#pragma omp parallel for
		for (uint32_t i = 0; i < vertsize; i++) {
			GlobalB(3 * i + 0) = (mass / (dt * dt)) * (TempPositionList[i].x - extPositionList[i].x);
			GlobalB(3 * i + 1) = (mass / (dt * dt)) * (TempPositionList[i].y - extPositionList[i].y);
			GlobalB(3 * i + 2) = (mass / (dt * dt)) * (TempPositionList[i].z - extPositionList[i].z);
		}

		for (uint32_t i = 0; i < elementsize / 4; i++) {

			const uint32_t& Ind0 = elementlist[4 * i + 0];
			const uint32_t& Ind1 = elementlist[4 * i + 1];
			const uint32_t& Ind2 = elementlist[4 * i + 2];
			const uint32_t& Ind3 = elementlist[4 * i + 3];

			const fvec3& X0 = RestPositionList[Ind0];
			const fvec3& X1 = RestPositionList[Ind1];
			const fvec3& X2 = RestPositionList[Ind2];
			const fvec3& X3 = RestPositionList[Ind3];

			const fvec3& x0 = TempPositionList[Ind0];
			const fvec3& x1 = TempPositionList[Ind1];
			const fvec3& x2 = TempPositionList[Ind2];
			const fvec3& x3 = TempPositionList[Ind3];

			const fvec3& s0 = extPositionList[Ind0];
			const fvec3& s1 = extPositionList[Ind1];
			const fvec3& s2 = extPositionList[Ind2];
			const fvec3& s3 = extPositionList[Ind3];

			const fmat3& A = AList[i];
			const fmat3 F  = mat3(x1 - x0, x2 - x0, x3 - x0) * A;
			const fmat3 E  = 0.5 * (F.transpose() * F - fmat3::indentity());

			const float& V = VList[i];

			float W;
			fvec3 dx0, dx1, dx2, dx3;

			fmat3 AFinv, AAt;
			float logJ;
			bool ValidEnergy;
			ValidEnergy = FemElasticHessianNeoHookean(F, E, A, V, lambda, mu, W, logJ, dx0, dx1, dx2, dx3, AFinv, AAt);

			if (!ValidEnergy) {
				ValidationList[i] = false;
				continue;
			}

			ValidationList[i] = true;

			//////

			ElementHGList[i].Gradient[0](0) = -dx0.x;
			ElementHGList[i].Gradient[0](1) = -dx0.y;
			ElementHGList[i].Gradient[0](2) = -dx0.z;

			ElementHGList[i].Gradient[1](0) = -dx1.x;
			ElementHGList[i].Gradient[1](1) = -dx1.y;
			ElementHGList[i].Gradient[1](2) = -dx1.z;

			ElementHGList[i].Gradient[2](0) = -dx2.x;
			ElementHGList[i].Gradient[2](1) = -dx2.y;
			ElementHGList[i].Gradient[2](2) = -dx2.z;

			ElementHGList[i].Gradient[3](0) = -dx3.x;
			ElementHGList[i].Gradient[3](1) = -dx3.y;
			ElementHGList[i].Gradient[3](2) = -dx3.z;
		}
	}

	void FemProjectGlobalSolve()
	{

		DecentDirecton = Solver.solve(GlobalB);

		float CurrentEnergy;
		EvaluateEnergy(0.0, CurrentEnergy);

		float DDoDx = 0.0;
		for (uint32_t i = 0; i < 3 * vertsize; i++)
			DDoDx += DecentDirecton(i) * GlobalB(i);

		//linesearch
		for (float alpha = 1.0; alpha > 0.0001; alpha *= 0.5) {
			float Energy;

			if (EvaluateEnergy(alpha, Energy) && Energy < (CurrentEnergy - 0.5 * alpha * DDoDx)) {
				for (uint32_t i = 0; i < vertsize; i++) {
					TempPositionList[i].x += alpha * DecentDirecton(3 * i + 0);
					TempPositionList[i].y += alpha * DecentDirecton(3 * i + 1);
					TempPositionList[i].z += alpha * DecentDirecton(3 * i + 2);
				}
				return;
			}
		}
	}

	void FemProjectXPBDGS()
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
		}
	}

	void FemProjectXPBDJC()
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

			if (MaterialInd == 0) {
				ValidEnergy = FemElasticDxStVenant(F, E, A, V, lambda, mu, W, dx0, dx1, dx2, dx3);
			} else if (MaterialInd == 1) {
				ValidEnergy = FemElasticDxNeoHookean(F, E, A, V, lambda, mu, W, dx0, dx1, dx2, dx3);

			} else if (MaterialInd == 2) {
				ValidEnergy = FemElasticDxCoRotational(F, E, A, qList[i], V, lambda, mu, W, dx0, dx1, dx2, dx3);
			}

			if (ValidEnergy && W > 0.0001) {
				float C = std::sqrt(2.0 * W);

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

	CM.FemFixInversionJC();

#pragma omp parallel for
	for (uint32_t i = 0; i < NodeSize; i++)
		CM.extPositionList[i] = CM.TempPositionList[i];

	if (solver == 0) {
		for (uint32_t x = 0; x < 2; x++) {
			for (uint32_t y = 0; y < 8; y++)
				CM.FemFixInversionGS();
			CM.EvaluateElementHG();
			CM.MergeHessian();
			CM.MergeGradient();
			CM.FemProjectGlobalSolve();
			CM.FixedProjection();
		}
	} else if (solver == 1) {
		CM.EvaluateElementHessian();
		CM.MergeHessian();
		for (uint32_t x = 0; x < 2; x++) {
			for (uint32_t y = 0; y < 8; y++)
				CM.FemFixInversionGS();
			CM.EvaluateElementGradient();
			CM.MergeGradient();
			CM.FemProjectGlobalSolve();
			CM.FixedProjection();
		}
	} else if (solver == 2) {
		for (uint32_t x = 0; x < 15; x++) {
			CM.FemProjectXPBDJC();
			CM.FemFixInversionJC();
			CM.FixedProjection();
		}
	} else if (solver == 3) {
		for (uint32_t x = 0; x < 6; x++) {
			CM.FemProjectXPBDGS();
			CM.FemFixInversionGS();
			CM.FixedProjection();
		}
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
	Eigen::setNbThreads(20);
	//Eigen::initParallel();

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

			//ImGui::Combo("Solver", &(Physics::solver), "Local Newton(Gauss-Seidel)\0Global Newton\0XPBD Jacobi\0XPBD Gauss-Seidel\0\0");
			ImGui::Combo("Solver", &(Physics::solver), "Full Newton Solver\0Semi-Implicit\0XPBD Jacobi\0XPBD Gauss-Seidel\0\0");

			//ImGui::Combo("Material", &(CM0.MaterialInd), "Saint Venant-Kirchhoff\0Neo-Hookean\0Co-Rotational\0\0");
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
