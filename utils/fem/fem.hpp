#pragma once

#include <cstdint>
#include <cmath>

#include "utils/mathfunc/mathfunc.hpp"

#include "utils/mathfunc/polardecompose.hpp"

void FemElasticDxNeoHookean(
    const fmat3& F,
    const fmat3& E,
    const fmat3& A,
    const float& Volume,
    const float& lambda,
    const float& mu, float& W,
    fvec3& dx0,
    fvec3& dx1,
    fvec3& dx2,
    fvec3& dx3)
{

	float J	   = F.det();
	float logJ = std::log(J);
	if (J < 0.00001)
		logJ = 0.0;
	W = Volume * (0.5 * mu * (F.sqlength() - 3) - mu * logJ + 0.5 * lambda * logJ * logJ);

	if (std::abs(F.det()) < 0.00001)
		W = 0.0;
	fmat3 P = Volume * (mu * F - mu * (F.inverse()).transpose() + lambda * logJ * (F.inverse()).transpose());

	fmat3 PAt = P * A.transpose();

	dx1 = fvec3(PAt.m[0], PAt.m[3], PAt.m[6]);
	dx2 = fvec3(PAt.m[1], PAt.m[4], PAt.m[7]);
	dx3 = fvec3(PAt.m[2], PAt.m[5], PAt.m[8]);
	dx0 = -(dx1 + dx2 + dx3);
}

void FemElasticDxStVenant(
    const fmat3& F,
    const fmat3& E,
    const fmat3& A,
    const float& Volume,
    const float& lambda,
    const float& mu,
    float& W,
    fvec3& dx0,
    fvec3& dx1,
    fvec3& dx2,
    fvec3& dx3)
{

	W	= Volume * (mu * E.sqlength() + 0.5 * lambda * E.trace() * E.trace());
	fmat3 P = Volume * (2 * mu * F * E + lambda * E.trace() * F);

	fmat3 PAt = P * A.transpose();

	dx1 = fvec3(PAt.m[0], PAt.m[3], PAt.m[6]);
	dx2 = fvec3(PAt.m[1], PAt.m[4], PAt.m[7]);
	dx3 = fvec3(PAt.m[2], PAt.m[5], PAt.m[8]);
	dx0 = -(dx1 + dx2 + dx3);
}

void FemElasticDxCoRotational(
    const fmat3& F,
    const fmat3& E,
    const fmat3& A,
    fquaternion& q,
    const float& Volume,
    const float& lambda,
    const float& mu,
    float& W,
    fvec3& dx0,
    fvec3& dx1,
    fvec3& dx2,
    fvec3& dx3)
{

	q = ExtractRotation(F, 3, q);

	fmat3 R	  = q.qtomat();
	fmat3 S	  = R.transpose() * F;
	float trS = S.trace();
	W	  = Volume * (0.5 * mu * (F.sqlength() - 2.0 * trS + 3) + 0.5 * lambda * (trS * trS - 6.0 * trS + 9.0));

	fmat3 P = Volume * (mu * (F - R) + lambda * (trS - 3) * R);

	fmat3 PAt = P * A.transpose();

	dx1 = fvec3(PAt.m[0], PAt.m[3], PAt.m[6]);
	dx2 = fvec3(PAt.m[1], PAt.m[4], PAt.m[7]);
	dx3 = fvec3(PAt.m[2], PAt.m[5], PAt.m[8]);
	dx0 = -(dx1 + dx2 + dx3);
}

///////

void FemElasticPNeoHookean(
    const fmat3& F,
    const fmat3& E,
    const fmat3& A,
    const float& Volume,
    const float& lambda,
    const float& mu,
    float& W,
    fmat3& P)
{

	float J	   = F.det();
	float logJ = std::log(J);
	if (J < 0.00001)
		logJ = 0.0;
	W = Volume * (0.5 * mu * (F.sqlength() - 3) - mu * logJ + 0.5 * lambda * logJ * logJ);

	if (std::abs(F.det()) < 0.00001)
		W = 0.0;
	P = Volume * (mu * F - mu * (F.inverse()).transpose() + lambda * logJ * (F.inverse()).transpose());
}

void FemElasticPStVenant(
    const fmat3& F,
    const fmat3& E,
    const fmat3& A,
    const float& Volume,
    const float& lambda,
    const float& mu,
    float& W,
    fmat3& P)
{

	W = Volume * (mu * E.sqlength() + 0.5 * lambda * E.trace() * E.trace());
	P = Volume * (2 * mu * F * E + lambda * E.trace() * F);
}

void FemElasticPCoRotational(
    const fmat3& F,
    const fmat3& E,
    const fmat3& A,
    fquaternion& q,
    const float& Volume,
    const float& lambda,
    const float& mu,
    float& W,
    fmat3& P)
{
	q = ExtractRotation(F, 3, q);

	fmat3 R	  = q.qtomat();
	fmat3 S	  = R.transpose() * F;
	float trS = S.trace();
	W	  = Volume * (0.5 * mu * (F.sqlength() - 2.0 * trS + 3) + 0.5 * lambda * (trS * trS - 6.0 * trS + 9.0));
	P	  = Volume * (mu * (F - R) + lambda * (trS - 3) * R);
}

//////

/*
void FemElasticHessianNeoHookean(
    const fmat3& F,
    const fmat3& E,
    const fmat3& A,
    const float& Volume,
    const float& lambda,
    const float& mu,
    float& W,
    fmat3& H0,
    fmat3& H1,
    fmat3& H2,
    fmat3& H3)
{

	float J	   = F.det();
	float logJ = std::log(J);
	if (J < 0.00001)
		logJ = 0.0;
	W = Volume * (0.5 * mu * (F.sqlength() - 3) - mu * logJ + 0.5 * lambda * logJ * logJ);

	if (std::abs(F.det()) < 0.00001) {
		W = 0.0;
		return;
	}

	fmat3 P = Volume * (mu * F - mu * (F.inverse()).transpose() + lambda * logJ * (F.inverse()).transpose());

	fmat3 PAt = P * A.transpose();

	fvec3 dx1 = fvec3(PAt.m[0], PAt.m[3], PAt.m[6]);
	fvec3 dx2 = fvec3(PAt.m[1], PAt.m[4], PAt.m[7]);
	fvec3 dx3 = fvec3(PAt.m[2], PAt.m[5], PAt.m[8]);
	fvec3 dx0 = -(dx1 + dx2 + dx3);

	fmat3 AFinv  = A * F.inverse();
	fmat3 AAt    = A * A.transpose();
	float AAtsum = AAt.m[0] + AAt.m[1] + AAt.m[2] + AAt.m[3] + AAt.m[4] + AAt.m[5] + AAt.m[6] + AAt.m[7] + AAt.m[8];

	fvec3 rowAFinv1 = fvec3(AFinv.m[0], AFinv.m[1], AFinv.m[2]);
	fvec3 rowAFinv2 = fvec3(AFinv.m[3], AFinv.m[4], AFinv.m[5]);
	fvec3 rowAFinv3 = fvec3(AFinv.m[6], AFinv.m[7], AFinv.m[8]);

	H1 = Volume * ((lambda + mu - lambda * logJ) * (rowAFinv1.tensorproduct(rowAFinv1)) + mu * AAt.m[0] * fmat3::indentity());
	H2 = Volume * ((lambda + mu - lambda * logJ) * (rowAFinv2.tensorproduct(rowAFinv2)) + mu * AAt.m[1] * fmat3::indentity());
	H3 = Volume * ((lambda + mu - lambda * logJ) * (rowAFinv3.tensorproduct(rowAFinv3)) + mu * AAt.m[2] * fmat3::indentity());
	H0 = Volume * ((lambda + mu - lambda * logJ) * ((rowAFinv1 + rowAFinv2 + rowAFinv3).tensorproduct(rowAFinv1 + rowAFinv2 + rowAFinv3)) + mu * AAtsum * fmat3::indentity());
}

void FemElasticDDxNeoHookean(
    const fmat3& F,
    const fmat3& E,
    const fmat3& A,
    const float& Volume,
    const float& lambda,
    const float& mu,
    float& W,
    fvec3& ddx0,
    fvec3& ddx1,
    fvec3& ddx2,
    fvec3& ddx3)
{

	float J	   = F.det();
	float logJ = std::log(J);
	if (J < 0.00001)
		logJ = 0.0;
	W = Volume * (0.5 * mu * (F.sqlength() - 3) - mu * logJ + 0.5 * lambda * logJ * logJ);

	if (std::abs(F.det()) < 0.00001) {
		W = 0.0;
		return;
	}

	fmat3 P = Volume * (mu * F - mu * (F.inverse()).transpose() + lambda * logJ * (F.inverse()).transpose());

	fmat3 PAt = P * A.transpose();

	fvec3 dx1 = fvec3(PAt.m[0], PAt.m[3], PAt.m[6]);
	fvec3 dx2 = fvec3(PAt.m[1], PAt.m[4], PAt.m[7]);
	fvec3 dx3 = fvec3(PAt.m[2], PAt.m[5], PAt.m[8]);
	fvec3 dx0 = -(dx1 + dx2 + dx3);

	fmat3 AFinv  = A * F.inverse();
	fmat3 AAt    = A * A.transpose();
	float AAtsum = AAt.m[0] + AAt.m[1] + AAt.m[2] + AAt.m[3] + AAt.m[4] + AAt.m[5] + AAt.m[6] + AAt.m[7] + AAt.m[8];

	fvec3 rowAFinv1 = fvec3(AFinv.m[0], AFinv.m[1], AFinv.m[2]);
	fvec3 rowAFinv2 = fvec3(AFinv.m[3], AFinv.m[4], AFinv.m[5]);
	fvec3 rowAFinv3 = fvec3(AFinv.m[6], AFinv.m[7], AFinv.m[8]);

	fmat3 H1 = Volume * ((lambda + mu - lambda * logJ) * (rowAFinv1.tensorproduct(rowAFinv1)) + mu * AAt.m[0] * fmat3::indentity());
	fmat3 H2 = Volume * ((lambda + mu - lambda * logJ) * (rowAFinv2.tensorproduct(rowAFinv2)) + mu * AAt.m[1] * fmat3::indentity());
	fmat3 H3 = Volume * ((lambda + mu - lambda * logJ) * (rowAFinv3.tensorproduct(rowAFinv3)) + mu * AAt.m[2] * fmat3::indentity());
	fmat3 H0 = Volume * ((lambda + mu - lambda * logJ) * ((rowAFinv1 + rowAFinv2 + rowAFinv3).tensorproduct(rowAFinv1 + rowAFinv2 + rowAFinv3)) + mu * AAtsum * fmat3::indentity());

	ddx0 = (H0 + 0.001 * fmat3::indentity()).inverse() * dx0;
	ddx1 = (H1 + 0.001 * fmat3::indentity()).inverse() * dx1;
	ddx2 = (H2 + 0.001 * fmat3::indentity()).inverse() * dx2;
	ddx3 = (H3 + 0.001 * fmat3::indentity()).inverse() * dx3;
}
*/

void FemElasticDxNeoHookean(
    const fmat2& F,
    const fmat2& E,
    const fmat2& A,
    const float& Volume,
    const float& lambda,
    const float& mu, float& W,
    fvec2& dx0,
    fvec2& dx1,
    fvec2& dx2)
{

	float J	   = F.det();
	float logJ = std::log(J);
	if (J < 0.00001)
		logJ = 0.0;
	W = Volume * (0.5 * mu * (F.sqlength() - 2) - mu * logJ + 0.5 * lambda * logJ * logJ);

	if (std::abs(F.det()) < 0.00001)
		W = 0.0;
	fmat2 P = Volume * (mu * F - mu * (F.inverse()).transpose() + lambda * logJ * (F.inverse()).transpose());

	fmat2 PAt = P * A.transpose();

	dx1 = fvec2(PAt.m[0], PAt.m[2]);
	dx2 = fvec2(PAt.m[1], PAt.m[3]);
	dx0 = -(dx1 + dx2);
}

void FemElasticDxStVenant(
    const fmat2& F,
    const fmat2& E,
    const fmat2& A,
    const float& Volume,
    const float& lambda,
    const float& mu,
    float& W,
    fvec2& dx0,
    fvec2& dx1,
    fvec2& dx2)
{

	W	= Volume * (mu * E.sqlength() + 0.5 * lambda * E.trace() * E.trace());
	fmat2 P = Volume * (2 * mu * F * E + lambda * E.trace() * F);

	fmat2 PAt = P * A.transpose();

	dx1 = fvec2(PAt.m[0], PAt.m[2]);
	dx2 = fvec2(PAt.m[1], PAt.m[3]);
	dx0 = -(dx1 + dx2);
}

void FemElasticDxCoRotational(
    const fmat2& F,
    const fmat2& E,
    const fmat2& A,
    float& omega,
    const float& Volume,
    const float& lambda,
    const float& mu,
    float& W,
    fvec2& dx0,
    fvec2& dx1,
    fvec2& dx2)
{
	omega = ExtractRotation(F, 3, omega);

	fmat2 R(omega);
	fmat2 S	  = R.transpose() * F;
	float trS = S.trace();
	W	  = Volume * (0.5 * mu * ((F - R).sqlength()) + 0.5 * lambda * (trS * trS - 4 * trS + 4));
	fmat2 P	  = Volume * (mu * (F - R) + lambda * (trS - 2) * R);

	fmat2 PAt = P * A.transpose();

	dx1 = fvec2(PAt.m[0], PAt.m[2]);
	dx2 = fvec2(PAt.m[1], PAt.m[3]);
	dx0 = -(dx1 + dx2);
}
