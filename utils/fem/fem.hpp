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
    const fquaternion& q,
    const float& Volume,
    const float& lambda,
    const float& mu,
    float& W,
    fmat3& P)
{

	fmat3 R	  = q.qtomat();
	fmat3 S	  = R.transpose() * F;
	float trS = S.trace();
	W	  = Volume * (0.5 * mu * (F.sqlength() - 2.0 * trS + 3) + 0.5 * lambda * (trS * trS - 6.0 * trS + 9.0));
	P	  = Volume * (mu * (F - R) + lambda * (trS - 3) * R);
}
