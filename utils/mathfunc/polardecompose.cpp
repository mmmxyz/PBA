#include <cmath>

#include "utils/mathfunc/mathfunc.hpp"
#include "utils/mathfunc/polardecompose.hpp"

/*
fquaternion ExtractRotation(const fmat3& A, const uint32_t& numitr, const fquaternion initq)
{

	//std::cout << A << std::endl;

	fquaternion q = initq;

	//std::cout << q << std::endl;

	fvec3 a0 = fvec3(A.m[0], A.m[3], A.m[6]);
	fvec3 a1 = fvec3(A.m[1], A.m[4], A.m[7]);
	fvec3 a2 = fvec3(A.m[2], A.m[5], A.m[8]);

	for (uint32_t i = 0; i < numitr; i++) {

		fmat3 R = q.qtomat();

		fvec3 r0 = fvec3(R.m[0], R.m[3], R.m[6]);
		fvec3 r1 = fvec3(R.m[1], R.m[4], R.m[7]);
		fvec3 r2 = fvec3(R.m[2], R.m[5], R.m[8]);

		fvec3 grad    = a0.cross(r0) + a1.cross(r1) + a2.cross(r2);
		fmat3 Hessian = fmat3::indentity() * (A.transpose() * R).trace() - 0.5 * R * (A.transpose()) - 0.5 * A * (R.transpose());
		//Hessian	      = 2.001 * fmat3::indentity();
		Hessian = Hessian + 0.001 * fmat3::indentity();

		//std::cout << grad << std::endl;
		//std::cout << Hessian << std::endl;

		q = fquaternion(-1.0 * Hessian.inverse() * grad) * q;
		q = q.normalize();

		//std::cout << q.qtomat() << std::endl;
		//std::cout << (A - q.qtomat()).sqlength() << std::endl;

		//std::cout << "-------------------" << std::endl;
	}

	//std::cout << std::endl;
	//std::cout << (A - q.qtomat()).sqlength() << std::endl;

	return q;
}
*/

fquaternion ExtractRotation(const fmat3& A, const uint32_t& numitr, const fquaternion initq)
{
	//std::cout << A << std::endl;

	fquaternion q = initq;

	//std::cout << q << std::endl;

	fvec3 a0 = fvec3(A.m[0], A.m[3], A.m[6]);
	fvec3 a1 = fvec3(A.m[1], A.m[4], A.m[7]);
	fvec3 a2 = fvec3(A.m[2], A.m[5], A.m[8]);

	for (uint32_t i = 0; i < numitr; i++) {

		fmat3 R = q.qtomat();

		fvec3 r0 = fvec3(R.m[0], R.m[3], R.m[6]);
		fvec3 r1 = fvec3(R.m[1], R.m[4], R.m[7]);
		fvec3 r2 = fvec3(R.m[2], R.m[5], R.m[8]);

		fvec3 omega = r0.cross(a0) + r1.cross(a1) + r2.cross(a2);
		omega	    = omega / (std::abs(a0.dot(r0) + a1.dot(r1) + a2.dot(r2)) + 0.000000001);

		q = fquaternion(omega) * q;
		q = q.normalize();

		//std::cout << q.qtomat() << std::endl;
		//std::cout << (A - q.qtomat()).sqlength() << std::endl;
	}
	//std::cout << (A - q.qtomat()).sqlength() << std::endl;
	return q;
}

float ExtractRotation(const fmat2& A, const uint32_t& numitr, const float& initomega)
{
	float omega = initomega;

	fvec2 a0 = fvec2(A.m[0], A.m[2]);
	fvec2 a1 = fvec2(A.m[1], A.m[3]);

	for (uint32_t i = 0; i < numitr; i++) {
		fmat2 R(omega);
		fvec2 r0 = fvec2(R.m[0], R.m[2]);
		fvec2 r1 = fvec2(R.m[1], R.m[3]);

		float domega = r0.cross(a0) + r1.cross(a1);
		domega	     = domega / (std::abs(a0.dot(r0) + a1.dot(r1)) + 0.00000001);

		omega = domega + omega;
	}
	omega;
}
