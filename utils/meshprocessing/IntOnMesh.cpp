
#include "utils/mathfunc/mathfunc.hpp"

void MeshCM(
    fvec3& Cm,
    float& Mass,
    const fvec3* const MVdata,
    const uint32_t& MVsize,
    const uint32_t* const MIlist,
    const uint32_t& Misize,
    const float& rho)
{
	Cm   = fvec3(0.0, 0.0, 0.0);
	Mass = 0.0;

	for (uint32_t i = 0; i < Misize / 3; i++) {

		fvec3 r0 = MVdata[MIlist[3 * i + 0]];
		fvec3 r1 = MVdata[MIlist[3 * i + 1]];
		fvec3 r2 = MVdata[MIlist[3 * i + 2]];

		Cm = Cm + rho * fvec3::STP(r0, r1, r2) * (1.0 / 24.0) * (r0 + r1 + r2);

		Mass += rho * (1.0 / 6.0) * fvec3::STP(r0, r1, r2);
	}

	Cm = Cm / Mass;
}

void MeshInertia(
    fmat3& Inertia,
    const fvec3& cm,
    const fvec3* const MVdata,
    const uint32_t& MVsize,
    const uint32_t* const MIlist,
    const uint32_t& Misize,
    const float& rho)
{

	Inertia = fmat3();

	for (uint32_t i = 0; i < Misize / 3; i++) {

		fvec3 r0 = MVdata[MIlist[3 * i + 0]] - cm;
		fvec3 r1 = MVdata[MIlist[3 * i + 1]] - cm;
		fvec3 r2 = MVdata[MIlist[3 * i + 2]] - cm;

		Inertia = Inertia + rho * fvec3::STP(r0, r1, r2) * (1.0 / 120.0) * (2.0 * (r0.sqlength() + r1.sqlength() + r2.sqlength() + r0.dot(r1) + r1.dot(r2) + r2.dot(r0)) * fmat3::indentity() - 2.0 * (r0.tensorproduct(r0) + r1.tensorproduct(r1) + r2.tensorproduct(r2)) - (r0.tensorproduct(r1) + r1.tensorproduct(r0) + r1.tensorproduct(r2) + r2.tensorproduct(r1) + r2.tensorproduct(r0) + r0.tensorproduct(r2)));
	}
}
