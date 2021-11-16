#include <iostream>
#include <cmath>
#include <cstdint>
#include <vector>
#include <deque>
#include <algorithm>
#include <random>

#include "opengl/visualizer.hpp"
#include "opengl/vertarray.hpp"
#include "opengl/drawobject.hpp"
#include "opengl/renderer3d.hpp"

#include "utils/mathfunc/mathfunc.hpp"
#include "utils/fileloader/OBJLoader.hpp"
#include "utils/collision/geometry.hpp"
#include "utils/collision/RigidObbBvh.hpp"
#include "utils/meshprocessing/MeshConv.hpp"
#include "utils/meshprocessing/IntOnMesh.hpp"

#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>

#include <omp.h>

constexpr uint32_t FPS = 60;
constexpr float dt     = 1.0 / FPS;
//glfwSwapInterval(1); なので60FPS以下になる。
//これはモニタやGPUに依るかもしれない。

constexpr float CoFs = 0.20f; //静止摩擦係数
constexpr float CoFk = 0.10f; //動摩擦係数
constexpr float CoE  = 0.8f;  //反発係数

class rigidbody {
    public:
	fvec3* PMeshVert;
	uint32_t PMeshVertSize;
	uint32_t* PMeshIlist;
	uint32_t PMeshIlistSize;

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

	float rho;

	fquaternion rotq;
	fvec3 cm;
	float mass;
	fmat3 Inertia;

	fvec3 omega;
	fvec3 velocity;

	RigidObbBvh Bvh;

	rigidbody(const char* Rfilename, const char* Pfilename, const fvec3& offset, const float& scale, const float& rho = 1.0)
	    : tva()
	    , lva()
	    , Bvh(this->rotq, this->cm)
	{

		LoadOBJtoPhysicsTriangleMesh(Pfilename, &PMeshVert, PMeshVertSize, &PMeshIlist, PMeshIlistSize, offset, scale * 1.001);

		MeshCM(cm, mass, PMeshVert, PMeshVertSize, PMeshIlist, PMeshIlistSize, rho);
		MeshInertia(Inertia, cm, PMeshVert, PMeshVertSize, PMeshIlist, PMeshIlistSize, rho);

		//std::cout << "Mass : " << mass << std::endl;
		//std::cout << "Center of Mass : " << cm << std::endl;
		//std::cout << "Inertia tensor : " << Inertia << std::endl;

		for (uint32_t i = 0; i < PMeshVertSize; i++)
			PMeshVert[i] = PMeshVert[i] - cm;

		ConvertPTMtoREM(PMeshVert, PMeshVertSize, PMeshIlist, PMeshIlistSize, &lvdata, lvsize, &lilist, lisize);
		lva.resetvertarray(lvsize, lvdata, lisize, lilist);
		lva.settype(0);

		LoadOBJtoRenderTriangleMesh(
		    Rfilename,
		    &tvdata, tvsize, &tilist, tisize, offset - cm, scale);
		tva.resetvertarray(tvsize, tvdata, tisize, tilist);
		tva.settype(2);

		Bvh.ConstructBvh(PMeshVert, PMeshVertSize, PMeshIlist, PMeshIlistSize);
	}
};

struct dvdo {
	fvec3 dV, dO;
};

class Staticcontact {
	fvec3 relr;
	fvec3 normal;
	rigidbody& RBody;
	float d;
	float lambda;

	fvec3& p;
	fquaternion& q;

    public:
	Staticcontact(fvec3 r, fvec3 n, rigidbody& rbody, float d, fvec3& p, fquaternion& q)
	    : relr(r)
	    , normal(n.normalize())
	    , RBody(rbody)
	    , d(d)
	    , p(p)
	    , q(q)
	    , lambda(0.0)
	{
	}

	void projection(const float stiff = 1.0)
	{
		fvec3 r = p + q.qtomat() * relr;
		float C = r.dot(normal) - d;
		if (C > -0.000000)
			return;

		float w1 = 1.0 / RBody.mass;

		fvec3 Rtn   = (q.qtomat().transpose()) * normal;
		fvec3 rxRtn = relr.cross(Rtn);

		float w2 = rxRtn.dot((RBody.Inertia.inverse()) * rxRtn);

		lambda += (-C / (w1 + w2));
		fvec3 dx = (-C / (w1 + w2)) * w1 * normal;
		fvec3 wd = (-C / (w1 + w2)) * (RBody.Inertia.inverse()) * rxRtn;
		fvec3 w	 = q.qtomat() * wd;

		p = p + stiff * dx;
		//q = q + q * (0.5 * fquaternion(wd, 0.0));
		q = q + (0.5 * fquaternion(stiff * w, 0.0) * q);
		q = q.normalize();

		//std::cout << "constraint" << std::endl;
	}

	dvdo restitution(fvec3 pvc, fvec3 pome, fvec3 cvc, fvec3 come)
	{
		fvec3 Rrelr = q.qtomat() * relr;
		fvec3 prevv = pvc + pome.cross(Rrelr);
		fvec3 corrv = cvc + come.cross(Rrelr);

		float C = CoE * std::min(normal.dot(prevv), 0.0f) + normal.dot(corrv);
		if (C > 0.0)
			return { fvec3(0.0), fvec3(0.0) };

		float w1 = 1.0 / RBody.mass;

		fvec3 Rtn   = (q.qtomat().transpose()) * normal;
		fvec3 rxRtn = relr.cross(Rtn);

		float w2 = rxRtn.dot(RBody.Inertia.inverse() * rxRtn);

		fvec3 dv  = (-C / (w1 + w2)) * w1 * normal;
		fvec3 dwd = (-C / (w1 + w2)) * RBody.Inertia.inverse() * rxRtn;
		fvec3 dw  = q.qtomat() * dwd;

		//std::cout << "restitution" << std::endl;
		//std::cout << dv << dw << std::endl;

		return { dv, dw };
	}

	dvdo friction(fvec3 pvc, fvec3 pome, fvec3 cvc, fvec3 come)
	{
		fvec3 Rrelr = q.qtomat() * relr;
		fvec3 prevv = pvc + pome.cross(Rrelr);
		fvec3 corrv = cvc + come.cross(Rrelr);

		fvec3 force = (lambda / (dt * dt)) * normal;

		fvec3 vtan = corrv - normal.dot(corrv) * normal;

		if (vtan.length() < 0.00000001)
			return { fvec3(0.0), 0.0 };

		fvec3 nvtan = vtan.normalize();

		float minusmin = -std::min(dt * CoFk * (lambda / (dt * dt)), vtan.length());

		float w1 = 1.0 / RBody.mass;

		fvec3 Rtnvtan	= (q.qtomat().transpose()) * nvtan;
		fvec3 rxRtnvtan = relr.cross(Rtnvtan);

		float w2 = rxRtnvtan.dot(RBody.Inertia.inverse() * rxRtnvtan);

		fvec3 dv  = (minusmin / (w1 + w2)) * w1 * nvtan;
		fvec3 dwd = (minusmin / (w1 + w2)) * RBody.Inertia.inverse() * rxRtnvtan;
		fvec3 dw  = q.qtomat() * dwd;

		//std::cout << "friction" << std::endl;
		//std::cout << dv << std::endl;
		//std::cout << dw << std::endl
		//	  << std::endl;

		return { dv, dw };
	}
};

struct dvdodvdo {
	fvec3 dv0, dw0, dv1, dw1;
};

class Dynamiccontact {
	rigidbody& rbody0;
	rigidbody& rbody1;

	fvec3& cm0;
	fvec3& cm1;

	fquaternion& q0;
	fquaternion& q1;

    public:
	fvec3 relr0, relr1;
	fvec3 normal; //v1 - v0
	float lambda;
	Dynamiccontact(const fvec3 r0, const fvec3 r1, fvec3 normal, rigidbody& rbody0, rigidbody& rbody1, fvec3& cm0, fvec3& cm1, fquaternion& q0, fquaternion& q1)
	    : relr0(r0)
	    , relr1(r1)
	    , normal(normal)
	    , rbody0(rbody0)
	    , rbody1(rbody1)
	    , cm0(cm0)
	    , cm1(cm1)
	    , q0(q0)
	    , q1(q1)
	    , lambda(0.0)
	{
	}

	void projection(const float stiff = 1.0)
	{

		//std::cout << "projection" << std::endl;

		fvec3 r0 = cm0 + q0.qtomat() * relr0;
		fvec3 r1 = cm1 + q1.qtomat() * relr1;

		float C = (r0 - r1).dot(normal);
		if (C > 0.01)
			return;

		//std::cout << C << std::endl;

		float w0l    = 1.0 / rbody0.mass;
		fvec3 Rtn0   = (q0.qtomat().transpose()) * normal;
		fvec3 rxRtn0 = relr0.cross(Rtn0);
		float w0r    = rxRtn0.dot((rbody0.Inertia.inverse()) * rxRtn0);

		float w1l    = 1.0 / rbody1.mass;
		fvec3 Rtn1   = (q1.qtomat().transpose()) * -normal;
		fvec3 rxRtn1 = relr1.cross(Rtn1);
		float w1r    = rxRtn1.dot((rbody1.Inertia.inverse()) * rxRtn1);

		lambda += (-C / (w0l + w0r + w1l + w1r));

		fvec3 dx0 = (-C / (w0l + w0r + w1l + w1r)) * w0l * normal;
		fvec3 wd0 = (-C / (w0l + w0r + w1l + w1r)) * (rbody0.Inertia.inverse()) * rxRtn0;

		fvec3 dx1 = (-C / (w0l + w0r + w1l + w1r)) * w1l * -normal;
		fvec3 wd1 = (-C / (w0l + w0r + w1l + w1r)) * (rbody1.Inertia.inverse()) * rxRtn1;

		fvec3 w0 = q0.qtomat() * wd0;
		fvec3 w1 = q1.qtomat() * wd1;

		cm0 = cm0 + stiff * dx0;
		cm1 = cm1 + stiff * dx1;

		q0 = q0 + (0.5 * fquaternion(stiff * w0, 0.0) * q0);
		q0 = q0.normalize();
		q1 = q1 + (0.5 * fquaternion(stiff * w1, 0.0) * q1);
		q1 = q1.normalize();

		//std::cout << dx0 << dx1 << w0 << w1 << std::endl;
	}

	dvdodvdo restitution(
	    fvec3 pvc0, fvec3 pome0, fvec3 cvc0, fvec3 come0,
	    fvec3 pvc1, fvec3 pome1, fvec3 cvc1, fvec3 come1)
	{
		fvec3 Rrelr0 = q0.qtomat() * relr0;
		fvec3 prevv0 = pvc0 + pome0.cross(Rrelr0);
		fvec3 corrv0 = cvc0 + come0.cross(Rrelr0);

		fvec3 Rrelr1 = q1.qtomat() * relr1;
		fvec3 prevv1 = pvc1 + pome1.cross(Rrelr1);
		fvec3 corrv1 = cvc1 + come1.cross(Rrelr1);

		fvec3 prevv = -prevv1 + prevv0;
		fvec3 corrv = -corrv1 + corrv0;

		//std::cout << prevv << " " << corrv << std::endl;
		float C = CoE * std::min(normal.dot(prevv), 0.0f) + normal.dot(corrv);
		if (C > 0.0)
			return { fvec3(0.0), fvec3(0.0), fvec3(0.0), fvec3(0.0) };

		float w0l = 1.0 / rbody0.mass;
		float w1l = 1.0 / rbody1.mass;

		fvec3 Rtn0   = (q0.qtomat().transpose()) * normal;
		fvec3 rxRtn0 = relr0.cross(Rtn0);

		fvec3 Rtn1   = (q1.qtomat().transpose()) * -normal;
		fvec3 rxRtn1 = relr1.cross(Rtn1);

		float w0r = rxRtn0.dot(rbody0.Inertia.inverse() * rxRtn0);
		float w1r = rxRtn1.dot(rbody1.Inertia.inverse() * rxRtn1);

		fvec3 dv0  = (-C / (w0l + w0r + w1l + w1r)) * w0l * normal;
		fvec3 dwd0 = (-C / (w0l + w0r + w1l + w1r)) * rbody0.Inertia.inverse() * rxRtn0;
		fvec3 dw0  = q0.qtomat() * dwd0;

		fvec3 dv1  = (-C / (w0l + w0r + w1l + w1r)) * w1l * -normal;
		fvec3 dwd1 = (-C / (w0l + w0r + w1l + w1r)) * rbody1.Inertia.inverse() * rxRtn1;
		fvec3 dw1  = q1.qtomat() * dwd1;

		//std::cout << rbody0.mass * dv0 + rbody1.mass * dv1 << std::endl;
		//std::cout << q0.qtomat() * (rbody0.Inertia * dwd0) + q1.qtomat() * (rbody1.Inertia * dwd1) << std::endl;

		//std::cout << "restitution" << std::endl;

		return { dv0, dw0, dv1, dw1 };
	}

	dvdodvdo friction(
	    fvec3 pvc0, fvec3 pome0, fvec3 cvc0, fvec3 come0,
	    fvec3 pvc1, fvec3 pome1, fvec3 cvc1, fvec3 come1)
	{
		fvec3 Rrelr0 = q0.qtomat() * relr0;
		fvec3 prevv0 = pvc0 + pome0.cross(Rrelr0);
		fvec3 corrv0 = cvc0 + come0.cross(Rrelr0);

		fvec3 Rrelr1 = q1.qtomat() * relr1;
		fvec3 prevv1 = pvc1 + pome1.cross(Rrelr1);
		fvec3 corrv1 = cvc1 + come1.cross(Rrelr1);

		fvec3 prevv = -prevv1 + prevv0;
		fvec3 corrv = -corrv1 + corrv0;

		fvec3 force = (lambda / (dt * dt)) * normal;

		fvec3 vtan = corrv - normal.dot(corrv) * normal;

		if (vtan.length() < 0.00000001)
			return { fvec3(0.0), fvec3(0.0), fvec3(0.0), fvec3(0.0) };

		fvec3 nvtan    = vtan.normalize();
		float minusmin = -std::min(dt * CoFk * (lambda / (dt * dt)), vtan.length());

		float w0l = 1.0 / rbody0.mass;
		float w1l = 1.0 / rbody1.mass;

		fvec3 Rtn0   = (q0.qtomat().transpose()) * nvtan;
		fvec3 rxRtn0 = relr0.cross(Rtn0);

		fvec3 Rtn1   = (q1.qtomat().transpose()) * -nvtan;
		fvec3 rxRtn1 = relr1.cross(Rtn1);

		float w0r = rxRtn0.dot(rbody0.Inertia.inverse() * rxRtn0);
		float w1r = rxRtn1.dot(rbody1.Inertia.inverse() * rxRtn1);

		fvec3 dv0  = (minusmin / (w0l + w0r + w1l + w1r)) * w0l * nvtan;
		fvec3 dwd0 = (minusmin / (w0l + w0r + w1l + w1r)) * rbody0.Inertia.inverse() * rxRtn0;
		fvec3 dw0  = q0.qtomat() * dwd0;

		fvec3 dv1  = (minusmin / (w0l + w0r + w1l + w1r)) * w1l * -nvtan;
		fvec3 dwd1 = (minusmin / (w0l + w0r + w1l + w1r)) * rbody1.Inertia.inverse() * rxRtn1;
		fvec3 dw1  = q1.qtomat() * dwd1;

		return { dv0, dw0, dv1, dw1 };
	}
};

class Staticfriction {
	fvec3 relr0, relr1;
	fvec3 normal; //v1 - v0
	rigidbody& rbody0;
	rigidbody& rbody1;
	float forcelambda;

	fvec3* cm0;
	fvec3* cm1;

	fquaternion* q0;
	fquaternion* q1;

    public:
	float lambda;
	Staticfriction(const fvec3 r0, const fvec3 r1, fvec3 normal, rigidbody& rbody0, rigidbody& rbody1, const float forcelambda)
	    : relr0(r0)
	    , relr1(r1)
	    , normal(normal)
	    , rbody0(rbody0)
	    , rbody1(rbody1)
	    , cm0(cm0)
	    , cm1(cm1)
	    , q0(q0)
	    , q1(q1)
	    , forcelambda(forcelambda)
	    , lambda(0.0)
	{
	}

	void setCorrd(fvec3* cm0, fvec3* cm1, fquaternion* q0, fquaternion* q1)
	{
		this->cm0 = cm0;
		this->cm1 = cm1;
		this->q0  = q0;
		this->q1  = q1;
	}

	void projection(const float stiff = 1.0)
	{
		fvec3 r0 = *cm0 + q0->qtomat() * relr0;
		fvec3 r1 = *cm1 + q1->qtomat() * relr1;

		fvec3 tr01 = (r0 - r1) - ((r0 - r1).dot(normal)) * normal;

		if (tr01.sqlength() < 0.0000001)
			return;

		fvec3 nr01 = tr01.normalize();

		float alphat = 1.0;
		float C	     = std::min(CoFs * (forcelambda / (dt * dt)), tr01.length());

		float w0l = 1.0 / rbody0.mass;
		float w1l = 1.0 / rbody1.mass;

		fvec3 Rtn0   = (q0->qtomat().transpose()) * nr01;
		fvec3 rxRtn0 = relr0.cross(Rtn0);

		fvec3 Rtn1   = (q1->qtomat().transpose()) * -nr01;
		fvec3 rxRtn1 = relr1.cross(Rtn1);

		float w0r = rxRtn0.dot(rbody0.Inertia.inverse() * rxRtn0);
		float w1r = rxRtn1.dot(rbody1.Inertia.inverse() * rxRtn1);

		lambda += (-C - alphat * lambda / (w0l + w0r + w1l + w1r + alphat * lambda));
		fvec3 dx0 = (-C - alphat * lambda / (w0l + w0r + w1l + w1r + alphat * lambda)) * w0l * nr01;
		fvec3 wd0 = (-C - alphat * lambda / (w0l + w0r + w1l + w1r + alphat * lambda)) * (rbody0.Inertia.inverse()) * rxRtn0;
		fvec3 dx1 = (-C - alphat * lambda / (w0l + w0r + w1l + w1r + alphat * lambda)) * w1l * -nr01;
		fvec3 wd1 = (-C - alphat * lambda / (w0l + w0r + w1l + w1r + alphat * lambda)) * (rbody1.Inertia.inverse()) * rxRtn1;

		fvec3 w0 = q0->qtomat() * wd0;
		fvec3 w1 = q1->qtomat() * wd1;

		std::cout << dx0 << " " << dx1 << std::endl;

		*cm0 = *cm0 + stiff * dx0;
		*cm1 = *cm1 + stiff * dx1;

		*q0 = *q0 + (0.5 * fquaternion(stiff * w0, 0.0) * *q0);
		*q0 = q0->normalize();
		*q1 = *q1 + (0.5 * fquaternion(stiff * w1, 0.0) * *q1);
		*q1 = q1->normalize();
	}
};

namespace Physics {
std::vector<rigidbody> rbodyvec;
std::vector<std::deque<Staticfriction>> StaticFrictionVec;

void SetRbody(rigidbody&& r)
{
	rbodyvec.emplace_back(r);
}

void timestep()
{

	uint32_t rsize = rbodyvec.size();

	std::vector<fvec3> tempp(rsize);
	std::vector<fquaternion> tempq(rsize);

	for (uint32_t i = 0; i < rsize; i++) {
		fvec3 velocity	 = rbodyvec[i].velocity + dt * fvec3(0.0, -9.8, 0.0);
		tempp[i]	 = rbodyvec[i].cm + dt * velocity;
		fmat3 Inertia	 = (rbodyvec[i].rotq.qtomat()) * (rbodyvec[i].Inertia) * (rbodyvec[i].rotq.qtomat().transpose());
		fmat3 invInertia = Inertia.inverse();
		fvec3 hogeomega	 = rbodyvec[i].omega - dt * invInertia * (rbodyvec[i].omega.cross(Inertia * rbodyvec[i].omega));
		tempq[i]	 = rbodyvec[i].rotq + (0.5 * dt * fquaternion(hogeomega, 0.0) * rbodyvec[i].rotq);
		tempq[i]	 = tempq[i].normalize();
	}

	std::vector<std::deque<Staticcontact>> Scontactlist(rsize);

#pragma omp parallel for
	for (uint32_t i = 0; i < rsize; i++) {
		for (uint32_t j = 0; j < rbodyvec[i].PMeshVertSize; j++) {

			fvec3 relr  = rbodyvec[i].PMeshVert[j];
			fvec3 edgep = tempp[i] + tempq[i].qtomat() * relr;

			if (edgep.x < -12.0f)
				Scontactlist[i].emplace_back(Staticcontact(relr, fvec3(1.0, 0.0, 0.0), rbodyvec[i], -12.0, tempp[i], tempq[i]));
			if (edgep.x > 12.0f)
				Scontactlist[i].emplace_back(Staticcontact(relr, fvec3(-1.0, 0.0, 0.0), rbodyvec[i], -12.0, tempp[i], tempq[i]));

			if (edgep.y < -12.0f)
				Scontactlist[i].emplace_back(Staticcontact(relr, fvec3(0.0, 1.0, 0.0), rbodyvec[i], -12.0, tempp[i], tempq[i]));
			if (edgep.y > 12.0f)
				Scontactlist[i].emplace_back(Staticcontact(relr, fvec3(0.0, -1.0, 0.0), rbodyvec[i], -12.0, tempp[i], tempq[i]));

			if (edgep.z < -12.0f)
				Scontactlist[i].emplace_back(Staticcontact(relr, fvec3(0.0, 0.0, 1.0), rbodyvec[i], -12.0, tempp[i], tempq[i]));
			if (edgep.z > 12.0f)
				Scontactlist[i].emplace_back(Staticcontact(relr, fvec3(0.0, 0.0, -1.0), rbodyvec[i], -12.0, tempp[i], tempq[i]));
		}
	}

	std::vector<vec2<uint32_t>> ContactPair;
	for (uint32_t i = 0; i < rsize; i++) {
		for (uint32_t j = i + 1; j < rsize; j++) {
			ContactPair.emplace_back(vec2<uint32_t>(i, j));
		}
	}

	uint32_t PairSize = ContactPair.size();

	std::vector<std::deque<TriangleInd>> PCList(PairSize);

#pragma omp parallel for
	for (uint32_t i = 0; i < PairSize; i++) {
		uint32_t Ind0 = ContactPair[i].x;
		uint32_t Ind1 = ContactPair[i].y;
		RigidObbBvhPotentialCollisionPair(rbodyvec[Ind0].Bvh, rbodyvec[Ind1].Bvh, PCList[i], tempp[Ind0], tempp[Ind1], tempq[Ind0], tempq[Ind1]);
	}

	//std::cout << PCList.size() << std::endl;

	//StaticFrictionVec.resize(((rsize - 1) * rsize) / 2);
	//for (uint32_t i = 0; i < ((rsize - 1) * rsize) / 2; i++) {
	//	uint32_t Ind0 = ContactPair[i].x;
	//	uint32_t Ind1 = ContactPair[i].y;
	//	for (auto& con : StaticFrictionVec[i]) {
	//		con.setCorrd(&(tempp[Ind0]), &(tempp[Ind1]), &(tempq[Ind0]), &(tempq[Ind1]));
	//	}
	//}

	std::vector<std::deque<Dynamiccontact>> Dcontactlist(PairSize);
	std::vector<std::deque<Dynamiccontact>> CurrentDcontactlist(PairSize);

	for (uint32_t i = 0; i < 15; i++) {

		//float stiff = std::max((i / 80.0), 1.0);
		float stiff = 1.0 / 15;

		//for (uint32_t i = 0; i < ((rsize - 1) * rsize) / 2; i++) {
		//	uint32_t Ind0 = ContactPair[i].x;
		//	uint32_t Ind1 = ContactPair[i].y;
		//	for (auto& con : StaticFrictionVec[i]) {
		//		con.projection(stiff);
		//	}
		//}

		for (uint32_t i = 0; i < rsize; i++) {
			for (auto& con : Scontactlist[i]) {
				con.projection(stiff);
			}
		}

		if (i % 5 == 0) {

#pragma omp parallel for
			for (uint32_t j = 0; j < PairSize; j++) {
				std::move(CurrentDcontactlist[j].begin(), CurrentDcontactlist[j].end(), std::back_inserter(Dcontactlist[j]));
				CurrentDcontactlist[j].clear();

				uint32_t Ind0 = ContactPair[j].x;
				uint32_t Ind1 = ContactPair[j].y;

				for (const auto& x : PCList[j]) {
					fvec3 v00 = tempp[Ind0] + tempq[Ind0].qtomat() * rbodyvec[Ind0].PMeshVert[x.v00];
					fvec3 v01 = tempp[Ind0] + tempq[Ind0].qtomat() * rbodyvec[Ind0].PMeshVert[x.v01];
					fvec3 v02 = tempp[Ind0] + tempq[Ind0].qtomat() * rbodyvec[Ind0].PMeshVert[x.v02];

					fvec3 v10 = tempp[Ind1] + tempq[Ind1].qtomat() * rbodyvec[Ind1].PMeshVert[x.v10];
					fvec3 v11 = tempp[Ind1] + tempq[Ind1].qtomat() * rbodyvec[Ind1].PMeshVert[x.v11];
					fvec3 v12 = tempp[Ind1] + tempq[Ind1].qtomat() * rbodyvec[Ind1].PMeshVert[x.v12];

					ContactFeature CF1;
					ContactFeature CF2;
					uint32_t hogec = contactontriangle(v00, v01, v02, v10, v11, v12, CF1, CF2);
					if (hogec == 1) {
						fvec3 r0rel = (tempq[Ind0].qtomat().transpose()) * (CF1.r0rel - tempp[Ind0]);
						fvec3 r1rel = (tempq[Ind1].qtomat().transpose()) * (CF1.r1rel - tempp[Ind1]);
						CurrentDcontactlist[j].push_back(Dynamiccontact(r0rel, r1rel, CF1.normal, rbodyvec[Ind0], rbodyvec[Ind1], tempp[Ind0], tempp[Ind1], tempq[Ind0], tempq[Ind1]));
					}
				}
			}
		}

		for (uint32_t j = 0; j < PairSize; j++)
			if (i % 2 == 1)
				for (auto& con : CurrentDcontactlist[j]) {
					con.projection(stiff);
				}
			else
				for (std::deque<Dynamiccontact>::reverse_iterator it = CurrentDcontactlist[j].rbegin(); it != CurrentDcontactlist[j].rend(); ++it)
					(*it).projection(stiff);
	}

	std::vector<fvec3> prevvelocity(rsize);
	std::vector<fvec3> prevomega(rsize);

	for (uint32_t i = 0; i < rsize; i++) {
		prevvelocity[i]	     = rbodyvec[i].velocity;
		rbodyvec[i].velocity = (tempp[i] - rbodyvec[i].cm) / dt;
		rbodyvec[i].cm	     = tempp[i];

		fquaternion prevq = rbodyvec[i].rotq;
		rbodyvec[i].rotq  = tempq[i];
		fvec3 prevomega	  = rbodyvec[i].omega;

		fquaternion dq	  = tempq[i] * (prevq.inverse());
		rbodyvec[i].omega = (2.0 / dt) * dq.getv();
	}

	//反発・摩擦は逐次的に計算する。並列にdv,dwを計算してから総和を取ると、速度が0にならず、物体が静止しない。
	for (uint32_t i = 0; i < 15; i++) {

		float stiff = 0.6;

		for (uint32_t j = 0; j < PairSize; j++) {
			uint32_t Ind0 = ContactPair[j].x;
			uint32_t Ind1 = ContactPair[j].y;
			for (auto& con : Dcontactlist[j]) {
				auto [tdvx, tdwx, tdvy, tdwy] = con.restitution(
				    prevvelocity[Ind0], prevomega[Ind0], rbodyvec[Ind0].velocity, rbodyvec[Ind0].omega,
				    prevvelocity[Ind1], prevomega[Ind1], rbodyvec[Ind1].velocity, rbodyvec[Ind1].omega);

				rbodyvec[Ind0].velocity = rbodyvec[Ind0].velocity + stiff * tdvx;
				rbodyvec[Ind0].omega	= rbodyvec[Ind0].omega + stiff * tdwx;
				rbodyvec[Ind1].velocity = rbodyvec[Ind1].velocity + stiff * tdvy;
				rbodyvec[Ind1].omega	= rbodyvec[Ind1].omega + stiff * tdwy;
			}
		}

		for (uint32_t j = 0; j < rsize; j++) {
			for (auto& con : Scontactlist[j]) {
				auto [tdv, tdw]	     = con.restitution(prevvelocity[j], prevomega[j], rbodyvec[j].velocity, rbodyvec[j].omega);
				rbodyvec[j].velocity = rbodyvec[j].velocity + stiff * tdv;
				rbodyvec[j].omega    = rbodyvec[j].omega + stiff * tdw;
			}
		}
	}

	for (uint32_t j = 0; j < rsize; j++) {
		for (auto& con : Scontactlist[j]) {
			auto [tdv, tdw]	     = con.friction(prevvelocity[j], prevomega[j], rbodyvec[j].velocity, rbodyvec[j].omega);
			rbodyvec[j].velocity = rbodyvec[j].velocity + tdv;
			rbodyvec[j].omega    = rbodyvec[j].omega + tdw;
		}
	}
	for (uint32_t j = 0; j < PairSize; j++) {
		uint32_t Ind0 = ContactPair[j].x;
		uint32_t Ind1 = ContactPair[j].y;
		for (auto& con : Dcontactlist[j]) {
			auto [tdvx, tdwx, tdvy, tdwy] = con.friction(
			    prevvelocity[Ind0], prevomega[Ind0], rbodyvec[Ind0].velocity, rbodyvec[Ind0].omega,
			    prevvelocity[Ind1], prevomega[Ind1], rbodyvec[Ind1].velocity, rbodyvec[Ind1].omega);

			rbodyvec[Ind0].velocity = rbodyvec[Ind0].velocity + tdvx;
			rbodyvec[Ind0].omega	= rbodyvec[Ind0].omega + tdwx;
			rbodyvec[Ind1].velocity = rbodyvec[Ind1].velocity + tdvy;
			rbodyvec[Ind1].omega	= rbodyvec[Ind1].omega + tdwy;
		}
	}

	for (auto& scon : StaticFrictionVec)
		scon.clear();

	//for (uint32_t i = 0; i < ((rsize - 1) * rsize) / 2; i++) {
	//	uint32_t Ind0 = ContactPair[i].x;
	//	uint32_t Ind1 = ContactPair[i].y;
	//	for (auto& con : Dcontactlist[i]) {
	//		fvec3 r0 = rbodyvec[Ind0].cm + rbodyvec[Ind0].rotq.qtomat() * con.relr0;
	//		fvec3 r1 = rbodyvec[Ind1].cm + rbodyvec[Ind1].rotq.qtomat() * con.relr1;
	//		fvec3 v0 = rbodyvec[Ind0].velocity + rbodyvec[Ind0].omega.cross(con.relr0);
	//		fvec3 v1 = rbodyvec[Ind1].velocity + rbodyvec[Ind1].omega.cross(con.relr1);

	//		if (std::abs((r0 - r1).dot(con.normal)) < 0.00001 && std::abs((v0 - v1).dot(con.normal)) < 0.01 && con.lambda > 0.001) {
	//			StaticFrictionVec[i].emplace_back(Staticfriction(con.relr0, con.relr1, con.normal, rbodyvec[Ind0], rbodyvec[Ind0], con.lambda));
	//		}
	//	}
	//}
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

	omp_set_num_threads(20);

	Renderer3D::Init();

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

	Physics::rbodyvec.reserve(10);

	//いろいろ
	/*
	uint32_t objectsize = 6;
	Physics::rbodyvec.emplace_back("../../../resource/Bunny.obj", "../../../resource/LightBunny.obj", fvec3(4.8, -5.1, -0.8), 4.00);
	Physics::rbodyvec[0].rotq     = fquaternion(fvec3(-0.000 * 3.1415, 2.00 * 3.1415, 0.0 * 3.1415));
	Physics::rbodyvec[0].omega    = fvec3(.0, .0, .0);
	Physics::rbodyvec[0].velocity = fvec3(.0, .0, .0);
	Physics::rbodyvec[0].tva.setcolor(0.4, 0.4, 0.8, 1.0);
	Physics::rbodyvec[0].tva.settype(2);
	Physics::rbodyvec[0].lva.setcolor(0.8, 0.8, 0.8, 1.0);

	Physics::rbodyvec.emplace_back("../../../resource/Cube.obj", "../../../resource/Cube.obj", fvec3(0.0, -8.0, 0.0), 4.00, 1.00);
	Physics::rbodyvec[1].rotq     = fquaternion(fvec3(.0 * 3.1415, 0.50 * 3.1415, 0.0));
	Physics::rbodyvec[1].omega    = fvec3(3.0, 0.0, 0.0);
	Physics::rbodyvec[1].velocity = fvec3(0.0, 0.0, -3.0);
	Physics::rbodyvec[1].tva.setcolor(0.8, 0.4, 0.3, 1.0);
	Physics::rbodyvec[1].tva.settype(1);
	Physics::rbodyvec[1].lva.setcolor(0.8, 0.8, 0.8, 1.0);

	Physics::rbodyvec.emplace_back("../../../resource/TeaPot.obj", "../../../resource/LightTeaPot.obj", fvec3(0.0, 1.0, 0.0), 0.04, 1.00);
	Physics::rbodyvec[2].rotq     = fquaternion(fvec3(.0 * 3.1415, 0.50 * 3.1415, 0.0));
	Physics::rbodyvec[2].omega    = fvec3(3.0, 3.0, 0.0);
	Physics::rbodyvec[2].velocity = fvec3(3.0, 0.0, -3.0);
	Physics::rbodyvec[2].tva.setcolor(0.8, 0.4, 0.3, 1.0);
	Physics::rbodyvec[2].tva.settype(1);
	Physics::rbodyvec[2].lva.setcolor(0.8, 0.8, 0.8, 1.0);

	Physics::rbodyvec.emplace_back("../../../resource/Serapis.obj", "../../../resource/LightSerapis.obj", fvec3(0.0, 10.0, 0.0), 0.18, 8.00);
	Physics::rbodyvec[3].rotq     = fquaternion(fvec3(.0 * 3.1415, 0.50 * 3.1415, 0.0));
	Physics::rbodyvec[3].omega    = fvec3(0.0, -8.0, 0.0);
	Physics::rbodyvec[3].velocity = fvec3(3.0, 3.0, -3.0);
	Physics::rbodyvec[3].tva.setcolor(0.8, 0.4, 0.3, 1.0);
	Physics::rbodyvec[3].tva.settype(1);
	Physics::rbodyvec[3].lva.setcolor(0.8, 0.8, 0.8, 1.0);

	Physics::rbodyvec.emplace_back("../../../resource/Dragon.obj", "../../../resource/LightDragon.obj", fvec3(0.0, 2.0, 6.0), 7.00, 1.00);
	Physics::rbodyvec[4].rotq     = fquaternion(fvec3(.0 * 3.1415, 0.50 * 3.1415, 0.0));
	Physics::rbodyvec[4].omega    = fvec3(3.0, 3.0, 3.0);
	Physics::rbodyvec[4].velocity = fvec3(0.0, 3.0, 10.0);
	Physics::rbodyvec[4].tva.setcolor(0.8, 0.4, 0.3, 1.0);
	Physics::rbodyvec[4].tva.settype(1);
	Physics::rbodyvec[4].lva.setcolor(0.8, 0.8, 0.8, 1.0);

	Physics::rbodyvec.emplace_back("../../../resource/Bunny.obj", "../../../resource/LightBunny.obj", fvec3(-3.0, 2.0, 0.0), 2.00, 1.00);
	Physics::rbodyvec[5].rotq     = fquaternion(fvec3(.0 * 3.1415, 0.50 * 3.1415, 0.0));
	Physics::rbodyvec[5].omega    = fvec3(3.0, 2.0, 2.0);
	Physics::rbodyvec[5].velocity = fvec3(3.0, 30.0, -10.0);
	Physics::rbodyvec[5].tva.setcolor(0.8, 0.4, 0.3, 1.0);
	Physics::rbodyvec[5].tva.settype(1);
	Physics::rbodyvec[5].lva.setcolor(0.8, 0.8, 0.8, 1.0);
	*/

	/////////////////////////////////////

	//cubeから滑り落ちるcubeとbunnyのstack
	uint32_t objectsize = 5;
	Physics::rbodyvec.emplace_back("../../../resource/Cube.obj", "../../../resource/Cube.obj", fvec3(-5.0, -8.0, 0.0), 4.00, 1.00);
	Physics::rbodyvec[0].rotq     = fquaternion(fvec3(.0 * 3.1415, 2.00 * 3.1415, 0.0));
	Physics::rbodyvec[0].omega    = fvec3(0.0, 0.0, 0.0);
	Physics::rbodyvec[0].velocity = fvec3(0.0, 0.0, 0.0);
	Physics::rbodyvec[0].tva.setcolor(0.8, 0.4, 0.3, 1.0);
	Physics::rbodyvec[0].tva.settype(2);
	Physics::rbodyvec[0].lva.setcolor(0.8, 0.8, 0.8, 1.0);

	Physics::rbodyvec.emplace_back("../../../resource/Cube.obj", "../../../resource/Cube.obj", fvec3(-1.0, -1.0, 0.0), 3.00, 1.00);
	Physics::rbodyvec[1].rotq     = fquaternion(fvec3(.0 * 3.1415, 2.00 * 3.1415, 0.0 * 3.1415));
	Physics::rbodyvec[1].omega    = fvec3(0.0, 0.0, 0.0);
	Physics::rbodyvec[1].velocity = fvec3(0.0, 0.0, 0.0);
	Physics::rbodyvec[1].tva.setcolor(0.6, 0.1, 0.7, 1.0);
	Physics::rbodyvec[1].tva.settype(2);
	Physics::rbodyvec[1].lva.setcolor(0.8, 0.8, 0.8, 1.0);

	Physics::rbodyvec.emplace_back("../../../resource/Cube.obj", "../../../resource/Cube.obj", fvec3(-1.0, 4.00, 0.0), 2.00, 1.00);
	Physics::rbodyvec[2].rotq     = fquaternion(fvec3(.0 * 3.1415, 2.00 * 3.1415, 0.0));
	Physics::rbodyvec[2].omega    = fvec3(0.0, 0.0, 0.0);
	Physics::rbodyvec[2].velocity = fvec3(0.0, 0.0, 1.0);
	Physics::rbodyvec[2].tva.setcolor(0.2, 0.8, 0.3, 1.0);
	Physics::rbodyvec[2].tva.settype(2);
	Physics::rbodyvec[2].lva.setcolor(0.8, 0.8, 0.8, 1.0);

	Physics::rbodyvec.emplace_back("../../../resource/Dragon.obj", "../../../resource/LightDragon.obj", fvec3(0.0, 8.2, 0.5), 8.00, 7.00);
	Physics::rbodyvec[3].rotq     = fquaternion(fvec3(.0 * 3.1415, 2.00 * 3.1415, 0.0));
	Physics::rbodyvec[3].omega    = fvec3(0.0, 0.0, 0.0);
	Physics::rbodyvec[3].velocity = fvec3(0.0, 0.0, 1.0);
	Physics::rbodyvec[3].tva.setcolor(0.3, 0.8, 0.8, 1.0);
	Physics::rbodyvec[3].tva.settype(1);
	Physics::rbodyvec[3].lva.setcolor(0.8, 0.8, 0.8, 1.0);

	Physics::rbodyvec.emplace_back("../../../resource/Bunny.obj", "../../../resource/LightBunny.obj", fvec3(8.5, 6.0, -0.5), 2.00, 6.00);
	Physics::rbodyvec[4].rotq     = fquaternion(fvec3(.0 * 3.1415, 0.25 * 3.1415, 0.0));
	Physics::rbodyvec[4].omega    = fvec3(4.0, 4.0, 5.0);
	Physics::rbodyvec[4].velocity = fvec3(-12.0, 3.0, 1.0);
	Physics::rbodyvec[4].tva.setcolor(0.8, 0.3, 0.3, 1.0);
	Physics::rbodyvec[4].tva.settype(1);
	Physics::rbodyvec[4].lva.setcolor(0.8, 0.8, 0.8, 1.0);

	////////////////////////////////////

	//cubeのstack
	//重心が支点の内側にあっても倒壊する
	/*
	uint32_t objectsize = 4;
	Physics::rbodyvec.emplace_back("../../../resource/Cube.obj", "../../../resource/Cube.obj", fvec3(-5.0, -8.0, 0.0), 4.00, 1.00);
	Physics::rbodyvec[0].rotq     = fquaternion(fvec3(.0 * 3.1415, 2.00 * 3.1415, 0.0));
	Physics::rbodyvec[0].omega    = fvec3(0.0, 0.0, 0.0);
	Physics::rbodyvec[0].velocity = fvec3(0.0, 0.0, 0.0);
	Physics::rbodyvec[0].tva.setcolor(0.8, 0.4, 0.3, 1.0);
	Physics::rbodyvec[0].tva.settype(2);
	Physics::rbodyvec[0].lva.setcolor(0.8, 0.8, 0.8, 1.0);

	Physics::rbodyvec.emplace_back("../../../resource/Cube.obj", "../../../resource/Cube.obj", fvec3(-1.2, -2.0, 1.0), 2.00, 1.00);
	Physics::rbodyvec[1].rotq     = fquaternion(fvec3(.0 * 3.1415, 2.00 * 3.1415, 0.0 * 3.1415));
	Physics::rbodyvec[1].omega    = fvec3(0.0, 0.0, 0.0);
	Physics::rbodyvec[1].velocity = fvec3(0.0, 0.0, 0.0);
	Physics::rbodyvec[1].tva.setcolor(0.6, 0.1, 0.7, 1.0);
	Physics::rbodyvec[1].tva.settype(2);
	Physics::rbodyvec[1].lva.setcolor(0.8, 0.8, 0.8, 1.0);

	Physics::rbodyvec.emplace_back("../../../resource/Cube.obj", "../../../resource/Cube.obj", fvec3(-1.2, 3.0, 1.0), 3.00, 1.00);
	Physics::rbodyvec[2].rotq     = fquaternion(fvec3(.0 * 3.1415, 2.00 * 3.1415, 0.0));
	Physics::rbodyvec[2].omega    = fvec3(0.0, 0.0, 0.0);
	Physics::rbodyvec[2].velocity = fvec3(0.0, 0.0, 1.0);
	Physics::rbodyvec[2].tva.setcolor(0.2, 0.8, 0.3, 1.0);
	Physics::rbodyvec[2].tva.settype(2);
	Physics::rbodyvec[2].lva.setcolor(0.8, 0.8, 0.8, 1.0);

	Physics::rbodyvec.emplace_back("../../../resource/Cube.obj", "../../../resource/Cube.obj", fvec3(-1.2, 8.0, 1.0), 2.00, 1.00);
	Physics::rbodyvec[3].rotq     = fquaternion(fvec3(.0 * 3.1415, 2.00 * 3.1415, 0.0));
	Physics::rbodyvec[3].omega    = fvec3(0.0, 0.0, 0.0);
	Physics::rbodyvec[3].velocity = fvec3(0.0, 0.0, 1.0);
	Physics::rbodyvec[3].tva.setcolor(0.8, 0.8, 0.8, 1.0);
	Physics::rbodyvec[3].tva.settype(2);
	Physics::rbodyvec[3].lva.setcolor(0.8, 0.8, 0.8, 1.0);
	*/

	////////////////////////////////////

	//一応安定するcubeとbunnyのstack
	//仮想時間で1分30秒くらいは耐える
	/*
	uint32_t objectsize = 5;
	Physics::rbodyvec.emplace_back("../../../resource/Cube.obj", "../../../resource/Cube.obj", fvec3(0.0, -10.0, 0.0), 2.00, 1.00);
	Physics::rbodyvec[0].rotq     = fquaternion(fvec3(.0 * 3.1415, 2.00 * 3.1415, 0.0));
	Physics::rbodyvec[0].omega    = fvec3(0.0, 0.0, 0.0);
	Physics::rbodyvec[0].velocity = fvec3(0.0, 0.0, 0.0);
	Physics::rbodyvec[0].tva.setcolor(0.8, 0.4, 0.3, 1.0);
	Physics::rbodyvec[0].tva.settype(2);
	Physics::rbodyvec[0].lva.setcolor(0.8, 0.8, 0.8, 1.0);

	Physics::rbodyvec.emplace_back("../../../resource/Cube.obj", "../../../resource/Cube.obj", fvec3(0.0, -6.0, 0.0), 2.00, 1.00);
	Physics::rbodyvec[1].rotq     = fquaternion(fvec3(.0 * 3.1415, 0.25 * 3.1415, 0.0 * 3.1415));
	Physics::rbodyvec[1].omega    = fvec3(0.0, 0.0, 0.0);
	Physics::rbodyvec[1].velocity = fvec3(0.0, 0.0, 0.0);
	Physics::rbodyvec[1].tva.setcolor(0.6, 0.1, 0.7, 1.0);
	Physics::rbodyvec[1].tva.settype(2);
	Physics::rbodyvec[1].lva.setcolor(0.8, 0.8, 0.8, 1.0);

	Physics::rbodyvec.emplace_back("../../../resource/Cube.obj", "../../../resource/Cube.obj", fvec3(0.0, -2.0, 0.0), 2.00, 1.00);
	Physics::rbodyvec[2].rotq     = fquaternion(fvec3(.0 * 3.1415, 2.00 * 3.1415, 0.0));
	Physics::rbodyvec[2].omega    = fvec3(0.0, 0.0, 0.0);
	Physics::rbodyvec[2].velocity = fvec3(0.0, 0.0, 1.0);
	Physics::rbodyvec[2].tva.setcolor(0.2, 0.8, 0.3, 1.0);
	Physics::rbodyvec[2].tva.settype(2);
	Physics::rbodyvec[2].lva.setcolor(0.8, 0.8, 0.8, 1.0);

	Physics::rbodyvec.emplace_back("../../../resource/Cube.obj", "../../../resource/Cube.obj", fvec3(0.0, 2.0, 0.0), 2.00, 1.00);
	Physics::rbodyvec[3].rotq     = fquaternion(fvec3(.0 * 3.1415, 0.25 * 3.1415, 0.0));
	Physics::rbodyvec[3].omega    = fvec3(0.0, 0.0, 0.0);
	Physics::rbodyvec[3].velocity = fvec3(0.0, 0.0, 1.0);
	Physics::rbodyvec[3].tva.setcolor(0.8, 0.8, 0.8, 1.0);
	Physics::rbodyvec[3].tva.settype(2);
	Physics::rbodyvec[3].lva.setcolor(0.8, 0.8, 0.8, 1.0);

	Physics::rbodyvec.emplace_back("../../../resource/Bunny.obj", "../../../resource/LightBunny.obj", fvec3(0.0, 4.0, 0.0), 2.00, 1.00);
	Physics::rbodyvec[4].rotq     = fquaternion(fvec3(.0 * 3.1415, 0.25 * 3.1415, 0.0));
	Physics::rbodyvec[4].omega    = fvec3(0.0, 0.0, 0.0);
	Physics::rbodyvec[4].velocity = fvec3(0.0, 0.0, 1.0);
	Physics::rbodyvec[4].tva.setcolor(0.3, 0.3, 0.3, 1.0);
	Physics::rbodyvec[4].tva.settype(1);
	Physics::rbodyvec[4].lva.setcolor(0.8, 0.8, 0.8, 1.0);
	*/

	////////////////////////////////////

	//衝突するBunny
	/*
	uint32_t objectsize = 4;
	Physics::rbodyvec.emplace_back("../../../resource/Bunny.obj", "../../../resource/LightBunny.obj", fvec3(0.0, -11.0, 0.0), 6.00, 5.00);
	Physics::rbodyvec[0].rotq     = fquaternion(fvec3(.0 * 3.1415, 2.00 * 3.1415, 0.0));
	Physics::rbodyvec[0].omega    = fvec3(0.0, 4.0, 0.0);
	Physics::rbodyvec[0].velocity = fvec3(0.0, 4.0, 0.0);
	Physics::rbodyvec[0].tva.setcolor(0.8, 0.4, 0.3, 1.0);
	Physics::rbodyvec[0].tva.settype(2);
	Physics::rbodyvec[0].lva.setcolor(0.8, 0.8, 0.8, 1.0);

	Physics::rbodyvec.emplace_back("../../../resource/Bunny.obj", "../../../resource/LightBunny.obj", fvec3(2.0, 2.0, 0.0), 3.00, 3.00);
	Physics::rbodyvec[1].rotq     = fquaternion(fvec3(.0 * 3.1415, 0.25 * 3.1415, 0.0 * 3.1415));
	Physics::rbodyvec[1].omega    = fvec3(0.0, 0.0, 0.0);
	Physics::rbodyvec[1].velocity = fvec3(0.0, 0.0, 0.0);
	Physics::rbodyvec[1].tva.setcolor(0.6, 0.1, 0.7, 1.0);
	Physics::rbodyvec[1].tva.settype(2);
	Physics::rbodyvec[1].lva.setcolor(0.8, 0.8, 0.8, 1.0);

	Physics::rbodyvec.emplace_back("../../../resource/Bunny.obj", "../../../resource/LightBunny.obj", fvec3(-2.0, 4.0, 2.0), 3.00, 3.00);
	Physics::rbodyvec[2].rotq     = fquaternion(fvec3(.0 * 3.1415, 0.25 * 3.1415, 0.0 * 3.1415));
	Physics::rbodyvec[2].omega    = fvec3(2.0, 0.0, 2.0);
	Physics::rbodyvec[2].velocity = fvec3(2.0, 0.0, -2.0);
	Physics::rbodyvec[2].tva.setcolor(0.6, 0.1, 0.7, 1.0);
	Physics::rbodyvec[2].tva.settype(2);
	Physics::rbodyvec[2].lva.setcolor(0.8, 0.8, 0.8, 1.0);

	Physics::rbodyvec.emplace_back("../../../resource/Bunny.obj", "../../../resource/LightBunny.obj", fvec3(-2.0, 6.0, 2.0), 2.00, 3.00);
	Physics::rbodyvec[3].rotq     = fquaternion(fvec3(.0 * 3.1415, 0.25 * 3.1415, 0.0 * 3.1415));
	Physics::rbodyvec[3].omega    = fvec3(2.0, 4.0, 0.0);
	Physics::rbodyvec[3].velocity = fvec3(2.0, -4.0, 0.0);
	Physics::rbodyvec[3].tva.setcolor(0.6, 0.1, 0.7, 1.0);
	Physics::rbodyvec[3].tva.settype(2);
	Physics::rbodyvec[3].lva.setcolor(0.8, 0.8, 0.8, 1.0);
	*/

	//init

	std::vector<fvec3> Initcm(Physics::rbodyvec.size());
	std::vector<fquaternion> Initrotq(Physics::rbodyvec.size());
	std::vector<fvec3> Initvelocity(Physics::rbodyvec.size());
	std::vector<fvec3> Initomega(Physics::rbodyvec.size());

	for (uint32_t i = 0; i < Physics::rbodyvec.size(); i++) {
		Initcm[i]	= Physics::rbodyvec[i].cm;
		Initrotq[i]	= Physics::rbodyvec[i].rotq;
		Initvelocity[i] = Physics::rbodyvec[i].velocity;
		Initomega[i]	= Physics::rbodyvec[i].omega;
	}

	std::vector<Renderer3D::drawobject> shadowlist;
	std::vector<Renderer3D::drawobject> edgelist;
	std::vector<Renderer3D::drawobject> renderlist;

	shadowlist.emplace_back(Renderer3D::drawobject { floor, nullptr, nullptr, nullptr });
	for (uint32_t i = 0; i < objectsize; i++)
		shadowlist.emplace_back(Renderer3D::drawobject { Physics::rbodyvec[i].tva, nullptr, &(Physics::rbodyvec[i].rotq), &(Physics::rbodyvec[i].cm) });

	//edgelist.emplace_back(Renderer3D::drawobject { r0.tva, nullptr, &r0.rotq, &r0.cm });
	//edgelist.emplace_back(Renderer3D::drawobject { r0.lva, nullptr, &r0.rotq, &r0.cm });

	renderlist.emplace_back(Renderer3D::drawobject { floor, &simasima0, nullptr, nullptr });
	renderlist.emplace_back(Renderer3D::drawobject { cage, nullptr, nullptr, nullptr });
	for (uint32_t i = 0; i < objectsize; i++) {
		renderlist.emplace_back(Renderer3D::drawobject { Physics::rbodyvec[i].tva, &simasima1, &(Physics::rbodyvec[i].rotq), &(Physics::rbodyvec[i].cm) });
		renderlist.emplace_back(Renderer3D::drawobject { Physics::rbodyvec[i].lva, nullptr, &(Physics::rbodyvec[i].rotq), &(Physics::rbodyvec[i].cm) });
	}

	//rendering loop

	double ctime = 0.0;
	double vtime = 0.0;

	while (!Visualizer::Is_Closed()) {
		//clear buffer
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
			Physics::timestep();
			vtime += dt;
			nextframe = false;
		}

		//imgui
		{

			ImGui::Begin("RigidBodyCollision");

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

			ImGui::Text("realtime = %.1f", ctime);
			ImGui::Text("virtualtime = %.1f", vtime);

			if (ImGui::Button("reset")) {

				vtime = 0.0;

				for (uint32_t i = 0; i < Physics::rbodyvec.size(); i++) {
					Physics::rbodyvec[i].cm	      = Initcm[i];
					Physics::rbodyvec[i].rotq     = Initrotq[i];
					Physics::rbodyvec[i].velocity = Initvelocity[i];
					Physics::rbodyvec[i].omega    = Initomega[i];
				}
				Physics::StaticFrictionVec.clear();
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

		//swap buffer
		Visualizer::SwapBuffer();
		//wait event
		Visualizer::PollEvent();
	}

	return 0;
}
