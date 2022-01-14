#pragma once

#include <cstdint>

#include "utils/mathfunc/mathfunc.hpp"

class bezier_curve {
	const uint32_t N; //N次
	fvec3* const P;	  //制御点
    public:
	bezier_curve(const uint32_t& n, const fvec3* p)
	    : N(n)
	    , P(new fvec3[n + 1])
	{
		//if(n < 0)

		if (p != nullptr) {
			for (uint32_t i = 0; i < N + 1; i++)
				P[i] = p[i];
		}
	}
	~bezier_curve()
	{
		delete[] P;
	}

	void setP(const uint32_t& ind, const fvec3& v)
	{
		P[ind] = v;
	}
	void setP(const fvec3* p)
	{
		if (p != nullptr) {
			for (uint32_t i = 0; i < N + 1; i++)
				P[i] = p[i];
		}
	}

	fvec3 operator()(const float& t);

	bezier_curve diff()
	{
		fvec3* Q = new fvec3[N];
		for (uint32_t i = 0; i < N; i++)
			Q[i] = N * (P[i + 1] - P[i]);

		return bezier_curve(N - 1, Q);
	}
};

class b_spline {
	uint32_t p;   //p次
	uint32_t n;   //制御点の数はn+1
	fvec3* P;     //制御点
	float* knots; //knotベクトル サイズはn+p+2，clampされていると仮定

    public:
	b_spline(const uint32_t& deg, const uint32_t& size, const fvec3* cp, const float* t)
	    : p(deg)
	    , n(size)
	    , P(new fvec3[n + 1])
	    , knots(new float[n + p + 2])
	{
		if (cp != nullptr)
			for (uint32_t i = 0; i < n + 1; i++)
				P[i] = cp[i];

		if (t != nullptr)
			for (uint32_t i = 0; i < n + p + 2; i++)
				knots[i] = t[i];
		else {
			for (uint32_t i = 0; i < p + 1; i++)
				knots[i] = 0.0;
			for (uint32_t i = p + 1; i < n + 1; i++)
				knots[i] = (i - p) / float(n - p + 1);
			for (uint32_t i = n + 1; i < n + p + 2; i++)
				knots[i] = 1.0;
		}
	}
	~b_spline()
	{

		delete[] P;
		delete[] knots;
	}

	b_spline& operator=(const b_spline& bs)
	{
		if (this == &bs)
			return *this;

		delete[] P;
		delete[] knots;

		p = bs.p;
		n = bs.n;

		P     = new fvec3[n + 1];
		knots = new float[n + p + 2];

		for (uint32_t i = 0; i <= n; i++)
			P[i] = bs.P[i];

		for (uint32_t i = 0; i < n + p + 2; i++)
			knots[i] = bs.knots[i];
	}

	void setdegree(const uint32_t& deg, const float* t)
	{

		delete[] knots;
		this->p = deg;
		knots	= new float[n + p + 2];

		if (t != nullptr)
			for (uint32_t i = 0; i < n + p + 2; i++)
				knots[i] = t[i];
		else {
			for (uint32_t i = 0; i < p + 1; i++)
				knots[i] = 0.0;
			for (uint32_t i = p + 1; i < n + 1; i++)
				knots[i] = (i - p) / float(n - p + 1);
			for (uint32_t i = n + 1; i < n + p + 2; i++)
				knots[i] = 1.0;
		}
	}

	void setP(const uint32_t& ind, const fvec3& v)
	{
		P[ind] = v;
	}
	void setP(const fvec3* p)
	{
		if (p != nullptr) {
			for (uint32_t i = 0; i < n + 1; i++)
				P[i] = p[i];
		}
	}

	fvec3 operator()(const float& t);

	b_spline diff();
};

class piecewise_cubic_hermite {
	const uint32_t n; //制御点の数
	fvec3* const P;	  //制御点 n+1個
	fvec3* const V;	  //制御点での勾配

    public:
	piecewise_cubic_hermite(const uint32_t& N, const fvec3* p, const fvec3* v)
	    : n(N)
	    , P(new fvec3[n + 1])
	    , V(new fvec3[n + 1])
	{

		if (p != nullptr) {
			for (uint32_t i = 0; i < n + 1; i++)
				P[i] = p[i];
		}
		if (v != nullptr) {
			for (uint32_t i = 0; i < n + 1; i++)
				V[i] = v[i];
		}
	}

	fvec3 operator()(const float& t);

	fvec3 derive(const float& t);
};

class catmull_rom {
	const uint32_t n; //制御点の数
	fvec3* const P;	  //制御点 n+1個
	fvec3* const V;	  //制御点での勾配

    public:
	catmull_rom(const uint32_t& N, const fvec3* p)
	    : n(N)
	    , P(new fvec3[n + 1])
	    , V(new fvec3[n + 1])
	{

		if (p != nullptr) {
			for (uint32_t i = 0; i < n + 1; i++)
				P[i] = p[i];
		}

		V[0] = (P[1] - P[0]) / 2.0;
		V[n] = (P[n] - P[n - 1]) / 2.0;

		for (uint32_t i = 1; i < n; i++)
			V[i] = (P[i + 1] - P[i - 1]) / 2.0;
	}

	fvec3 operator()(const float& t);

	fvec3 derive(const float& t);
};
