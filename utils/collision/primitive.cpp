#include <cmath>
//#include <iostream>

#include "utils/mathfunc/mathfunc.hpp"
#include "utils/collision/primitive.hpp"

static float epsilon = 0.000001;

ClosestDV DistLinePoint(const fvec3& a0, const fvec3& a1, const fvec3& b, float* const pt)
{
	fvec3 n = a1 - a0;
	float t = -(a0 - b).dot(n) / n.dot(n);

	if (pt != nullptr)
		*pt = t;

	fvec3 v	   = t * n + a0;
	float dist = (v - b).length();
	return { dist, v };
}

ClosestDV DistSegmentPoint(const fvec3& a0, const fvec3& a1, const fvec3& b)
{
	fvec3 n = a1 - a0;
	float t = -(a0 - b).dot(n) / n.dot(n);

	if (0 <= t && t <= 1.0) {
		fvec3 v	   = t * n + a0;
		float dist = (v - b).length();
		return { dist, v };
	} else if (t < 0.0) {
		fvec3 v	   = a0;
		float dist = (v - b).length();
		return { dist, v };
	} else {
		fvec3 v	   = a1;
		float dist = (v - b).length();
		return { dist, v };
	}
}

ClosestDVV DistLineLine(const fvec3& a0, const fvec3& a1, const fvec3& b0, const fvec3& b1, float* const pt, float* const ps)
{

	fvec3 na = a1 - a0;
	fvec3 nb = b1 - b0;

	fvec3 a, b;
	float dist;

	if (na.sqlength() * nb.sqlength() - (na.dot(nb)) * (na.dot(nb)) < epsilon) {
		float s = (a0.dot(na) - b0.dot(na)) / na.dot(nb);
		a	= a0;
		b	= s * (b1 - b0) + a0;

		if (pt != nullptr && ps != nullptr) {
			*pt == 0.0;
			*ps = s;
		}

		dist = (a - b).length();
		return { dist, a, b };
	}

	fmat2 hoge({ na.sqlength(), -na.dot(nb), na.dot(nb), -nb.sqlength() });
	hoge	 = hoge.inverse();
	fvec2 ts = hoge * fvec2(-a0.dot(na) + b0.dot(na), -a0.dot(nb) + b0.dot(nb));

	a = ts.x * (a1 - a0) + a0;
	b = ts.y * (b1 - b0) + b0;

	if (pt != nullptr && ps != nullptr) {
		*pt = ts.x;
		*ps = ts.y;
	}

	dist = (a - b).length();

	return { dist, a, b };
}

ClosestDVV DistLineSegment(const fvec3& a0, const fvec3& a1, const fvec3& b0, const fvec3& b1, float* const pt, float* const ps)
{
	fvec3 na = a1 - a0;
	fvec3 nb = b1 - b0;

	if (na.sqlength() * nb.sqlength() - (na.dot(nb)) * (na.dot(nb)) < epsilon) {
		auto [dist0, v0] = DistLinePoint(a0, a1, b0);
		return { dist0, v0, b0 };
	}

	fmat2 hoge({ na.sqlength(), -na.dot(nb), na.dot(nb), -nb.sqlength() });
	hoge	 = hoge.inverse();
	fvec2 ts = hoge * fvec2(-a0.dot(na) + b0.dot(na), -a0.dot(nb) + b0.dot(nb));

	fvec3 a = ts.x * (a1 - a0) + a0;
	fvec3 b = ts.y * (b1 - b0) + b0;

	if (0.0 <= ts.y && ts.y <= 1.0) {
		if (pt != nullptr && ps != nullptr) {
			*pt = ts.x;
			*ps = ts.y;
		}

		return { (a - b).length(), a, b };
	}

	if (ts.y <= 0.0) {
		if (pt != nullptr && ps != nullptr)
			*ps = 0.0;
		b = b0;
	} else {
		if (pt != nullptr && ps != nullptr)
			*ps = 1.0;
		b = b1;
	}

	auto [dist, va] = DistLinePoint(a0, a1, b, pt);
	return { dist, va, b };
}

ClosestDVV DistSegmentSegment(const fvec3& a0, const fvec3& a1, const fvec3& b0, const fvec3& b1)
{

	fvec3 na = a1 - a0;
	fvec3 nb = b1 - b0;

	fvec3 a, b;
	float dist;

	if (na.sqlength() * nb.sqlength() - (na.dot(nb)) * (na.dot(nb)) < epsilon) {

		auto [dist0, v0] = DistSegmentPoint(a0, a1, b0);
		auto [dist1, v1] = DistSegmentPoint(a0, a1, b1);
		auto [dist2, v2] = DistSegmentPoint(b0, b1, a0);
		auto [dist3, v3] = DistSegmentPoint(b0, b1, a1);

		if (dist0 <= dist1 + epsilon && dist0 <= dist2 + epsilon && dist0 <= dist3 + epsilon)
			return { dist0, v0, b0 };
		else if (dist1 <= dist0 + epsilon && dist1 <= dist2 + epsilon && dist1 <= dist3 + epsilon)
			return { dist1, v1, b1 };
		else if (dist2 <= dist0 + epsilon && dist2 <= dist1 + epsilon && dist2 <= dist3 + epsilon)
			return { dist2, a0, v2 };
		else
			return { dist3, a1, v3 };
	}

	fmat2 hoge({ na.sqlength(), -na.dot(nb), na.dot(nb), -nb.sqlength() });
	hoge	 = hoge.inverse();
	fvec2 ts = hoge * fvec2(-a0.dot(na) + b0.dot(na), -a0.dot(nb) + b0.dot(nb));
	//std::cout << ts << std::endl;

	if (0.0 <= ts.x && ts.x <= 1.0 && 0.0 <= ts.y && ts.y <= 1.0) {
		a = ts.x * (a1 - a0) + a0;
		b = ts.y * (b1 - b0) + b0;
	} else if (0.0 <= ts.x && ts.x <= 1.0) {
		if (ts.y < 0.0) {
			auto [dist, v] = DistSegmentPoint(a0, a1, b0);
			return { dist, v, b0 };
		} else {
			auto [dist, v] = DistSegmentPoint(a0, a1, b1);
			return { dist, v, b0 };
		}
	} else if (0.0 <= ts.y && ts.y <= 1.0) {
		if (ts.x < 0.0) {
			auto [dist, v] = DistSegmentPoint(b0, b1, a0);
			return { dist, a0, v };
		} else {
			auto [dist, v] = DistSegmentPoint(b0, b1, a1);
			return { dist, a1, v };
		}
	} else {
		if (ts.x < 0.0)
			a = a0;
		else
			a = a1;

		if (ts.y < 0.0)
			b = b0;
		else
			b = b1;
	}

	auto [dist1, v1] = DistSegmentPoint(a0, a1, b);
	auto [dist2, v2] = DistSegmentPoint(b0, b1, a);

	if (dist1 < dist2)
		return { dist1, v1, b };
	else
		return { dist2, a, v2 };
}

ClosestDV DistTrianglePoint(const fvec3& a0, const fvec3& a1, const fvec3& a2, const fvec3& b)
{
	fvec3 normal = ((a1 - a0).cross(a2 - a0)).normalize();
	float D	     = ((a1 - a0).cross(a2 - a0)).length();

	float alpha0 = ((a1 - b).cross(a2 - b)).dot(normal);
	alpha0 /= D;
	float alpha1 = ((a2 - b).cross(a0 - b)).dot(normal);
	alpha1 /= D;
	float alpha2 = ((a0 - b).cross(a1 - b)).dot(normal);
	alpha2 /= D;

	//if (0.0 <= alpha0 && alpha0 <= 1.0 && 0.0 <= alpha1 && alpha1 <= 1.0 && 0.0 <= alpha2 && alpha2 <= 1.0) {
	if (0.0 <= alpha0 && 0.0 <= alpha1 && 0.0 <= alpha2) {
		fvec3 v	   = alpha0 * a0 + alpha1 * a1 + alpha2 * a2;
		float dist = (v - b).length();
		return { dist, v };
	} else if ((b - a1).dot(a2 - a1) >= 0.0 && (b - a2).dot(a1 - a2) >= 0.0 && alpha0 < 0.0) {
		auto [dist, v] = DistSegmentPoint(a1, a2, b);
		return { dist, v };
	} else if ((b - a2).dot(a0 - a2) >= 0.0 && (b - a0).dot(a2 - a0) >= 0.0 && alpha1 < 0.0) {
		auto [dist, v] = DistSegmentPoint(a2, a0, b);
		return { dist, v };
	} else if ((b - a0).dot(a1 - a0) >= 0.0 && (b - a1).dot(a0 - a1) >= 0.0 && alpha2 < 0.0) {
		auto [dist, v] = DistSegmentPoint(a0, a1, b);
		return { dist, v };
	} else if ((b - a0).dot(a1 - a0) < 0.0 && (b - a0).dot(a2 - a0) < 0.0) {
		return { (a0 - b).length(), a0 };
	} else if ((b - a1).dot(a2 - a1) < 0.0 && (b - a1).dot(a0 - a1) < 0.0) {
		return { (a1 - b).length(), a1 };
	} else if ((b - a2).dot(a0 - a2) < 0.0 && (b - a2).dot(a1 - a2) < 0.0) {
		return { (a2 - b).length(), a2 };
	}

	return { 0.0, fvec3(0.0) };
}

ClosestDVV DistTriangleLine(const fvec3& a0, const fvec3& a1, const fvec3& a2, const fvec3& b0, const fvec3& b1, float* const ps)
{
	fvec3 normal = ((a1 - a0).cross(a2 - a0)).normalize();
	float D	     = ((a1 - a0).cross(a2 - a0)).length();

	if (std::abs((b1 - b0).dot(normal)) >= epsilon) {

		float s = ((a0 - b0).dot(normal) / (b1 - b0).dot(normal));
		fvec3 b = s * (b1 - b0) + b0;

		if (ps != nullptr)
			*ps = s;

		auto [dist, v] = DistTrianglePoint(a0, a1, a2, b);
		//std::cout << dist << std::endl;
		if (dist <= epsilon) {
			return { dist, v, b };
		}
	}

	float s;
	float p01, p12, p20;
	auto [dist01, vb01, va01] = DistLineSegment(b0, b1, a0, a1, &p01, &s);
	auto [dist12, vb12, va12] = DistLineSegment(b0, b1, a1, a2, &p12, &s);
	auto [dist20, vb20, va20] = DistLineSegment(b0, b1, a2, a0, &p20, &s);

	//ここではいずれかのパラメータを返す。
	if (dist01 <= dist12 + epsilon && dist01 <= dist20 + epsilon) {
		if (ps != nullptr)
			*ps = p01;
		return { dist01, va01, vb01 };
	} else if (dist12 <= dist01 + epsilon && dist12 <= dist20 + epsilon) {
		if (ps != nullptr)
			*ps = p12;
		return { dist12, va12, vb12 };
	} else {
		if (ps != nullptr)
			*ps = p20;
		return { dist20, va20, vb20 };
	}
}

ClosestDVV DistTriangleSegment(const fvec3& a0, const fvec3& a1, const fvec3& a2, const fvec3& b0, const fvec3& b1)
{
	fvec3 normal = ((a1 - a0).cross(a2 - a0)).normalize();
	float D	     = ((a1 - a0).cross(a2 - a0)).length();

	if (std::abs((b1 - b0).dot(normal)) >= epsilon) {

		float s;
		auto [dist, va, vb] = DistTriangleLine(a0, a1, a2, b0, b1, &s);

		if (s <= 0.0) {
			auto [dist, a] = DistTrianglePoint(a0, a1, a2, b0);
			return { dist, a, b0 };
		} else if (1.0 <= s) {
			auto [dist, a] = DistTrianglePoint(a0, a1, a2, b1);
			return { dist, a, b1 };
		} else {
			return { dist, va, vb };
		}

	} else {
		float s;
		auto [dist, va, vb] = DistTriangleLine(a0, a1, a2, b0, b1, &s);
		if (0.0 <= s && s <= 1.0) {
			return { dist, va, vb };
		} else {

			auto [dist0, va0] = DistTrianglePoint(a0, a1, a2, b0);
			auto [dist1, va1] = DistTrianglePoint(a0, a1, a2, b1);

			if (dist0 <= dist1)
				return { dist0, va0, b0 };
			else
				return { dist1, va1, b1 };
		}
	}
}

ClosestDVV DistTriangleTriangle(const fvec3& a0, const fvec3& a1, const fvec3& a2, const fvec3& b0, const fvec3& b1, const fvec3& b2)
{

	ClosestDVV result[6];

	result[0] = DistTriangleSegment(a0, a1, a2, b0, b1);
	result[1] = DistTriangleSegment(a0, a1, a2, b1, b2);
	result[2] = DistTriangleSegment(a0, a1, a2, b2, b0);

	result[3] = DistTriangleSegment(b0, b1, b2, a0, a1);
	result[4] = DistTriangleSegment(b0, b1, b2, a1, a2);
	result[5] = DistTriangleSegment(b0, b1, b2, a2, a0);

	for (uint32_t i = 3; i < 6; i++) {
		fvec3 hoge   = result[i].v0;
		result[i].v0 = result[i].v1;
		result[i].v1 = hoge;
	}

	float mindist	= result[0].dist;
	uint32_t minind = 0;

	for (uint32_t i = 1; i < 6; i++) {
		//std::cout << result[i].dist << std::endl;
		if (mindist >= result[i].dist) {
			mindist = result[i].dist;
			minind	= i;
		}
	}
	//std::cout << minind << std::endl;
	//std::cout << std::endl;

	return result[minind];
}
