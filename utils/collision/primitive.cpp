#include <cmath>
//#include <iostream>

#include "utils/mathfunc/mathfunc.hpp"
#include "utils/collision/primitive.hpp"

//TODO
//アプリケーションによって要求される情報は異なる
//衝突の有無，最短距離，修正ベクトル，最短距離を構成する頂点．．
//primitive -> basic として，その点を修正する

constexpr float epsilon = 0.0001;

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

bool Is_CollideTetraTetra(
    const fvec3& a0,
    const fvec3& a1,
    const fvec3& a2,
    const fvec3& a3,
    const fvec3& b0,
    const fvec3& b1,
    const fvec3& b2,
    const fvec3& b3)
{

	//sat

	fvec3 normallist[8 + 6 * 6];

	normallist[0] = (a1 - a0).cross(a2 - a0);
	normallist[1] = (a2 - a1).cross(a3 - a1);
	normallist[2] = (a3 - a2).cross(a0 - a2);
	normallist[3] = (a0 - a3).cross(a1 - a3);

	normallist[4] = (b1 - b0).cross(b2 - b0);
	normallist[5] = (b2 - b1).cross(b3 - b1);
	normallist[6] = (b3 - b2).cross(b0 - b2);
	normallist[7] = (b0 - b3).cross(b1 - b3);

	normallist[8 + 0 * 6 + 0] = (a1 - a0).cross(b1 - b0);
	normallist[8 + 0 * 6 + 1] = (a1 - a0).cross(b2 - b1);
	normallist[8 + 0 * 6 + 2] = (a1 - a0).cross(b3 - b2);
	normallist[8 + 0 * 6 + 3] = (a1 - a0).cross(b0 - b3);
	normallist[8 + 0 * 6 + 4] = (a1 - a0).cross(b2 - b0);
	normallist[8 + 0 * 6 + 5] = (a1 - a0).cross(b3 - b1);

	normallist[8 + 1 * 6 + 0] = (a2 - a1).cross(b1 - b0);
	normallist[8 + 1 * 6 + 1] = (a2 - a1).cross(b2 - b1);
	normallist[8 + 1 * 6 + 2] = (a2 - a1).cross(b3 - b2);
	normallist[8 + 1 * 6 + 3] = (a2 - a1).cross(b0 - b3);
	normallist[8 + 1 * 6 + 4] = (a2 - a1).cross(b2 - b0);
	normallist[8 + 1 * 6 + 5] = (a2 - a1).cross(b3 - b1);

	normallist[8 + 2 * 6 + 0] = (a3 - a2).cross(b1 - b0);
	normallist[8 + 2 * 6 + 1] = (a3 - a2).cross(b2 - b1);
	normallist[8 + 2 * 6 + 2] = (a3 - a2).cross(b3 - b2);
	normallist[8 + 2 * 6 + 3] = (a3 - a2).cross(b0 - b3);
	normallist[8 + 2 * 6 + 4] = (a3 - a2).cross(b2 - b0);
	normallist[8 + 2 * 6 + 5] = (a3 - a2).cross(b3 - b1);

	normallist[8 + 3 * 6 + 0] = (a0 - a3).cross(b1 - b0);
	normallist[8 + 3 * 6 + 1] = (a0 - a3).cross(b2 - b1);
	normallist[8 + 3 * 6 + 2] = (a0 - a3).cross(b3 - b2);
	normallist[8 + 3 * 6 + 3] = (a0 - a3).cross(b0 - b3);
	normallist[8 + 3 * 6 + 4] = (a0 - a3).cross(b2 - b0);
	normallist[8 + 3 * 6 + 5] = (a0 - a3).cross(b3 - b1);

	normallist[8 + 4 * 6 + 0] = (a2 - a0).cross(b1 - b0);
	normallist[8 + 4 * 6 + 1] = (a2 - a0).cross(b2 - b1);
	normallist[8 + 4 * 6 + 2] = (a2 - a0).cross(b3 - b2);
	normallist[8 + 4 * 6 + 3] = (a2 - a0).cross(b0 - b3);
	normallist[8 + 4 * 6 + 4] = (a2 - a0).cross(b2 - b0);
	normallist[8 + 4 * 6 + 5] = (a2 - a0).cross(b3 - b1);

	normallist[8 + 5 * 6 + 0] = (a3 - a1).cross(b1 - b0);
	normallist[8 + 5 * 6 + 1] = (a3 - a1).cross(b2 - b1);
	normallist[8 + 5 * 6 + 2] = (a3 - a1).cross(b3 - b2);
	normallist[8 + 5 * 6 + 3] = (a3 - a1).cross(b0 - b3);
	normallist[8 + 5 * 6 + 4] = (a3 - a1).cross(b2 - b0);
	normallist[8 + 5 * 6 + 5] = (a3 - a1).cross(b3 - b1);

	for (uint32_t i = 0; i < 44; i++) {

		if (normallist[i].sqlength() < 0.00001)
			continue;

		float maxa = std::max(std::max(std::max(a0.dot(normallist[i]), a1.dot(normallist[i])), a2.dot(normallist[i])), a3.dot(normallist[i]));
		float mina = std::min(std::min(std::min(a0.dot(normallist[i]), a1.dot(normallist[i])), a2.dot(normallist[i])), a3.dot(normallist[i]));
		float maxb = std::max(std::max(std::max(b0.dot(normallist[i]), b1.dot(normallist[i])), b2.dot(normallist[i])), b3.dot(normallist[i]));
		float minb = std::min(std::min(std::min(b0.dot(normallist[i]), b1.dot(normallist[i])), b2.dot(normallist[i])), b3.dot(normallist[i]));

		if (maxa - epsilon < minb || maxb - epsilon < mina)
			return false;
	}

	return true;
}

//ClosestPair3D DistTetraTetra(
//    const fvec3& a0,
//    const fvec3& a1,
//    const fvec3& a2,
//    const fvec3& a3,
//    const fvec3& b0,
//    const fvec3& b1,
//    const fvec3& b2,
//    const fvec3& b3)
//{
//	//sat + 最も遠い頂点のうち，最短のものを選択
//
//	//epaのほうが簡単そうだが，面と面が向かい合うときに面倒
//
//	fvec3 normallist[8 + 6 * 6];
//
//	normallist[0] = (a1 - a0).cross(a2 - a0);
//	normallist[1] = (a2 - a1).cross(a3 - a1);
//	normallist[2] = (a3 - a2).cross(a0 - a2);
//	normallist[3] = (a0 - a3).cross(a1 - a3);
//
//	normallist[4] = (b1 - b0).cross(b2 - b0);
//	normallist[5] = (b2 - b1).cross(b3 - b1);
//	normallist[6] = (b3 - b2).cross(b0 - b2);
//	normallist[7] = (b0 - b3).cross(b1 - b3);
//
//	normallist[8 + 0 * 6 + 0] = (a1 - a0).cross(b1 - b0);
//	normallist[8 + 0 * 6 + 1] = (a1 - a0).cross(b2 - b1);
//	normallist[8 + 0 * 6 + 2] = (a1 - a0).cross(b3 - b2);
//	normallist[8 + 0 * 6 + 3] = (a1 - a0).cross(b0 - b3);
//	normallist[8 + 0 * 6 + 4] = (a1 - a0).cross(b2 - b0);
//	normallist[8 + 0 * 6 + 5] = (a1 - a0).cross(b3 - b1);
//
//	normallist[8 + 1 * 6 + 0] = (a2 - a1).cross(b1 - b0);
//	normallist[8 + 1 * 6 + 1] = (a2 - a1).cross(b2 - b1);
//	normallist[8 + 1 * 6 + 2] = (a2 - a1).cross(b3 - b2);
//	normallist[8 + 1 * 6 + 3] = (a2 - a1).cross(b0 - b3);
//	normallist[8 + 1 * 6 + 4] = (a2 - a1).cross(b2 - b0);
//	normallist[8 + 1 * 6 + 5] = (a2 - a1).cross(b3 - b1);
//
//	normallist[8 + 2 * 6 + 0] = (a3 - a2).cross(b1 - b0);
//	normallist[8 + 2 * 6 + 1] = (a3 - a2).cross(b2 - b1);
//	normallist[8 + 2 * 6 + 2] = (a3 - a2).cross(b3 - b2);
//	normallist[8 + 2 * 6 + 3] = (a3 - a2).cross(b0 - b3);
//	normallist[8 + 2 * 6 + 4] = (a3 - a2).cross(b2 - b0);
//	normallist[8 + 2 * 6 + 5] = (a3 - a2).cross(b3 - b1);
//
//	normallist[8 + 3 * 6 + 0] = (a0 - a3).cross(b1 - b0);
//	normallist[8 + 3 * 6 + 1] = (a0 - a3).cross(b2 - b1);
//	normallist[8 + 3 * 6 + 2] = (a0 - a3).cross(b3 - b2);
//	normallist[8 + 3 * 6 + 3] = (a0 - a3).cross(b0 - b3);
//	normallist[8 + 3 * 6 + 4] = (a0 - a3).cross(b2 - b0);
//	normallist[8 + 3 * 6 + 5] = (a0 - a3).cross(b3 - b1);
//
//	normallist[8 + 4 * 6 + 0] = (a2 - a0).cross(b1 - b0);
//	normallist[8 + 4 * 6 + 1] = (a2 - a0).cross(b2 - b1);
//	normallist[8 + 4 * 6 + 2] = (a2 - a0).cross(b3 - b2);
//	normallist[8 + 4 * 6 + 3] = (a2 - a0).cross(b0 - b3);
//	normallist[8 + 4 * 6 + 4] = (a2 - a0).cross(b2 - b0);
//	normallist[8 + 4 * 6 + 5] = (a2 - a0).cross(b3 - b1);
//
//	normallist[8 + 5 * 6 + 0] = (a3 - a1).cross(b1 - b0);
//	normallist[8 + 5 * 6 + 1] = (a3 - a1).cross(b2 - b1);
//	normallist[8 + 5 * 6 + 2] = (a3 - a1).cross(b3 - b2);
//	normallist[8 + 5 * 6 + 3] = (a3 - a1).cross(b0 - b3);
//	normallist[8 + 5 * 6 + 4] = (a3 - a1).cross(b2 - b0);
//	normallist[8 + 5 * 6 + 5] = (a3 - a1).cross(b3 - b1);
//
//	//
//
//	float dist     = 99999.0;
//	int32_t status = -1;
//	// 0 - 63 face v.s. vertex
//	//64 -135 (0-71)
//	//
//
//	//min,maxを2回計算している todo
//	for (uint32_t i = 0; i < 4; i++) {
//		if (normallist[i].sqlength() < 0.00001)
//			continue;
//		normallist[i] = normallist[i].normalize();
//		//aのface v.s. bのvertex
//
//		float maxa = std::max(std::max(std::max(a0.dot(normallist[i]), a1.dot(normallist[i])), a2.dot(normallist[i])), a3.dot(normallist[i]));
//		float mina = std::min(std::min(std::min(a0.dot(normallist[i]), a1.dot(normallist[i])), a2.dot(normallist[i])), a3.dot(normallist[i]));
//		float maxb = std::max(std::max(std::max(b0.dot(normallist[i]), b1.dot(normallist[i])), b2.dot(normallist[i])), b3.dot(normallist[i]));
//		float minb = std::min(std::min(std::min(b0.dot(normallist[i]), b1.dot(normallist[i])), b2.dot(normallist[i])), b3.dot(normallist[i]));
//
//		if (maxa - epsilon < minb || maxb - epsilon < mina)
//			return { 0.0, -1 };
//
//		if (maxa - minb < maxb - mina) {
//			//normalは外を向いている
//
//			if (maxa - minb < dist) {
//
//				dist = maxa - minb;
//
//				if (
//				    b0.dot(normallist[i]) < b1.dot(normallist[i]) + epsilon && b0.dot(normallist[i]) < b2.dot(normallist[i]) + epsilon && b0.dot(normallist[i]) < b3.dot(normallist[i]) + epsilon)
//					status = 8 * i + 0;
//				else if (
//				    b1.dot(normallist[i]) < b2.dot(normallist[i]) + epsilon && b1.dot(normallist[i]) < b3.dot(normallist[i]) + epsilon && b1.dot(normallist[i]) < b0.dot(normallist[i]) + epsilon)
//					status = 8 * i + 1;
//				else if (
//				    b2.dot(normallist[i]) < b3.dot(normallist[i]) + epsilon && b2.dot(normallist[i]) < b0.dot(normallist[i]) + epsilon && b2.dot(normallist[i]) < b1.dot(normallist[i]) + epsilon)
//					status = 8 * i + 2;
//				else
//					status = 8 * i + 3;
//			}
//
//		} else {
//			//normalを内を向いている
//
//			if (maxb - mina < dist) {
//
//				dist = maxb - mina;
//
//				if (
//				    b0.dot(normallist[i]) > b1.dot(normallist[i]) - epsilon && b0.dot(normallist[i]) > b2.dot(normallist[i]) - epsilon && b0.dot(normallist[i]) > b3.dot(normallist[i]) - epsilon)
//					status = 8 * i + 4;
//				else if (
//				    b1.dot(normallist[i]) > b2.dot(normallist[i]) - epsilon && b1.dot(normallist[i]) > b3.dot(normallist[i]) - epsilon && b1.dot(normallist[i]) > b0.dot(normallist[i]) - epsilon)
//					status = 8 * i + 5;
//				else if (
//				    b2.dot(normallist[i]) > b3.dot(normallist[i]) - epsilon && b2.dot(normallist[i]) > b0.dot(normallist[i]) - epsilon && b2.dot(normallist[i]) > b1.dot(normallist[i]) - epsilon)
//					status = 8 * i + 6;
//				else
//					status = 8 * i + 7;
//			}
//		}
//	}
//
//	for (uint32_t i = 4; i < 8; i++) {
//		if (normallist[i].sqlength() < 0.00001)
//			continue;
//		normallist[i] = normallist[i].normalize();
//
//		//bのface v.s. aのvertex
//
//		float maxa = std::max(std::max(std::max(a0.dot(normallist[i]), a1.dot(normallist[i])), a2.dot(normallist[i])), a3.dot(normallist[i]));
//		float mina = std::min(std::min(std::min(a0.dot(normallist[i]), a1.dot(normallist[i])), a2.dot(normallist[i])), a3.dot(normallist[i]));
//		float maxb = std::max(std::max(std::max(b0.dot(normallist[i]), b1.dot(normallist[i])), b2.dot(normallist[i])), b3.dot(normallist[i]));
//		float minb = std::min(std::min(std::min(b0.dot(normallist[i]), b1.dot(normallist[i])), b2.dot(normallist[i])), b3.dot(normallist[i]));
//
//		if (maxa - epsilon < minb || maxb - epsilon < mina)
//			return { 0.0, -1 };
//
//		if (maxb - mina < maxa - minb) {
//			//normalは外を向いている
//
//			if (maxb - mina < dist) {
//
//				dist = maxb - mina;
//
//				if (
//				    a0.dot(normallist[i]) < a1.dot(normallist[i]) + epsilon && a0.dot(normallist[i]) < a2.dot(normallist[i]) + epsilon && a0.dot(normallist[i]) < a3.dot(normallist[i]) + epsilon)
//					status = 8 * i + 0;
//				else if (
//				    a1.dot(normallist[i]) < a2.dot(normallist[i]) + epsilon && a1.dot(normallist[i]) < a3.dot(normallist[i]) + epsilon && a1.dot(normallist[i]) < a0.dot(normallist[i]) + epsilon)
//					status = 8 * i + 1;
//				else if (
//				    a2.dot(normallist[i]) < a3.dot(normallist[i]) + epsilon && a2.dot(normallist[i]) < a0.dot(normallist[i]) + epsilon && a2.dot(normallist[i]) < a1.dot(normallist[i]) + epsilon)
//					status = 8 * i + 2;
//				else
//					status = 8 * i + 3;
//			}
//
//		} else {
//			//normalを内を向いている
//
//			if (maxa - minb < dist) {
//
//				dist = maxa - minb;
//
//				if (
//				    a0.dot(normallist[i]) > a1.dot(normallist[i]) - epsilon && a0.dot(normallist[i]) > a2.dot(normallist[i]) - epsilon && a0.dot(normallist[i]) > a3.dot(normallist[i]) - epsilon)
//					status = 8 * i + 0;
//				else if (
//				    a1.dot(normallist[i]) > a2.dot(normallist[i]) - epsilon && a1.dot(normallist[i]) > a3.dot(normallist[i]) - epsilon && a1.dot(normallist[i]) > a0.dot(normallist[i]) - epsilon)
//					status = 8 * i + 1;
//				else if (
//				    a2.dot(normallist[i]) > a3.dot(normallist[i]) - epsilon && a2.dot(normallist[i]) > a0.dot(normallist[i]) - epsilon && a2.dot(normallist[i]) > a1.dot(normallist[i]) - epsilon)
//					status = 8 * i + 2;
//				else
//					status = 8 * i + 3;
//			}
//		}
//	}
//
//	//同じ法線を生成する組がある場合に，適切な最短の辺の組を計算できない
//	//iに応じて辺を構成する頂点のみ選択する
//	for (uint32_t i = 8; i < 44; i++) {
//		if (normallist[i].sqlength() < 0.00001)
//			continue;
//
//		normallist[i] = normallist[i].normalize();
//
//		//edge v.s. edge
//		//サーチする際に浸透が最も小さい組がedge v.s. edgeの場合である．
//		//gjkでミンコフスキー差の外周に現れる頂点がどちらの組で構成されるかの違いである
//		//normallist[i]を基準にしてa,bの辺のどちらが上にあるかを確認する
//
//		uint32_t Inda0, Inda1, Indb0, Indb1;
//
//		if ((i - 8) / 6 == 0) {
//			Inda0 = 0;
//			Inda1 = 1;
//		} else if ((i - 8) / 6 == 1) {
//			Inda0 = 1;
//			Inda1 = 2;
//		} else if ((i - 8) / 6 == 2) {
//			Inda0 = 2;
//			Inda1 = 3;
//		} else if ((i - 8) / 6 == 3) {
//			Inda0 = 3;
//			Inda1 = 0;
//		} else if ((i - 8) / 6 == 4) {
//			Inda0 = 0;
//			Inda1 = 2;
//		} else if ((i - 8) / 6 == 5) {
//			Inda0 = 1;
//			Inda1 = 3;
//		}
//
//		if ((i - 8) % 6 == 0) {
//			Indb0 = 0;
//			Indb1 = 1;
//		} else if ((i - 8) % 6 == 1) {
//			Indb0 = 1;
//			Indb1 = 2;
//		} else if ((i - 8) % 6 == 2) {
//			Indb0 = 2;
//			Indb1 = 3;
//		} else if ((i - 8) % 6 == 3) {
//			Indb0 = 3;
//			Indb1 = 0;
//		} else if ((i - 8) % 6 == 4) {
//			Indb0 = 0;
//			Indb1 = 2;
//		} else if ((i - 8) % 6 == 5) {
//			Indb0 = 1;
//			Indb1 = 3;
//		}
//
//		float maxa;
//		float mina;
//		float maxb;
//		float minb;
//
//		uint32_t Indmina, Indmaxa, Indminb, Indmaxb;
//
//		if (a0.dot(normallist[i]) < a1.dot(normallist[i]) + epsilon && a0.dot(normallist[i]) < a2.dot(normallist[i]) + epsilon && a0.dot(normallist[i]) < a3.dot(normallist[i]) + epsilon) {
//			Indmina = 0;
//			mina	= a0.dot(normallist[i]);
//		} else if (a1.dot(normallist[i]) < a2.dot(normallist[i]) + epsilon && a1.dot(normallist[i]) < a3.dot(normallist[i]) + epsilon && a1.dot(normallist[i]) < a0.dot(normallist[i]) + epsilon) {
//			Indmina = 1;
//			mina	= a1.dot(normallist[i]);
//		} else if (a2.dot(normallist[i]) < a3.dot(normallist[i]) + epsilon && a2.dot(normallist[i]) < a0.dot(normallist[i]) + epsilon && a2.dot(normallist[i]) < a1.dot(normallist[i]) + epsilon) {
//			Indmina = 2;
//			mina	= a2.dot(normallist[i]);
//		} else {
//			Indmina = 3;
//			mina	= a3.dot(normallist[i]);
//		}
//
//		if (a0.dot(normallist[i]) > a1.dot(normallist[i]) - epsilon && a0.dot(normallist[i]) > a2.dot(normallist[i]) - epsilon && a0.dot(normallist[i]) > a3.dot(normallist[i]) - epsilon) {
//			Indmaxa = 0;
//			maxa	= a0.dot(normallist[i]);
//		} else if (a1.dot(normallist[i]) > a2.dot(normallist[i]) - epsilon && a1.dot(normallist[i]) > a3.dot(normallist[i]) - epsilon && a1.dot(normallist[i]) > a0.dot(normallist[i]) - epsilon) {
//			Indmaxa = 1;
//			maxa	= a1.dot(normallist[i]);
//		} else if (a2.dot(normallist[i]) > a3.dot(normallist[i]) - epsilon && a2.dot(normallist[i]) > a0.dot(normallist[i]) - epsilon && a2.dot(normallist[i]) > a1.dot(normallist[i]) - epsilon) {
//			Indmaxa = 2;
//			maxa	= a2.dot(normallist[i]);
//		} else {
//			Indmaxa = 3;
//			maxa	= a3.dot(normallist[i]);
//		}
//
//		if (b0.dot(normallist[i]) < b1.dot(normallist[i]) + epsilon && b0.dot(normallist[i]) < b2.dot(normallist[i]) + epsilon && b0.dot(normallist[i]) < b3.dot(normallist[i]) + epsilon) {
//			Indminb = 0;
//			minb	= b0.dot(normallist[i]);
//		} else if (b1.dot(normallist[i]) < b2.dot(normallist[i]) + epsilon && b1.dot(normallist[i]) < b3.dot(normallist[i]) + epsilon && b1.dot(normallist[i]) < b0.dot(normallist[i]) + epsilon) {
//			Indminb = 1;
//			minb	= b1.dot(normallist[i]);
//		} else if (b2.dot(normallist[i]) < b3.dot(normallist[i]) + epsilon && b2.dot(normallist[i]) < b0.dot(normallist[i]) + epsilon && b2.dot(normallist[i]) < b1.dot(normallist[i]) + epsilon) {
//			Indminb = 2;
//			minb	= b2.dot(normallist[i]);
//		} else {
//			Indminb = 3;
//			minb	= b3.dot(normallist[i]);
//		}
//
//		if (b0.dot(normallist[i]) > b1.dot(normallist[i]) - epsilon && b0.dot(normallist[i]) > b2.dot(normallist[i]) - epsilon && b0.dot(normallist[i]) > b3.dot(normallist[i]) - epsilon) {
//			Indmaxb = 0;
//			maxb	= b0.dot(normallist[i]);
//		} else if (b1.dot(normallist[i]) > b2.dot(normallist[i]) - epsilon && b1.dot(normallist[i]) > b3.dot(normallist[i]) - epsilon && b1.dot(normallist[i]) > b0.dot(normallist[i]) - epsilon) {
//			Indmaxb = 1;
//			maxb	= b1.dot(normallist[i]);
//		} else if (b2.dot(normallist[i]) > b3.dot(normallist[i]) - epsilon && b2.dot(normallist[i]) > b0.dot(normallist[i]) - epsilon && b2.dot(normallist[i]) > b1.dot(normallist[i]) - epsilon) {
//			Indmaxb = 2;
//			maxb	= b2.dot(normallist[i]);
//		} else {
//			Indmaxb = 3;
//			maxb	= b3.dot(normallist[i]);
//		}
//
//		if (maxa - epsilon < minb || maxb - epsilon < mina)
//			return { 0.0, -1 };
//
//		if (maxa - minb < maxb - mina) {
//			//bのほうが高い
//			if (maxa - minb < dist) {
//				if ((Inda0 == Indmaxa || Inda1 == Indmaxa) && (Indb0 == Indminb || Indb1 == Indminb)) {
//					dist   = maxa - minb;
//					status = 64 + 2 * (i - 8) + 0;
//				}
//			}
//		} else {
//			//aのほうが高い
//			if (maxb - mina < dist) {
//				if ((Inda0 == Indmina || Inda1 == Indmina) && (Indb0 == Indmaxb || Indb1 == Indmaxb)) {
//					dist   = maxb - mina;
//					status = 64 + 2 * (i - 8) + 1;
//				}
//			}
//		}
//	}
//
//	return { dist, status };
//}

ClosestPair3D DistTetraTetra(
    const fvec3& a0,
    const fvec3& a1,
    const fvec3& a2,
    const fvec3& a3,
    const fvec3& b0,
    const fvec3& b1,
    const fvec3& b2,
    const fvec3& b3)
{
	//sat + 最も遠い頂点のうち，最短のものを選択

	//epaのほうが簡単そうだが，面と面が向かい合うときに面倒

	fvec3 normallist[8 + 6 * 6];

	normallist[0] = (a1 - a0).cross(a2 - a0);
	normallist[1] = (a2 - a1).cross(a3 - a1);
	normallist[2] = (a3 - a2).cross(a0 - a2);
	normallist[3] = (a0 - a3).cross(a1 - a3);

	normallist[4] = (b1 - b0).cross(b2 - b0);
	normallist[5] = (b2 - b1).cross(b3 - b1);
	normallist[6] = (b3 - b2).cross(b0 - b2);
	normallist[7] = (b0 - b3).cross(b1 - b3);

	normallist[8 + 0 * 6 + 0] = (a1 - a0).cross(b1 - b0);
	normallist[8 + 0 * 6 + 1] = (a1 - a0).cross(b2 - b1);
	normallist[8 + 0 * 6 + 2] = (a1 - a0).cross(b3 - b2);
	normallist[8 + 0 * 6 + 3] = (a1 - a0).cross(b0 - b3);
	normallist[8 + 0 * 6 + 4] = (a1 - a0).cross(b2 - b0);
	normallist[8 + 0 * 6 + 5] = (a1 - a0).cross(b3 - b1);

	normallist[8 + 1 * 6 + 0] = (a2 - a1).cross(b1 - b0);
	normallist[8 + 1 * 6 + 1] = (a2 - a1).cross(b2 - b1);
	normallist[8 + 1 * 6 + 2] = (a2 - a1).cross(b3 - b2);
	normallist[8 + 1 * 6 + 3] = (a2 - a1).cross(b0 - b3);
	normallist[8 + 1 * 6 + 4] = (a2 - a1).cross(b2 - b0);
	normallist[8 + 1 * 6 + 5] = (a2 - a1).cross(b3 - b1);

	normallist[8 + 2 * 6 + 0] = (a3 - a2).cross(b1 - b0);
	normallist[8 + 2 * 6 + 1] = (a3 - a2).cross(b2 - b1);
	normallist[8 + 2 * 6 + 2] = (a3 - a2).cross(b3 - b2);
	normallist[8 + 2 * 6 + 3] = (a3 - a2).cross(b0 - b3);
	normallist[8 + 2 * 6 + 4] = (a3 - a2).cross(b2 - b0);
	normallist[8 + 2 * 6 + 5] = (a3 - a2).cross(b3 - b1);

	normallist[8 + 3 * 6 + 0] = (a0 - a3).cross(b1 - b0);
	normallist[8 + 3 * 6 + 1] = (a0 - a3).cross(b2 - b1);
	normallist[8 + 3 * 6 + 2] = (a0 - a3).cross(b3 - b2);
	normallist[8 + 3 * 6 + 3] = (a0 - a3).cross(b0 - b3);
	normallist[8 + 3 * 6 + 4] = (a0 - a3).cross(b2 - b0);
	normallist[8 + 3 * 6 + 5] = (a0 - a3).cross(b3 - b1);

	normallist[8 + 4 * 6 + 0] = (a2 - a0).cross(b1 - b0);
	normallist[8 + 4 * 6 + 1] = (a2 - a0).cross(b2 - b1);
	normallist[8 + 4 * 6 + 2] = (a2 - a0).cross(b3 - b2);
	normallist[8 + 4 * 6 + 3] = (a2 - a0).cross(b0 - b3);
	normallist[8 + 4 * 6 + 4] = (a2 - a0).cross(b2 - b0);
	normallist[8 + 4 * 6 + 5] = (a2 - a0).cross(b3 - b1);

	normallist[8 + 5 * 6 + 0] = (a3 - a1).cross(b1 - b0);
	normallist[8 + 5 * 6 + 1] = (a3 - a1).cross(b2 - b1);
	normallist[8 + 5 * 6 + 2] = (a3 - a1).cross(b3 - b2);
	normallist[8 + 5 * 6 + 3] = (a3 - a1).cross(b0 - b3);
	normallist[8 + 5 * 6 + 4] = (a3 - a1).cross(b2 - b0);
	normallist[8 + 5 * 6 + 5] = (a3 - a1).cross(b3 - b1);

	//

	float dist     = 99999.0;
	int32_t status = -1;
	// 0 - 63 face v.s. vertex
	//64 -135 (0-71)
	//

	//min,maxを2回計算している todo
	for (uint32_t i = 0; i < 4; i++) {
		//aのface v.s. bのvertex
		if (normallist[i].sqlength() < 0.00001)
			continue;
		normallist[i] = normallist[i].normalize();

		uint32_t Inda0, Inda1, Inda2;
		bool is_fliped = false;

		if (i == 0) {
			Inda0 = 0;
			Inda1 = 1;
			Inda2 = 2;

			//法線を内向きにする
			if (a3.dot(normallist[i]) < a0.dot(normallist[i])) {
				normallist[i] = -1.0 * normallist[i];
				is_fliped     = true;
			}

		} else if (i == 1) {
			Inda0 = 1;
			Inda1 = 2;
			Inda2 = 3;

			if (a0.dot(normallist[i]) < a1.dot(normallist[i])) {
				normallist[i] = -1.0 * normallist[i];
				is_fliped     = true;
			}
		} else if (i == 2) {
			Inda0 = 2;
			Inda1 = 3;
			Inda2 = 0;

			if (a1.dot(normallist[i]) < a2.dot(normallist[i])) {
				normallist[i] = -1.0 * normallist[i];
				is_fliped     = true;
			}
		} else if (i == 3) {
			Inda0 = 3;
			Inda1 = 0;
			Inda2 = 1;

			if (a2.dot(normallist[i]) < a0.dot(normallist[i])) {
				normallist[i] = -1.0 * normallist[i];
				is_fliped     = true;
			}
		}

		float a0on = a0.dot(normallist[i]);
		float a1on = a1.dot(normallist[i]);
		float a2on = a2.dot(normallist[i]);
		float a3on = a3.dot(normallist[i]);
		float b0on = b0.dot(normallist[i]);
		float b1on = b1.dot(normallist[i]);
		float b2on = b2.dot(normallist[i]);
		float b3on = b3.dot(normallist[i]);

		float maxa = std::max(std::max(std::max(a0on, a1on), a2on), a3on);
		float mina = std::min(std::min(std::min(a0on, a1on), a2on), a3on);

		float maxb = std::max(std::max(std::max(b0on, b1on), b2on), b3on);
		float minb = std::min(std::min(std::min(b0on, b1on), b2on), b3on);

		if (maxb < mina + epsilon)
			return { 0.0, -1 };

		if (maxb - mina < dist) {

			dist = maxb - mina;

			if (
			    b0on > b1on - epsilon && b0on > b2on - epsilon && b0on > b3on - epsilon) {
				status = 8 * i + 0;
				if (is_fliped)
					status += 4;
			} else if (
			    b1on > b2on - epsilon && b1on > b3on - epsilon && b1on > b0on - epsilon) {
				status = 8 * i + 1;
				if (is_fliped)
					status += 4;
			} else if (
			    b2on > b3on + epsilon && b2on > b0on - epsilon && b2on > b1on - epsilon) {
				status = 8 * i + 2;
				if (is_fliped)
					status += 4;
			} else {
				status = 8 * i + 3;
				if (is_fliped)
					status += 4;
			}
		}
	}

	for (uint32_t i = 4; i < 8; i++) {
		//bのface v.s. aのvertex
		if (normallist[i].sqlength() < 0.00001)
			continue;
		normallist[i] = normallist[i].normalize();

		uint32_t Indb0, Indb1, Indb2;
		bool is_fliped = false;

		if (i == 0) {
			Indb0 = 0;
			Indb1 = 1;
			Indb2 = 2;

			//法線を内向きにする
			if (b3.dot(normallist[i]) < b0.dot(normallist[i])) {
				normallist[i] = -1.0 * normallist[i];
				is_fliped     = true;
			}

		} else if (i == 1) {
			Indb0 = 1;
			Indb1 = 2;
			Indb2 = 3;

			if (b0.dot(normallist[i]) < b1.dot(normallist[i])) {
				normallist[i] = -1.0 * normallist[i];
				is_fliped     = true;
			}
		} else if (i == 2) {
			Indb0 = 2;
			Indb1 = 3;
			Indb2 = 0;

			if (b1.dot(normallist[i]) < b2.dot(normallist[i])) {
				normallist[i] = -1.0 * normallist[i];
				is_fliped     = true;
			}
		} else if (i == 3) {
			Indb0 = 3;
			Indb1 = 0;
			Indb2 = 1;

			if (b2.dot(normallist[i]) < b0.dot(normallist[i])) {
				normallist[i] = -1.0 * normallist[i];
				is_fliped     = true;
			}
		}

		float a0on = a0.dot(normallist[i]);
		float a1on = a1.dot(normallist[i]);
		float a2on = a2.dot(normallist[i]);
		float a3on = a3.dot(normallist[i]);
		float b0on = b0.dot(normallist[i]);
		float b1on = b1.dot(normallist[i]);
		float b2on = b2.dot(normallist[i]);
		float b3on = b3.dot(normallist[i]);

		float maxa = std::max(std::max(std::max(a0on, a1on), a2on), a3on);
		float mina = std::min(std::min(std::min(a0on, a1on), a2on), a3on);

		float maxb = std::max(std::max(std::max(b0on, b1on), b2on), b3on);
		float minb = std::min(std::min(std::min(b0on, b1on), b2on), b3on);

		if (maxa < minb + epsilon)
			return { 0.0, -1 };

		if (maxa - minb < dist) {

			dist = maxa - minb;

			if (
			    a0on > a1on - epsilon && a0on > a2on - epsilon && a0on > a3on - epsilon) {
				status = 8 * i + 0;
				if (is_fliped)
					status += 4;
			} else if (
			    a1on > a2on - epsilon && a1on > a3on - epsilon && a1on > a0on - epsilon) {
				status = 8 * i + 1;
				if (is_fliped)
					status += 4;
			} else if (
			    a2on > a3on + epsilon && a2on > a0on - epsilon && a2on > a1on - epsilon) {
				status = 8 * i + 2;
				if (is_fliped)
					status += 4;
			} else {
				status = 8 * i + 3;
				if (is_fliped)
					status += 4;
			}
		}
	}

	//同じ法線を生成する組がある場合に，適切な最短の辺の組を計算できない
	//iに応じて辺を構成する頂点のみ選択する
	for (uint32_t i = 8; i < 44; i++) {
		if (normallist[i].sqlength() < 0.00001)
			continue;

		normallist[i] = normallist[i].normalize();

		//edge v.s. edge
		//サーチする際に浸透が最も小さい組がedge v.s. edgeの場合である．
		//gjkでミンコフスキー差の外周に現れる頂点がどちらの組で構成されるかの違いである
		//normallist[i]を基準にしてa,bの辺のどちらが上にあるかを確認する

		uint32_t Inda0, Inda1, Indb0, Indb1;

		if ((i - 8) / 6 == 0) {
			Inda0 = 0;
			Inda1 = 1;
		} else if ((i - 8) / 6 == 1) {
			Inda0 = 1;
			Inda1 = 2;
		} else if ((i - 8) / 6 == 2) {
			Inda0 = 2;
			Inda1 = 3;
		} else if ((i - 8) / 6 == 3) {
			Inda0 = 3;
			Inda1 = 0;
		} else if ((i - 8) / 6 == 4) {
			Inda0 = 0;
			Inda1 = 2;
		} else if ((i - 8) / 6 == 5) {
			Inda0 = 1;
			Inda1 = 3;
		}

		if ((i - 8) % 6 == 0) {
			Indb0 = 0;
			Indb1 = 1;
		} else if ((i - 8) % 6 == 1) {
			Indb0 = 1;
			Indb1 = 2;
		} else if ((i - 8) % 6 == 2) {
			Indb0 = 2;
			Indb1 = 3;
		} else if ((i - 8) % 6 == 3) {
			Indb0 = 3;
			Indb1 = 0;
		} else if ((i - 8) % 6 == 4) {
			Indb0 = 0;
			Indb1 = 2;
		} else if ((i - 8) % 6 == 5) {
			Indb0 = 1;
			Indb1 = 3;
		}

		float maxa;
		float mina;
		float maxb;
		float minb;

		uint32_t Indmina, Indmaxa, Indminb, Indmaxb;

		if (a0.dot(normallist[i]) < a1.dot(normallist[i]) + epsilon && a0.dot(normallist[i]) < a2.dot(normallist[i]) + epsilon && a0.dot(normallist[i]) < a3.dot(normallist[i]) + epsilon) {
			Indmina = 0;
			mina	= a0.dot(normallist[i]);
		} else if (a1.dot(normallist[i]) < a2.dot(normallist[i]) + epsilon && a1.dot(normallist[i]) < a3.dot(normallist[i]) + epsilon && a1.dot(normallist[i]) < a0.dot(normallist[i]) + epsilon) {
			Indmina = 1;
			mina	= a1.dot(normallist[i]);
		} else if (a2.dot(normallist[i]) < a3.dot(normallist[i]) + epsilon && a2.dot(normallist[i]) < a0.dot(normallist[i]) + epsilon && a2.dot(normallist[i]) < a1.dot(normallist[i]) + epsilon) {
			Indmina = 2;
			mina	= a2.dot(normallist[i]);
		} else {
			Indmina = 3;
			mina	= a3.dot(normallist[i]);
		}

		if (a0.dot(normallist[i]) > a1.dot(normallist[i]) - epsilon && a0.dot(normallist[i]) > a2.dot(normallist[i]) - epsilon && a0.dot(normallist[i]) > a3.dot(normallist[i]) - epsilon) {
			Indmaxa = 0;
			maxa	= a0.dot(normallist[i]);
		} else if (a1.dot(normallist[i]) > a2.dot(normallist[i]) - epsilon && a1.dot(normallist[i]) > a3.dot(normallist[i]) - epsilon && a1.dot(normallist[i]) > a0.dot(normallist[i]) - epsilon) {
			Indmaxa = 1;
			maxa	= a1.dot(normallist[i]);
		} else if (a2.dot(normallist[i]) > a3.dot(normallist[i]) - epsilon && a2.dot(normallist[i]) > a0.dot(normallist[i]) - epsilon && a2.dot(normallist[i]) > a1.dot(normallist[i]) - epsilon) {
			Indmaxa = 2;
			maxa	= a2.dot(normallist[i]);
		} else {
			Indmaxa = 3;
			maxa	= a3.dot(normallist[i]);
		}

		if (b0.dot(normallist[i]) < b1.dot(normallist[i]) + epsilon && b0.dot(normallist[i]) < b2.dot(normallist[i]) + epsilon && b0.dot(normallist[i]) < b3.dot(normallist[i]) + epsilon) {
			Indminb = 0;
			minb	= b0.dot(normallist[i]);
		} else if (b1.dot(normallist[i]) < b2.dot(normallist[i]) + epsilon && b1.dot(normallist[i]) < b3.dot(normallist[i]) + epsilon && b1.dot(normallist[i]) < b0.dot(normallist[i]) + epsilon) {
			Indminb = 1;
			minb	= b1.dot(normallist[i]);
		} else if (b2.dot(normallist[i]) < b3.dot(normallist[i]) + epsilon && b2.dot(normallist[i]) < b0.dot(normallist[i]) + epsilon && b2.dot(normallist[i]) < b1.dot(normallist[i]) + epsilon) {
			Indminb = 2;
			minb	= b2.dot(normallist[i]);
		} else {
			Indminb = 3;
			minb	= b3.dot(normallist[i]);
		}

		if (b0.dot(normallist[i]) > b1.dot(normallist[i]) - epsilon && b0.dot(normallist[i]) > b2.dot(normallist[i]) - epsilon && b0.dot(normallist[i]) > b3.dot(normallist[i]) - epsilon) {
			Indmaxb = 0;
			maxb	= b0.dot(normallist[i]);
		} else if (b1.dot(normallist[i]) > b2.dot(normallist[i]) - epsilon && b1.dot(normallist[i]) > b3.dot(normallist[i]) - epsilon && b1.dot(normallist[i]) > b0.dot(normallist[i]) - epsilon) {
			Indmaxb = 1;
			maxb	= b1.dot(normallist[i]);
		} else if (b2.dot(normallist[i]) > b3.dot(normallist[i]) - epsilon && b2.dot(normallist[i]) > b0.dot(normallist[i]) - epsilon && b2.dot(normallist[i]) > b1.dot(normallist[i]) - epsilon) {
			Indmaxb = 2;
			maxb	= b2.dot(normallist[i]);
		} else {
			Indmaxb = 3;
			maxb	= b3.dot(normallist[i]);
		}

		if (maxa - epsilon < minb || maxb - epsilon < mina)
			return { 0.0, -1 };

		if (maxa - minb < maxb - mina) {
			//bのほうが高い
			if (maxa - minb < dist) {
				if ((Inda0 == Indmaxa || Inda1 == Indmaxa) && (Indb0 == Indminb || Indb1 == Indminb)) {
					dist   = maxa - minb;
					status = 64 + 2 * (i - 8) + 0;
				}
			}
		} else {
			//aのほうが高い
			if (maxb - mina < dist) {
				if ((Inda0 == Indmina || Inda1 == Indmina) && (Indb0 == Indmaxb || Indb1 == Indmaxb)) {
					dist   = maxb - mina;
					status = 64 + 2 * (i - 8) + 1;
				}
			}
		}
	}

	return { dist, status };
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

bool Is_CollideTriTri(
    const fvec3& a0,
    const fvec3& a1,
    const fvec3& a2,
    const fvec3& b0,
    const fvec3& b1,
    const fvec3& b2)
{
	//sat

	{
		fvec3 normal = (a1 - a0).cross(a2 - a0);
		if (normal.sqlength() > 0.0000001) {
			float maxa = std::max(std::max(a0.dot(normal), a1.dot(normal)), a2.dot(normal));
			float mina = std::min(std::min(a0.dot(normal), a1.dot(normal)), a2.dot(normal));
			float maxb = std::max(std::max(b0.dot(normal), b1.dot(normal)), b2.dot(normal));
			float minb = std::min(std::min(b0.dot(normal), b1.dot(normal)), b2.dot(normal));

			if (maxa - epsilon < minb || maxb - epsilon < mina)
				return false;
		}
	}
	{
		fvec3 normal = (b1 - b0).cross(b2 - b0);
		if (normal.sqlength() > 0.0000001) {
			float maxa = std::max(std::max(a0.dot(normal), a1.dot(normal)), a2.dot(normal));
			float mina = std::min(std::min(a0.dot(normal), a1.dot(normal)), a2.dot(normal));
			float maxb = std::max(std::max(b0.dot(normal), b1.dot(normal)), b2.dot(normal));
			float minb = std::min(std::min(b0.dot(normal), b1.dot(normal)), b2.dot(normal));

			if (maxa - epsilon < minb || maxb - epsilon < mina)
				return false;
		}
	}

	fvec3 normallist[9];

	normallist[0] = (a1 - a0).cross(b1 - b0);
	normallist[1] = (a1 - a0).cross(b2 - b1);
	normallist[2] = (a1 - a0).cross(b0 - b2);

	normallist[3] = (a2 - a1).cross(b1 - b0);
	normallist[4] = (a2 - a1).cross(b2 - b1);
	normallist[5] = (a2 - a1).cross(b0 - b2);

	normallist[6] = (a0 - a2).cross(b1 - b0);
	normallist[7] = (a0 - a2).cross(b2 - b1);
	normallist[8] = (a0 - a2).cross(b0 - b2);

	for (uint32_t i = 0; i < 9; i++) {
		fvec3& normal = normallist[i];
		if (normal.sqlength() > 0.0000001) {

			float maxa = std::max(std::max(a0.dot(normal), a1.dot(normal)), a2.dot(normal));
			float mina = std::min(std::min(a0.dot(normal), a1.dot(normal)), a2.dot(normal));
			float maxb = std::max(std::max(b0.dot(normal), b1.dot(normal)), b2.dot(normal));
			float minb = std::min(std::min(b0.dot(normal), b1.dot(normal)), b2.dot(normal));

			if (maxa - epsilon < minb || maxb - epsilon < mina)
				return false;
		}
	}

	return true;
}

ClosestPair3D DistTriTri(
    const fvec3& a0,
    const fvec3& a1,
    const fvec3& a2,
    const fvec3& b0,
    const fvec3& b1,
    const fvec3& b2)
{
	float dist     = 99999999.9;
	int32_t status = -1;

	// 0-11 face v.s. vertex
	//12-29 (0-17) edge v.s. edge

	{
		fvec3 normal = (a1 - a0).cross(a2 - a0);
		if (normal.sqlength() > 0.0000001) {

			normal = normal.normalize();

			float a0on = a0.dot(normal);
			float a1on = a1.dot(normal);
			float a2on = a2.dot(normal);
			float b0on = b0.dot(normal);
			float b1on = b1.dot(normal);
			float b2on = b2.dot(normal);

			float maxa = std::max(std::max(a0on, a1on), a2on);
			float mina = std::min(std::min(a0on, a1on), a2on);
			float maxb = std::max(std::max(b0on, b1on), b2on);
			float minb = std::min(std::min(b0on, b1on), b2on);

			if (maxa - epsilon < minb || maxb - epsilon < mina)
				return { 0.0, -1 };

			if (maxa - minb < maxb - mina) {
				if (maxa - minb < dist) {
					dist = maxa - minb;
					if (b0on < b1on + epsilon && b0on < b2on + epsilon) {
						status = 0;
					} else if (b1on < b2on + epsilon && b1on < b0on + epsilon) {
						status = 1;
					} else {
						status = 2;
					}
				}
			} else {
				if (maxb - mina < dist) {
					dist = maxb - mina;
					if (b0on > b1on - epsilon && b0on > b2on - epsilon) {
						status = 3;
					} else if (b1on > b2on - epsilon && b1on > b0on - epsilon) {
						status = 4;
					} else {
						status = 5;
					}
				}
			}
		}
	}
	{
		fvec3 normal = (b1 - b0).cross(b2 - b0);
		if (normal.sqlength() > 0.0000001) {
			normal = normal.normalize();

			float a0on = a0.dot(normal);
			float a1on = a1.dot(normal);
			float a2on = a2.dot(normal);
			float b0on = b0.dot(normal);
			float b1on = b1.dot(normal);
			float b2on = b2.dot(normal);

			float maxa = std::max(std::max(a0on, a1on), a2on);
			float mina = std::min(std::min(a0on, a1on), a2on);
			float maxb = std::max(std::max(b0on, b1on), b2on);
			float minb = std::min(std::min(b0on, b1on), b2on);

			if (maxa - epsilon < minb || maxb - epsilon < mina)
				return { 0.0, -1 };

			if (maxa - minb < maxb - mina) {
				if (maxa - minb < dist) {
					dist = maxa - minb;
					if (a0on > a1on - epsilon && a0on > a2on - epsilon) {
						status = 6;
					} else if (a1on > a2on - epsilon && a1on > a0on - epsilon) {
						status = 7;
					} else {
						status = 8;
					}
				}
			} else {
				if (maxb - mina < dist) {
					dist = maxb - mina;
					if (a0on < a1on + epsilon && a0on < a2on + epsilon) {
						status = 9;
					} else if (a1on < a2on + epsilon && a1on < a0on + epsilon) {
						status = 10;
					} else {
						status = 11;
					}
				}
			}
		}
	}

	fvec3 normallist[9];

	normallist[0] = (a1 - a0).cross(b1 - b0);
	normallist[1] = (a1 - a0).cross(b2 - b1);
	normallist[2] = (a1 - a0).cross(b0 - b2);

	normallist[3] = (a2 - a1).cross(b1 - b0);
	normallist[4] = (a2 - a1).cross(b2 - b1);
	normallist[5] = (a2 - a1).cross(b0 - b2);

	normallist[6] = (a0 - a2).cross(b1 - b0);
	normallist[7] = (a0 - a2).cross(b2 - b1);
	normallist[8] = (a0 - a2).cross(b0 - b2);

	for (uint32_t i = 0; i < 9; i++) {
		fvec3 normal = normallist[i];
		if (normal.sqlength() < 0.0000001)
			continue;

		normal = normal.normalize();

		uint32_t Inda0 = i / 3;
		uint32_t Inda1 = ((i / 3) + 1) % 3;
		uint32_t Indb0 = i % 3;
		uint32_t Indb1 = ((i % 3) + 1) % 3;

		float a0on = a0.dot(normal);
		float a1on = a1.dot(normal);
		float a2on = a2.dot(normal);
		float b0on = b0.dot(normal);
		float b1on = b1.dot(normal);
		float b2on = b2.dot(normal);

		uint32_t Indmaxa, Indmina, Indmaxb, Indminb;
		float maxa, mina, maxb, minb;

		if (a0on > a1on - epsilon && a0on > a2on - epsilon) {
			Indmaxa = 0;
			maxa	= a0on;
		} else if (a1on > a2on - epsilon && a1on > a0on - epsilon) {
			Indmaxa = 1;
			maxa	= a1on;
		} else {
			Indmaxa = 2;
			maxa	= a2on;
		}

		if (b0on > b1on - epsilon && b0on > b2on - epsilon) {
			Indmaxb = 0;
			maxb	= b0on;
		} else if (b1on > b2on - epsilon && b1on > b0on - epsilon) {
			Indmaxb = 1;
			maxb	= b1on;
		} else {
			Indmaxb = 2;
			maxb	= b2on;
		}

		if (a0on < a1on + epsilon && a0on < a2on + epsilon) {
			Indmina = 0;
			mina	= a0on;
		} else if (a1on < a2on + epsilon && a1on < a0on + epsilon) {
			Indmina = 1;
			mina	= a1on;
		} else {
			Indmina = 2;
			mina	= a2on;
		}

		if (b0on < b1on + epsilon && b0on < b2on + epsilon) {
			Indminb = 0;
			minb	= b0on;
		} else if (b1on < b2on + epsilon && b1on < b0on + epsilon) {
			Indminb = 1;
			minb	= b1on;
		} else {
			Indminb = 2;
			minb	= b2on;
		}

		if (maxa - epsilon < minb || maxb - epsilon < mina)
			return { 0.0, -1 };

		if (maxa - minb < maxb - mina) {
			//aのほうが高い
			if (maxa - minb < dist) {
				if ((Inda0 == Indmaxa || Inda1 == Indmaxa) && (Indb0 == Indminb || Indb1 == Indminb)) {
					dist   = maxa - minb;
					status = 12 + 2 * i + 0;
				}
			}
		} else {
			if (maxb - mina < dist) {
				if ((Inda0 == Indmina || Inda1 == Indmina) && (Indb0 == Indmaxb || Indb1 == Indmaxb)) {
					dist   = maxb - mina;
					status = 12 + 2 * i + 1;
				}
			}
		}
	}

	return { dist, status };
}
