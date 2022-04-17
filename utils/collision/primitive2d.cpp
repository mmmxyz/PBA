#include "utils/mathfunc/mathfunc.hpp"
#include "utils/collision/primitive2d.hpp"

float DistLinePoint(
    const fvec2& a0,
    const fvec2& a1,
    const fvec2& b)
{
	fvec2 n = a1 - a0;
	float t = -(a0 - b).dot(n) / n.sqlength();

	fvec2 v = t * n + a0;
	return (v - b).length();
}

ClosestPair DistSegmentPoint(
    const fvec2& a0,
    const fvec2& a1,
    const fvec2& b)
{
	fvec2 n = a1 - a0;
	float t = -(a0 - b).dot(n) / n.sqlength();

	if (0.0 <= t && t <= 1.0) {
		fvec2 v = t * n + a0;
		return { (v - b).length(), 0 };
	} else if (t < 0.0) {
		fvec2 v = a0;
		return { (v - b).length(), 1 };
	} else {
		fvec2 v = a1;
		return { (v - b).length(), 2 };
	}
}

ClosestPair DistTrianglePoint(
    const fvec2& a0,
    const fvec2& a1,
    const fvec2& a2,
    const fvec2& b)
{
	float D = (a1 - a0).cross(a2 - a0);

	float alpha0 = (a1 - b).cross(a2 - b) / D;
	float alpha1 = (a2 - b).cross(a0 - b) / D;
	float alpha2 = (a0 - b).cross(a1 - b) / D;

	if (0.0 <= alpha0 && alpha0 <= 1.0 && 0.0 <= alpha1 && alpha1 <= 1.0 && 0.0 <= alpha2 && alpha2 <= 1.0) {
		auto [dist01, status01] = DistSegmentPoint(a0, a1, b);
		auto [dist12, status12] = DistSegmentPoint(a1, a2, b);
		auto [dist20, status20] = DistSegmentPoint(a2, a0, b);

		if (dist01 < dist12 && dist01 < dist20) {
			return { dist01, 0 };
		} else if (dist12 < dist20 && dist12 < dist01) {
			return { dist12, 1 };
		} else {
			return { dist20, 2 };
		}

	} else {
		return { 0.0, -1 };
	}
}

ClosestPair DistTriangleTriangle(
    const fvec2& a0,
    const fvec2& a1,
    const fvec2& a2,
    const fvec2& b0,
    const fvec2& b1,
    const fvec2& b2)
{
	//(a1 - a0)x(a2 - a0) > 0

	//形式的には，sat+最も遠い頂点の近い頂点(min max)を戻す

	//todo 内積は事前に計算する

	fvec2 normallist[6];

	normallist[0] = (a1 - a0).rot();
	normallist[1] = (a2 - a1).rot();
	normallist[2] = (a0 - a2).rot();

	normallist[3] = (b1 - b0).rot();
	normallist[4] = (b2 - b1).rot();
	normallist[5] = (b0 - b2).rot();

	float maxa[6];
	float maxb[6];
	float mina[6];
	float minb[6];

	float dist  = 0.0;
	int32_t Ind = -1;

	for (uint32_t i = 0; i < 3; i++) {
		mina[i] = std::min(std::min(a0.dot(normallist[i]), a1.dot(normallist[i])), a2.dot(normallist[i]));
		minb[i] = std::min(std::min(b0.dot(normallist[i]), b1.dot(normallist[i])), b2.dot(normallist[i]));

		maxa[i] = std::max(std::max(a0.dot(normallist[i]), a1.dot(normallist[i])), a2.dot(normallist[i]));
		maxb[i] = std::max(std::max(b0.dot(normallist[i]), b1.dot(normallist[i])), b2.dot(normallist[i]));

		if (maxa[i] < minb[i])
			return { 0.0, -1 };
		if (mina[i] > maxb[i])
			return { 0.0, -1 };

		if (Ind == -1 || maxb[i] - mina[i] < dist) {
			dist = maxb[i] - mina[i];
			Ind  = i;
		}
	}

	for (uint32_t i = 3; i < 6; i++) {
		mina[i] = std::min(std::min(a0.dot(normallist[i]), a1.dot(normallist[i])), a2.dot(normallist[i]));
		minb[i] = std::min(std::min(b0.dot(normallist[i]), b1.dot(normallist[i])), b2.dot(normallist[i]));

		maxa[i] = std::max(std::max(a0.dot(normallist[i]), a1.dot(normallist[i])), a2.dot(normallist[i]));
		maxb[i] = std::max(std::max(b0.dot(normallist[i]), b1.dot(normallist[i])), b2.dot(normallist[i]));

		if (maxa[i] < minb[i])
			return { 0.0, -1 };
		if (mina[i] > maxb[i])
			return { 0.0, -1 };

		if (Ind == -1 || maxa[i] - minb[i] < dist) {
			dist = maxa[i] - minb[i];
			Ind  = i;
		}
	}

	if (Ind < 3) {
		float dist0 = b0.dot(normallist[Ind]);
		float dist1 = b1.dot(normallist[Ind]);
		float dist2 = b2.dot(normallist[Ind]);

		if (dist0 > dist1 && dist0 > dist2)
			return { dist0 - mina[Ind], 3 * Ind + 0 };
		else if (dist1 > dist2 && dist1 > dist0)
			return { dist1 - mina[Ind], 3 * Ind + 1 };
		else
			return { dist2 - mina[Ind], 3 * Ind + 2 };
	} else {
		float dist0 = a0.dot(normallist[Ind]);
		float dist1 = a1.dot(normallist[Ind]);
		float dist2 = a2.dot(normallist[Ind]);

		if (dist0 > dist1 && dist0 > dist2)
			return { dist0 - minb[Ind], 3 * Ind + 0 };
		else if (dist1 > dist2 && dist1 > dist0)
			return { dist1 - minb[Ind], 3 * Ind + 1 };
		else
			return { dist2 - minb[Ind], 3 * Ind + 2 };
	}
}

bool Is_CollideTriangleTriangle(
    const fvec2& a0,
    const fvec2& a1,
    const fvec2& a2,
    const fvec2& b0,
    const fvec2& b1,
    const fvec2& b2)
{
	//sat

	fvec2 normallist[6];

	normallist[0] = (a1 - a0).rot();
	normallist[1] = (a2 - a1).rot();
	normallist[2] = (a0 - a2).rot();
	normallist[3] = (b1 - b0).rot();
	normallist[4] = (b2 - b1).rot();
	normallist[5] = (b0 - b2).rot();

	for (uint32_t i = 0; i < 6; i++) {
		float mina = std::min(std::min(a0.dot(normallist[i]), a1.dot(normallist[i])), a2.dot(normallist[i]));
		float minb = std::min(std::min(b0.dot(normallist[i]), b1.dot(normallist[i])), b2.dot(normallist[i]));

		float maxa = std::max(std::max(a0.dot(normallist[i]), a1.dot(normallist[i])), a2.dot(normallist[i]));
		float maxb = std::max(std::max(b0.dot(normallist[i]), b1.dot(normallist[i])), b2.dot(normallist[i]));

		if (maxa < minb)
			return false;
		if (mina > maxb)
			return false;
	}

	return true;
}
