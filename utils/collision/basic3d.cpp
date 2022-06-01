#include <cmath>
//#include <iostream>

#include "utils/mathfunc/mathfunc.hpp"
#include "utils/collision/basic3d.hpp"

constexpr float epsilon = 0.000001;

//TODO 関数内部でより単純な問題を解く場合(Segment to SegmentからSegment to Pointを呼び出すなど)．関数を呼び出すのではなく，実装をコピペするほうが無駄がなさそう

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

		if (normallist[i].sqlength() < epsilon)
			continue;

		const float maxa = std::max(std::max(std::max(a0.dot(normallist[i]), a1.dot(normallist[i])), a2.dot(normallist[i])), a3.dot(normallist[i]));
		const float mina = std::min(std::min(std::min(a0.dot(normallist[i]), a1.dot(normallist[i])), a2.dot(normallist[i])), a3.dot(normallist[i]));
		const float maxb = std::max(std::max(std::max(b0.dot(normallist[i]), b1.dot(normallist[i])), b2.dot(normallist[i])), b3.dot(normallist[i]));
		const float minb = std::min(std::min(std::min(b0.dot(normallist[i]), b1.dot(normallist[i])), b2.dot(normallist[i])), b3.dot(normallist[i]));

		if (maxa - epsilon < minb || maxb - epsilon < mina)
			return false;
	}

	return true;
}

bool Is_CollideTriTri(
    const fvec3& a0,
    const fvec3& a1,
    const fvec3& a2,
    const fvec3& b0,
    const fvec3& b1,
    const fvec3& b2)
{
	{
		const fvec3 normal = (a1 - a0).cross(a2 - a0);
		if (normal.sqlength() > epsilon) {
			const float maxa = std::max(std::max(a0.dot(normal), a1.dot(normal)), a2.dot(normal));
			const float mina = std::min(std::min(a0.dot(normal), a1.dot(normal)), a2.dot(normal));
			const float maxb = std::max(std::max(b0.dot(normal), b1.dot(normal)), b2.dot(normal));
			const float minb = std::min(std::min(b0.dot(normal), b1.dot(normal)), b2.dot(normal));

			if (maxa - epsilon < minb || maxb - epsilon < mina)
				return false;
		}
	}
	{
		const fvec3 normal = (b1 - b0).cross(b2 - b0);
		if (normal.sqlength() > epsilon) {
			const float maxa = std::max(std::max(a0.dot(normal), a1.dot(normal)), a2.dot(normal));
			const float mina = std::min(std::min(a0.dot(normal), a1.dot(normal)), a2.dot(normal));
			const float maxb = std::max(std::max(b0.dot(normal), b1.dot(normal)), b2.dot(normal));
			const float minb = std::min(std::min(b0.dot(normal), b1.dot(normal)), b2.dot(normal));

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
		const fvec3& normal = normallist[i];
		if (normal.sqlength() > epsilon) {

			const float maxa = std::max(std::max(a0.dot(normal), a1.dot(normal)), a2.dot(normal));
			const float mina = std::min(std::min(a0.dot(normal), a1.dot(normal)), a2.dot(normal));
			const float maxb = std::max(std::max(b0.dot(normal), b1.dot(normal)), b2.dot(normal));
			const float minb = std::min(std::min(b0.dot(normal), b1.dot(normal)), b2.dot(normal));

			if (maxa - epsilon < minb || maxb - epsilon < mina)
				return false;
		}
	}

	return true;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

IntersectP4 IntersectTriangleRay(const fvec3& a0, const fvec3& a1, const fvec3& a2, const fvec3& b0, const fvec3& b1)
{
	// b0 == b1 でないことが保証される

	const fvec3 normal = ((a1 - a0).cross(a2 - a0)).normalize();
	const float D	   = ((a1 - a0).cross(a2 - a0)).length();

	if (D < epsilon) {
		//Triangleが潰れている場合は無視
		return { false };
	}

	if (std::abs((b1 - b0).dot(normal)) >= epsilon) {
		// aとbが平行でない

		const float s = ((a0 - b0).dot(normal) / (b1 - b0).dot(normal));
		const fvec3 b = s * (b1 - b0) + b0;

		const float alpha0 = ((a1 - b).cross(a2 - b)).dot(normal) / D;
		const float alpha1 = ((a2 - b).cross(a0 - b)).dot(normal) / D;
		const float alpha2 = ((a0 - b).cross(a1 - b)).dot(normal) / D;

		if (0.0 <= alpha0 && 0.0 <= alpha1 && 0.0 <= alpha2 && 0.0 <= s) {
			return { true, alpha0, alpha1, alpha2, 1.0f - s };
		} else {
			return { false };
		}
	} else {

		auto [dist01, p0] = DistLineSegment(b0, b1, a0, a1);
		auto [dist12, p1] = DistLineSegment(b0, b1, a1, a2);
		auto [dist20, p2] = DistLineSegment(b0, b1, a2, a0);

		if (dist01 <= dist12 + epsilon && dist01 <= dist20 + epsilon && dist01 < epsilon && p0[0] <= 1.0) {
			return { true, p0[1], 1.0f - p0[1], 0.0, p0[0] };
		} else if (dist12 <= dist01 + epsilon && dist12 <= dist20 + epsilon && dist12 < epsilon && p1[0] <= 1.0) {
			return { true, 0.0, p1[1], 1.0f - p1[1], p1[0] };
		} else if (dist20 <= dist01 + epsilon && dist20 <= dist12 + epsilon && dist20 < epsilon && p2[0] <= 1.0) {
			return { dist20, 1.0f - p2[1], 0.0, p2[1], p2[0] };
		} else {
			return { false };
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

DistP1 DistLinePoint(const fvec3& a0, const fvec3& a1, const fvec3& b)
{

	const fvec3 n = a1 - a0;

	if (n.sqlength() < epsilon)
		return { (b - a0).length(), 0.0 };

	const float t = (b - a0).dot(n) / n.dot(n);

	const fvec3 v = t * n + a0;

	return { (v - b).length(), 1.0f - t };
}

DistP1 DistSegmentPoint(const fvec3& a0, const fvec3& a1, const fvec3& b)
{

	const fvec3 n = a1 - a0;

	if (n.sqlength() < epsilon)
		return { (b - a0).length(), 0.0 };

	const float t = (b - a0).dot(n) / n.dot(n);

	if (0.0 <= t && t <= 1.0) {
		const fvec3 v = t * n + a0;
		return { (v - b).length(), 1.0f - t };
	} else if (t < 0.0) {
		const fvec3 v = a0;
		return { (v - b).length(), 1.0 };
	} else {
		const fvec3 v = a1;
		return { (v - b).length(), 0.0 };
	}
}

DistP2 DistLineLine(const fvec3& a0, const fvec3& a1, const fvec3& b0, const fvec3& b1)
{
	//a1==a0 || b1 == b0 でないことが保証される

	const fvec3 na = a1 - a0;
	const fvec3 nb = b1 - b0;

	// na//nb
	if (std::abs(na.sqlength() * nb.sqlength() - (na.dot(nb)) * (na.dot(nb))) < epsilon) {
		const float s = (a0.dot(na) - b0.dot(na)) / na.dot(nb);
		const fvec3 a = a0;
		const fvec3 b = s * (b1 - b0) + b0;

		return { (a - b).length(), 1.0, 1.0f - s };
	}

	const fmat2 A  = fmat2({ na.sqlength(), -na.dot(nb), na.dot(nb), -nb.sqlength() }).inverse();
	const fvec2 ts = A * fvec2(-a0.dot(na) + b0.dot(na), -a0.dot(nb) + b0.dot(nb));

	const fvec3 a = ts.x * na + a0;
	const fvec3 b = ts.y * nb + b0;

	return { (a - b).length(), 1.0f - ts.x, 1.0f - ts.y };
}

DistP2 DistLineSegment(const fvec3& a0, const fvec3& a1, const fvec3& b0, const fvec3& b1)
{

	// a1 == a0 でないことが保証される

	const fvec3 na = a1 - a0;
	const fvec3 nb = b1 - b0;

	// na//nb
	// b1 == b0
	if (na.sqlength() * nb.sqlength() - (na.dot(nb)) * (na.dot(nb)) < epsilon) {
		const auto [dist, pt] = DistLinePoint(a0, a1, b0);
		return { dist, pt, 1.0 };
	}

	const fmat2 A  = fmat2({ na.sqlength(), -na.dot(nb), na.dot(nb), -nb.sqlength() }).inverse();
	const fvec2 ts = A * fvec2(-a0.dot(na) + b0.dot(na), -a0.dot(nb) + b0.dot(nb));

	const fvec3 a = ts.x * na + a0;

	if (0.0 <= ts.y && ts.y <= 1.0) {
		const fvec3 b = ts.y * nb + b0;
		return { (a - b).length(), 1.0f - ts.x, 1.0f - ts.y };
	} else if (ts.y < 0.0) {
		const auto [dist, pt] = DistLinePoint(a0, a1, b0);
		return { dist, pt, 1.0 };
	} else {
		const auto [dist, pt] = DistLinePoint(a0, a1, b1);
		return { dist, pt, 0.0 };
	}
}

DistP2 DistSegmentSegment(const fvec3& a0, const fvec3& a1, const fvec3& b0, const fvec3& b1)
{

	const fvec3 na = a1 - a0;
	const fvec3 nb = b1 - b0;

	// na//nb
	// a1 == a0
	// b1 == b0
	if (na.sqlength() * nb.sqlength() - (na.dot(nb)) * (na.dot(nb)) < epsilon) {

		const auto [dist0, pt0] = DistSegmentPoint(a0, a1, b0);
		const auto [dist1, pt1] = DistSegmentPoint(a0, a1, b1);
		const auto [dist2, pt2] = DistSegmentPoint(b0, b1, a0);
		const auto [dist3, pt3] = DistSegmentPoint(b0, b1, a1);

		if (dist0 <= dist1 + epsilon && dist0 <= dist2 + epsilon && dist0 <= dist3 + epsilon)
			return { dist0, pt0, 1.0 };
		else if (dist1 <= dist0 + epsilon && dist1 <= dist2 + epsilon && dist1 <= dist3 + epsilon)
			return { dist1, pt1, 0.0 };
		else if (dist2 <= dist0 + epsilon && dist2 <= dist1 + epsilon && dist2 <= dist3 + epsilon)
			return { dist2, 1.0, pt2 };
		else
			return { dist3, 0.0, pt3 };
	}

	const fmat2 A  = fmat2({ na.sqlength(), -na.dot(nb), na.dot(nb), -nb.sqlength() }).inverse();
	const fvec2 ts = A * fvec2(-a0.dot(na) + b0.dot(na), -a0.dot(nb) + b0.dot(nb));

	if (0.0 <= ts.x && ts.x <= 1.0 && 0.0 <= ts.y && ts.y <= 1.0) {
		const fvec3 a = ts.x * (a1 - a0) + a0;
		const fvec3 b = ts.y * (b1 - b0) + b0;
		return { (a - b).length(), 1.0f - ts.x, 1.0f - ts.y };
	} else if (0.0 <= ts.x && ts.x <= 1.0) {
		if (ts.y < 0.0) {
			const auto [dist, pt] = DistSegmentPoint(a0, a1, b0);
			return { dist, pt, 1.0 };
		} else {
			const auto [dist, pt] = DistSegmentPoint(a0, a1, b1);
			return { dist, pt, 0.0 };
		}
	} else if (0.0 <= ts.y && ts.y <= 1.0) {
		if (ts.x < 0.0) {
			const auto [dist, pt] = DistSegmentPoint(b0, b1, a0);
			return { dist, 1.0, pt };
		} else {
			const auto [dist, pt] = DistSegmentPoint(b0, b1, a1);
			return { dist, 0.0, pt };
		}
	} else {

		if (ts.x < 0.0) {
			if (ts.y < 0.0) {
				const auto [dist1, pt1] = DistSegmentPoint(a0, a1, b0);
				const auto [dist2, pt2] = DistSegmentPoint(b0, b1, a0);

				if (dist1 < dist2)
					return { dist1, pt1, 1.0 };
				else
					return { dist2, 1.0, pt2 };
			} else {
				const auto [dist1, pt1] = DistSegmentPoint(a0, a1, b1);
				const auto [dist2, pt2] = DistSegmentPoint(b0, b1, a0);

				if (dist1 < dist2)
					return { dist1, pt1, 0.0 };
				else
					return { dist2, 1.0, pt2 };
			}

		} else {
			if (ts.y < 0.0) {
				const auto [dist1, pt1] = DistSegmentPoint(a0, a1, b0);
				const auto [dist2, pt2] = DistSegmentPoint(b0, b1, a1);

				if (dist1 < dist2)
					return { dist1, pt1, 1.0 };
				else
					return { dist2, 0.0, pt2 };

			} else {
				const auto [dist1, pt1] = DistSegmentPoint(a0, a1, b1);
				const auto [dist2, pt2] = DistSegmentPoint(b0, b1, a1);

				if (dist1 < dist2)
					return { dist1, pt1, 0.0 };
				else
					return { dist2, 0.0, pt2 };
			}
		}
	}
}

DistP3 DistTPlanePoint(const fvec3& a0, const fvec3& a1, const fvec3& a2, const fvec3& b)
{
	const fvec3 normal = ((a1 - a0).cross(a2 - a0)).normalize();
	const float D	   = ((a1 - a0).cross(a2 - a0)).length();

	if (D < epsilon) {
		//Triangleが潰れている場合は無視
		return { 0.0, 0.0, 0.0, 0.0 };
	}

	const float alpha0 = ((a1 - b).cross(a2 - b)).dot(normal) / D;
	const float alpha1 = ((a2 - b).cross(a0 - b)).dot(normal) / D;
	const float alpha2 = ((a0 - b).cross(a1 - b)).dot(normal) / D;

	const fvec3 v = alpha0 * a0 + alpha1 * a1 + alpha2 * a2;
	return { (v - b).length(), alpha0, alpha1, alpha2 };
}

DistP3 DistTrianglePoint(const fvec3& a0, const fvec3& a1, const fvec3& a2, const fvec3& b)
{
	const fvec3 normal = ((a1 - a0).cross(a2 - a0)).normalize();
	const float D	   = ((a1 - a0).cross(a2 - a0)).length();

	if (D < epsilon) {
		//Triangleが潰れている場合は無視
		return { 0.0, 0.0, 0.0, 0.0 };
	}

	const float alpha0 = ((a1 - b).cross(a2 - b)).dot(normal) / D;
	const float alpha1 = ((a2 - b).cross(a0 - b)).dot(normal) / D;
	const float alpha2 = ((a0 - b).cross(a1 - b)).dot(normal) / D;

	if (0.0 <= alpha0 && 0.0 <= alpha1 && 0.0 <= alpha2) {
		const fvec3 v = alpha0 * a0 + alpha1 * a1 + alpha2 * a2;
		return { (v - b).length(), alpha0, alpha1, alpha2 };
	} else if ((b - a1).dot(a2 - a1) > 0.0 - epsilon && (b - a2).dot(a1 - a2) > 0.0 - epsilon && alpha0 < epsilon) {
		auto [dist, pt] = DistSegmentPoint(a1, a2, b);
		return { dist, 0.0f, pt, 1.0f - pt };
	} else if ((b - a2).dot(a0 - a2) > 0.0 - epsilon && (b - a0).dot(a2 - a0) > 0.0 - epsilon && alpha1 < epsilon) {
		auto [dist, pt] = DistSegmentPoint(a2, a0, b);
		return { dist, 1.0f - pt, 0.0f, pt };
	} else if ((b - a0).dot(a1 - a0) > 0.0 - epsilon && (b - a1).dot(a0 - a1) > 0.0 - epsilon && alpha2 < epsilon) {
		auto [dist, pt] = DistSegmentPoint(a0, a1, b);
		return { dist, pt, 1.0f - pt, 0.0f };
	} else if ((b - a0).dot(a1 - a0) < 0.0 && (b - a0).dot(a2 - a0) < 0.0) {
		return { (a0 - b).length(), 1.0f, 0.0f, 0.0f };
	} else if ((b - a1).dot(a2 - a1) < 0.0 && (b - a1).dot(a0 - a1) < 0.0) {
		return { (a1 - b).length(), 0.0f, 1.0f, 0.0f };
		//} else if ((b - a2).dot(a0 - a2) < 0.0 && (b - a2).dot(a1 - a2) < 0.0) {
	} else {
		return { (a1 - b).length(), 0.0f, 0.0f, 1.0f };
	}

	//return { 0.0, fvec3(0.0) };
}

DistP4 DistTriangleLine(const fvec3& a0, const fvec3& a1, const fvec3& a2, const fvec3& b0, const fvec3& b1)
{

	// b0 == b1 でないことが保証される

	const fvec3 normal = ((a1 - a0).cross(a2 - a0)).normalize();
	const float D	   = ((a1 - a0).cross(a2 - a0)).length();

	if (D < epsilon) {
		//Triangleが潰れている場合は無視
		return { 0.0, 0.0, 0.0, 0.0 };
	}

	if (std::abs((b1 - b0).dot(normal)) >= epsilon) {
		// aとbが平行でない

		//Triangleの内部で交差する場合

		const float s = ((a0 - b0).dot(normal) / (b1 - b0).dot(normal));
		const fvec3 b = s * (b1 - b0) + b0;

		const float alpha0 = ((a1 - b).cross(a2 - b)).dot(normal) / D;
		const float alpha1 = ((a2 - b).cross(a0 - b)).dot(normal) / D;
		const float alpha2 = ((a0 - b).cross(a1 - b)).dot(normal) / D;

		if (0.0 <= alpha0 && 0.0 <= alpha1 && 0.0 <= alpha2) {
			return { 0.0, alpha0, alpha1, alpha2, 1.0f - s };
		}
	}

	auto [dist01, p0] = DistLineSegment(b0, b1, a0, a1);
	auto [dist12, p1] = DistLineSegment(b0, b1, a1, a2);
	auto [dist20, p2] = DistLineSegment(b0, b1, a2, a0);

	if (dist01 <= dist12 + epsilon && dist01 <= dist20 + epsilon) {
		return { dist01, p0[1], 1.0f - p0[1], 0.0, p0[0] };
	} else if (dist12 <= dist01 + epsilon && dist12 <= dist20 + epsilon) {
		return { dist12, 0.0, p1[1], 1.0f - p1[1], p1[0] };
	} else {
		return { dist20, 1.0f - p2[1], 0.0, p2[1], p2[0] };
	}
}

DistP4 DistTriangleSegment(const fvec3& a0, const fvec3& a1, const fvec3& a2, const fvec3& b0, const fvec3& b1)
{
	const fvec3 normal = ((a1 - a0).cross(a2 - a0)).normalize();
	const float D	   = ((a1 - a0).cross(a2 - a0)).length();

	if (D < epsilon) {
		//Triangleが潰れている場合は無視
		return { 0.0, 0.0, 0.0, 0.0 };
	}

	if (std::abs((b1 - b0).dot(normal)) >= epsilon) {
		// aとbが平行でない

		const auto [dist, p0] = DistTriangleLine(a0, a1, a2, b0, b1);

		if (1.0 <= p0[3]) {
			const auto [dist, p] = DistTrianglePoint(a0, a1, a2, b0);
			return { dist, p[0], p[1], p[2], 1.0 };
		} else if (p0[3] <= 0.0) {
			const auto [dist, p] = DistTrianglePoint(a0, a1, a2, b1);
			return { dist, p[0], p[1], p[2], 0.0 };
		} else {
			return { dist, p0[0], p0[1], p0[2], p0[3] };
		}

	} else {
		// aとbが平行
		const auto [dist, p] = DistTriangleLine(a0, a1, a2, b0, b1);

		if (0.0 <= p[3] && p[3] <= 1.0) {
			return { dist, p[0], p[1], p[2], p[3] };
		} else {

			const auto [dist0, p0] = DistTrianglePoint(a0, a1, a2, b0);
			const auto [dist1, p1] = DistTrianglePoint(a0, a1, a2, b1);

			if (dist0 <= dist1)
				return { dist0, p0[0], p0[1], p0[2], 1.0 };
			else
				return { dist1, p1[0], p1[1], p1[2], 0.0 };
		}
	}
}

DistP6 DistTriangleTriangle(const fvec3& a0, const fvec3& a1, const fvec3& a2, const fvec3& b0, const fvec3& b1, const fvec3& b2)
{

	//TODO 述べたとおり，関数内でDistTriangleSegmentを呼び出すのは筋が悪いかもしれない
	//TODO 交差する場合に最短距離の組として不正確な場所を返す
	//衝突解消のためであれば，CollisionTriangleTriangleを使うべき

	DistP4 result[6];

	result[0] = DistTriangleSegment(a0, a1, a2, b0, b1);
	result[1] = DistTriangleSegment(a0, a1, a2, b1, b2);
	result[2] = DistTriangleSegment(a0, a1, a2, b2, b0);

	result[3] = DistTriangleSegment(b0, b1, b2, a0, a1);
	result[4] = DistTriangleSegment(b0, b1, b2, a1, a2);
	result[5] = DistTriangleSegment(b0, b1, b2, a2, a0);

	float mindist	= result[0].dist;
	uint32_t minind = 0;

	for (uint32_t i = 1; i < 6; i++) {
		//std::cout << result[i].dist << std::endl;
		if (mindist >= result[i].dist) {
			mindist = result[i].dist;
			minind	= i;
		}
	}

	if (minind == 0) {
		return {
			result[0].dist,
			result[0].p[0],
			result[0].p[1],
			result[0].p[2],
			result[0].p[3],
			1.0f - result[0].p[3],
			0.0f
		};
	} else if (minind == 1) {
		return {
			result[1].dist,
			result[1].p[0],
			result[1].p[1],
			result[1].p[2],
			0.0f,
			result[1].p[3],
			1.0f - result[1].p[3]
		};
	} else if (minind == 2) {
		return {
			result[2].dist,
			result[2].p[0],
			result[2].p[1],
			result[2].p[2],
			1.0f - result[2].p[3],
			0.0f,
			result[2].p[3]
		};
	} else if (minind == 3) {
		return {
			result[3].dist,
			result[3].p[3],
			1.0f - result[3].p[3],
			0.0f,
			result[3].p[0],
			result[3].p[1],
			result[3].p[2]
		};
	} else if (minind == 4) {
		return {
			result[4].dist,
			0.0f,
			result[4].p[3],
			1.0f - result[4].p[3],
			result[4].p[0],
			result[4].p[1],
			result[4].p[2],
		};
	} else {
		return {
			result[5].dist,
			1.0f - result[5].p[3],
			0.0f,
			result[5].p[3],
			result[5].p[0],
			result[5].p[1],
			result[5].p[2]
		};
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

CollisionP6 CollisionTriangleTriangle(
    const fvec3& a0,
    const fvec3& a1,
    const fvec3& a2,
    const fvec3& b0,
    const fvec3& b1,
    const fvec3& b2)
{
	float dist     = 99999999.9;
	int32_t status = -1;
	float pt[6];

	const fvec3* const list[6] = { &a0, &a1, &a2, &b0, &b1, &b2 };

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
						const auto [hoge, p] = DistTrianglePoint(a0, a1, a2, b0);
						pt[0]		     = p[0];
						pt[1]		     = p[1];
						pt[2]		     = p[2];
						pt[3]		     = 1.0;
						pt[4]		     = 0.0;
						pt[5]		     = 0.0;
						status		     = 0;
					} else if (b1on < b2on + epsilon && b1on < b0on + epsilon) {
						const auto [hoge, p] = DistTrianglePoint(a0, a1, a2, b1);
						pt[0]		     = p[0];
						pt[1]		     = p[1];
						pt[2]		     = p[2];
						pt[3]		     = 0.0;
						pt[4]		     = 1.0;
						pt[5]		     = 0.0;
						status		     = 1;
					} else {
						const auto [hoge, p] = DistTrianglePoint(a0, a1, a2, b2);
						pt[0]		     = p[0];
						pt[1]		     = p[1];
						pt[2]		     = p[2];
						pt[3]		     = 0.0;
						pt[4]		     = 0.0;
						pt[5]		     = 1.0;
						status		     = 2;
					}
				}
			} else {
				if (maxb - mina < dist) {
					dist = maxb - mina;
					if (b0on > b1on - epsilon && b0on > b2on - epsilon) {
						const auto [hoge, p] = DistTrianglePoint(a0, a1, a2, b0);
						pt[0]		     = p[0];
						pt[1]		     = p[1];
						pt[2]		     = p[2];
						pt[3]		     = 1.0;
						pt[4]		     = 0.0;
						pt[5]		     = 0.0;
						status		     = 3;
					} else if (b1on > b2on - epsilon && b1on > b0on - epsilon) {
						const auto [hoge, p] = DistTrianglePoint(a0, a1, a2, b1);
						pt[0]		     = p[0];
						pt[1]		     = p[1];
						pt[2]		     = p[2];
						pt[3]		     = 0.0;
						pt[4]		     = 1.0;
						pt[5]		     = 0.0;
						status		     = 4;
					} else {
						const auto [hoge, p] = DistTrianglePoint(a0, a1, a2, b2);
						pt[0]		     = p[0];
						pt[1]		     = p[1];
						pt[2]		     = p[2];
						pt[3]		     = 0.0;
						pt[4]		     = 0.0;
						pt[5]		     = 1.0;
						status		     = 5;
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
						const auto [hoge, p] = DistTrianglePoint(b0, b1, b2, a0);
						pt[0]		     = 1.0;
						pt[1]		     = 0.0;
						pt[2]		     = 0.0;
						pt[3]		     = p[0];
						pt[4]		     = p[1];
						pt[5]		     = p[2];
						status		     = 6;
					} else if (a1on > a2on - epsilon && a1on > a0on - epsilon) {
						const auto [hoge, p] = DistTrianglePoint(b0, b1, b2, a1);
						pt[0]		     = 0.0;
						pt[1]		     = 1.0;
						pt[2]		     = 0.0;
						pt[3]		     = p[0];
						pt[4]		     = p[1];
						pt[5]		     = p[2];
						status		     = 7;
					} else {
						const auto [hoge, p] = DistTrianglePoint(b0, b1, b2, a2);
						pt[0]		     = 0.0;
						pt[1]		     = 0.0;
						pt[2]		     = 1.0;
						pt[3]		     = p[0];
						pt[4]		     = p[1];
						pt[5]		     = p[2];
						status		     = 8;
					}
				}
			} else {
				if (maxb - mina < dist) {
					dist = maxb - mina;
					if (a0on < a1on + epsilon && a0on < a2on + epsilon) {
						const auto [hoge, p] = DistTrianglePoint(b0, b1, b2, a0);
						pt[0]		     = 1.0;
						pt[1]		     = 0.0;
						pt[2]		     = 0.0;
						pt[3]		     = p[0];
						pt[4]		     = p[1];
						pt[5]		     = p[2];
						status		     = 9;
					} else if (a1on < a2on + epsilon && a1on < a0on + epsilon) {
						const auto [hoge, p] = DistTrianglePoint(b0, b1, b2, a1);
						pt[0]		     = 0.0;
						pt[1]		     = 1.0;
						pt[2]		     = 0.0;
						pt[3]		     = p[0];
						pt[4]		     = p[1];
						pt[5]		     = p[2];
						status		     = 10;
					} else {
						const auto [hoge, p] = DistTrianglePoint(b0, b1, b2, a2);
						pt[0]		     = 0.0;
						pt[1]		     = 0.0;
						pt[2]		     = 1.0;
						pt[3]		     = p[0];
						pt[4]		     = p[1];
						pt[5]		     = p[2];
						status		     = 11;
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

		const uint32_t Inda0 = i / 3;
		const uint32_t Inda1 = ((i / 3) + 1) % 3;
		const uint32_t Indb0 = i % 3;
		const uint32_t Indb1 = ((i % 3) + 1) % 3;

		const float a0on = a0.dot(normal);
		const float a1on = a1.dot(normal);
		const float a2on = a2.dot(normal);
		const float b0on = b0.dot(normal);
		const float b1on = b1.dot(normal);
		const float b2on = b2.dot(normal);

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
					const auto [hoge, p] = DistSegmentSegment(*list[Inda0], *list[Inda1], *list[3 + Indb0], *list[3 + Indb1]);

					pt[Inda0]		= p[0];
					pt[Inda1]		= 1.0 - p[0];
					pt[(Inda1 + 1) % 3]	= 0.0;
					pt[3 + Indb0]		= p[1];
					pt[3 + Indb1]		= 1.0 - p[1];
					pt[3 + (Indb1 + 1) % 3] = 0.0;

					dist   = maxa - minb;
					status = 12 + 2 * i + 0;
				}
			}
		} else {
			if (maxb - mina < dist) {
				if ((Inda0 == Indmina || Inda1 == Indmina) && (Indb0 == Indmaxb || Indb1 == Indmaxb)) {
					const auto [hoge, p] = DistSegmentSegment(*list[Inda0], *list[Inda1], *list[3 + Indb0], *list[3 + Indb1]);

					pt[Inda0]		= p[0];
					pt[Inda1]		= 1.0 - p[0];
					pt[(Inda1 + 1) % 3]	= 0.0;
					pt[3 + Indb0]		= p[1];
					pt[3 + Indb1]		= 1.0 - p[1];
					pt[3 + (Indb1 + 1) % 3] = 0.0;

					dist   = maxb - mina;
					status = 12 + 2 * i + 1;
				}
			}
		}
	}

	return { dist, status, pt[0], pt[1], pt[2], pt[3], pt[4], pt[5] };
}

CollisionP8 CollisionTetraTetra(const fvec3& a0, const fvec3& a1, const fvec3& a2, const fvec3& a3, const fvec3& b0, const fvec3& b1, const fvec3& b2, const fvec3& b3)
{

	constexpr uint32_t normalindlist[(8 + 6 * 6) * 4] = {
		0, 1, 2, 3,
		1, 2, 3, 0,
		3, 2, 0, 1,
		3, 0, 1, 2,

		0, 1, 2, 3,
		1, 2, 3, 0,
		3, 2, 0, 1,
		3, 0, 1, 2,

		1, 0, 1, 0,
		1, 0, 2, 1,
		1, 0, 3, 2,
		1, 0, 0, 3,
		1, 0, 2, 0,
		1, 0, 3, 1,

		2, 1, 1, 0,
		2, 1, 2, 1,
		2, 1, 3, 2,
		2, 1, 0, 3,
		2, 1, 2, 0,
		2, 1, 3, 1,

		3, 2, 1, 0,
		3, 2, 2, 1,
		3, 2, 3, 2,
		3, 2, 0, 3,
		3, 2, 2, 0,
		3, 2, 3, 1,

		0, 3, 1, 0,
		0, 3, 2, 1,
		0, 3, 3, 2,
		0, 3, 0, 3,
		0, 3, 2, 0,
		0, 3, 3, 1,

		2, 0, 1, 0,
		2, 0, 2, 1,
		2, 0, 3, 2,
		2, 0, 0, 3,
		2, 0, 2, 0,
		2, 0, 3, 1,

		3, 1, 1, 0,
		3, 1, 2, 1,
		3, 1, 3, 2,
		3, 1, 0, 3,
		3, 1, 2, 0,
		3, 1, 3, 1
	};

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

	const fvec3* const list[8] = { &a0, &a1, &a2, &a3, &b0, &b1, &b2, &b3 };

	//

	float dist     = 99999.0;
	int32_t status = -1;
	float pt[8];
	// 0 - 63 face v.s. vertex
	//64 -135 (0-71)
	//

	//min,maxを2回計算している todo
	for (uint32_t i = 0; i < 4; i++) {
		//aのface v.s. bのvertex
		if (normallist[i].sqlength() < 0.00001)
			continue;
		normallist[i] = normallist[i].normalize();

		const uint32_t& Inda0 = normalindlist[4 * i + 0];
		const uint32_t& Inda1 = normalindlist[4 * i + 1];
		const uint32_t& Inda2 = normalindlist[4 * i + 2];
		const uint32_t& Inda3 = normalindlist[4 * i + 3];

		bool is_fliped = false;

		//法線を内向きにする
		if ((*list[Inda3]).dot(normallist[i]) < (*list[Inda0]).dot(normallist[i])) {
			normallist[i] = -1.0 * normallist[i];
			is_fliped     = true;
		}

		const float a0on = a0.dot(normallist[i]);
		const float a1on = a1.dot(normallist[i]);
		const float a2on = a2.dot(normallist[i]);
		const float a3on = a3.dot(normallist[i]);
		const float b0on = b0.dot(normallist[i]);
		const float b1on = b1.dot(normallist[i]);
		const float b2on = b2.dot(normallist[i]);
		const float b3on = b3.dot(normallist[i]);

		const float maxa = std::max(std::max(std::max(a0on, a1on), a2on), a3on);
		const float mina = std::min(std::min(std::min(a0on, a1on), a2on), a3on);

		const float maxb = std::max(std::max(std::max(b0on, b1on), b2on), b3on);
		const float minb = std::min(std::min(std::min(b0on, b1on), b2on), b3on);

		if (maxb < mina + epsilon)
			return { 0.0, -1 };

		if (maxb - mina < dist) {

			dist = maxb - mina;

			if (
			    b0on > b1on - epsilon && b0on > b2on - epsilon && b0on > b3on - epsilon) {
				status = 8 * i + 0;
				if (is_fliped)
					status += 4;

				auto [dist, p]	    = DistTPlanePoint(*list[Inda0], *list[Inda1], *list[Inda2], b0);
				pt[Inda0]	    = p[0];
				pt[Inda1]	    = p[1];
				pt[Inda2]	    = p[2];
				pt[(Inda2 + 1) % 4] = 0.0;
				pt[4 + 0]	    = 1.0;
				pt[4 + 1]	    = 0.0;
				pt[4 + 2]	    = 0.0;
				pt[4 + 3]	    = 0.0;

			} else if (
			    b1on > b2on - epsilon && b1on > b3on - epsilon && b1on > b0on - epsilon) {
				status = 8 * i + 1;
				if (is_fliped)
					status += 4;

				auto [dist, p]	    = DistTPlanePoint(*list[Inda0], *list[Inda1], *list[Inda2], b1);
				pt[Inda0]	    = p[0];
				pt[Inda1]	    = p[1];
				pt[Inda2]	    = p[2];
				pt[(Inda2 + 1) % 4] = 0.0;
				pt[4 + 0]	    = 0.0;
				pt[4 + 1]	    = 1.0;
				pt[4 + 2]	    = 0.0;
				pt[4 + 3]	    = 0.0;
			} else if (
			    b2on > b3on + epsilon && b2on > b0on - epsilon && b2on > b1on - epsilon) {
				status = 8 * i + 2;
				if (is_fliped)
					status += 4;

				auto [dist, p]	    = DistTPlanePoint(*list[Inda0], *list[Inda1], *list[Inda2], b2);
				pt[Inda0]	    = p[0];
				pt[Inda1]	    = p[1];
				pt[Inda2]	    = p[2];
				pt[(Inda2 + 1) % 4] = 0.0;
				pt[4 + 0]	    = 0.0;
				pt[4 + 1]	    = 0.0;
				pt[4 + 2]	    = 1.0;
				pt[4 + 3]	    = 0.0;
			} else {
				status = 8 * i + 3;
				if (is_fliped)
					status += 4;

				auto [dist, p]	    = DistTPlanePoint(*list[Inda0], *list[Inda1], *list[Inda2], b3);
				pt[Inda0]	    = p[0];
				pt[Inda1]	    = p[1];
				pt[Inda2]	    = p[2];
				pt[(Inda2 + 1) % 4] = 0.0;
				pt[4 + 0]	    = 0.0;
				pt[4 + 1]	    = 0.0;
				pt[4 + 2]	    = 0.0;
				pt[4 + 3]	    = 1.0;
			}
		}
	}

	for (uint32_t i = 4; i < 8; i++) {
		//bのface v.s. aのvertex
		if (normallist[i].sqlength() < 0.00001)
			continue;
		normallist[i] = normallist[i].normalize();

		const uint32_t& Indb0 = normalindlist[4 * i + 0];
		const uint32_t& Indb1 = normalindlist[4 * i + 1];
		const uint32_t& Indb2 = normalindlist[4 * i + 2];
		const uint32_t& Indb3 = normalindlist[4 * i + 3];

		bool is_fliped = false;

		//法線を内向きにする
		if ((*list[4 + Indb3]).dot(normallist[i]) < (*list[4 + Indb0]).dot(normallist[i])) {
			normallist[i] = -1.0 * normallist[i];
			is_fliped     = true;
		}

		const float a0on = a0.dot(normallist[i]);
		const float a1on = a1.dot(normallist[i]);
		const float a2on = a2.dot(normallist[i]);
		const float a3on = a3.dot(normallist[i]);
		const float b0on = b0.dot(normallist[i]);
		const float b1on = b1.dot(normallist[i]);
		const float b2on = b2.dot(normallist[i]);
		const float b3on = b3.dot(normallist[i]);

		const float maxa = std::max(std::max(std::max(a0on, a1on), a2on), a3on);
		const float mina = std::min(std::min(std::min(a0on, a1on), a2on), a3on);
		const float maxb = std::max(std::max(std::max(b0on, b1on), b2on), b3on);
		const float minb = std::min(std::min(std::min(b0on, b1on), b2on), b3on);

		if (maxa < minb + epsilon)
			return { 0.0, -1 };

		if (maxa - minb < dist) {

			dist = maxa - minb;

			if (
			    a0on > a1on - epsilon && a0on > a2on - epsilon && a0on > a3on - epsilon) {
				status = 8 * i + 0;
				if (is_fliped)
					status += 4;

				auto [dist, p]		= DistTPlanePoint(*list[4 + Indb0], *list[4 + Indb1], *list[4 + Indb2], a0);
				pt[0]			= 1.0;
				pt[1]			= 0.0;
				pt[2]			= 0.0;
				pt[3]			= 0.0;
				pt[4 + Indb0]		= p[0];
				pt[4 + Indb1]		= p[1];
				pt[4 + Indb2]		= p[2];
				pt[4 + (Indb2 + 1) % 4] = 0.0;
			} else if (
			    a1on > a2on - epsilon && a1on > a3on - epsilon && a1on > a0on - epsilon) {
				status = 8 * i + 1;
				if (is_fliped)
					status += 4;

				auto [dist, p]		= DistTPlanePoint(*list[4 + Indb0], *list[4 + Indb1], *list[4 + Indb2], a1);
				pt[0]			= 0.0;
				pt[1]			= 1.0;
				pt[2]			= 0.0;
				pt[3]			= 0.0;
				pt[4 + Indb0]		= p[0];
				pt[4 + Indb1]		= p[1];
				pt[4 + Indb2]		= p[2];
				pt[4 + (Indb2 + 1) % 4] = 0.0;
			} else if (
			    a2on > a3on + epsilon && a2on > a0on - epsilon && a2on > a1on - epsilon) {
				status = 8 * i + 2;
				if (is_fliped)
					status += 4;

				auto [dist, p]		= DistTPlanePoint(*list[4 + Indb0], *list[4 + Indb1], *list[4 + Indb2], a2);
				pt[0]			= 0.0;
				pt[1]			= 0.0;
				pt[2]			= 1.0;
				pt[3]			= 0.0;
				pt[4 + Indb0]		= p[0];
				pt[4 + Indb1]		= p[1];
				pt[4 + Indb2]		= p[2];
				pt[4 + (Indb2 + 1) % 4] = 0.0;
			} else {
				status = 8 * i + 3;
				if (is_fliped)
					status += 4;

				auto [dist, p]		= DistTPlanePoint(*list[4 + Indb0], *list[4 + Indb1], *list[4 + Indb2], a3);
				pt[0]			= 0.0;
				pt[1]			= 0.0;
				pt[2]			= 0.0;
				pt[3]			= 1.0;
				pt[4 + Indb0]		= p[0];
				pt[4 + Indb1]		= p[1];
				pt[4 + Indb2]		= p[2];
				pt[4 + (Indb2 + 1) % 4] = 0.0;
			}
		}
	}

	//同じ法線を生成する組がある場合に，適切な最短の辺の組を計算できない
	//iに応じて辺を構成する頂点のみ選択する
	for (uint32_t i = 8; i < 44; i++) {

		//edge v.s. edge
		//サーチする際に浸透が最も小さい組がedge v.s. edgeの場合である．
		//gjkでミンコフスキー差の外周に現れる頂点がどちらの組で構成されるかの違いである
		//normallist[i]を基準にしてa,bの辺のどちらが上にあるかを確認する

		const uint32_t Inda0 = normalindlist[4 * i + 0];
		const uint32_t Inda1 = normalindlist[4 * i + 1];
		const uint32_t Indb0 = normalindlist[4 * i + 2];
		const uint32_t Indb1 = normalindlist[4 * i + 3];

		//a0 == a1
		//b0 == b1
		//(a0,a1) // (b0,b1)
		if (normallist[i].sqlength() < 0.00001) {

			if (((*list[Inda1]) - (*list[Inda0])).sqlength() > epsilon && ((*list[Indb1]) - (*list[Indb0])).sqlength() > epsilon) {

				fvec3 n = (*list[Indb0]) - (*list[Inda0]);
				if (n.sqlength() < epsilon)
					continue;
				fvec3 dir     = (*list[Inda1]) - (*list[Inda0]);
				normallist[i] = n - (n.dot(dir) / dir.sqlength()) * dir;
			} else {
				continue;
			}
		}

		normallist[i] = normallist[i].normalize();

		const float a0on = a0.dot(normallist[i]);
		const float a1on = a1.dot(normallist[i]);
		const float a2on = a2.dot(normallist[i]);
		const float a3on = a3.dot(normallist[i]);
		const float b0on = b0.dot(normallist[i]);
		const float b1on = b1.dot(normallist[i]);
		const float b2on = b2.dot(normallist[i]);
		const float b3on = b3.dot(normallist[i]);

		float maxa;
		float mina;
		float maxb;
		float minb;

		uint32_t Indmina, Indmaxa, Indminb, Indmaxb;

		if (a0on < a1on + epsilon && a0on < a2on + epsilon && a0on < a3on + epsilon) {
			Indmina = 0;
			mina	= a0on;
		} else if (a1on < a2on + epsilon && a1on < a3on + epsilon && a1on < a0on + epsilon) {
			Indmina = 1;
			mina	= a1on;
		} else if (a2on < a3on + epsilon && a2on < a0on + epsilon && a2on < a1on + epsilon) {
			Indmina = 2;
			mina	= a2on;
		} else {
			Indmina = 3;
			mina	= a3on;
		}

		if (a0on > a1on - epsilon && a0on > a2on - epsilon && a0on > a3on - epsilon) {
			Indmaxa = 0;
			maxa	= a0on;
		} else if (a1on > a2on - epsilon && a1on > a3on - epsilon && a1on > a0on - epsilon) {
			Indmaxa = 1;
			maxa	= a1on;
		} else if (a2on > a3on - epsilon && a2on > a0on - epsilon && a2on > a1on - epsilon) {
			Indmaxa = 2;
			maxa	= a2on;
		} else {
			Indmaxa = 3;
			maxa	= a3on;
		}

		if (b0on < b1on + epsilon && b0on < b2on + epsilon && b0on < b3on + epsilon) {
			Indminb = 0;
			minb	= b0on;
		} else if (b1on < b2on + epsilon && b1on < b3on + epsilon && b1on < b0on + epsilon) {
			Indminb = 1;
			minb	= b1on;
		} else if (b2on < b3on + epsilon && b2on < b0on + epsilon && b2on < b1on + epsilon) {
			Indminb = 2;
			minb	= b2on;
		} else {
			Indminb = 3;
			minb	= b3on;
		}

		if (b0on > b1on - epsilon && b0on > b2on - epsilon && b0on > b3on - epsilon) {
			Indmaxb = 0;
			maxb	= b0on;
		} else if (b1on > b2on - epsilon && b1on > b3on - epsilon && b1on > b0on - epsilon) {
			Indmaxb = 1;
			maxb	= b1on;
		} else if (b2on > b3on - epsilon && b2on > b0on - epsilon && b2on > b1on - epsilon) {
			Indmaxb = 2;
			maxb	= b2on;
		} else {
			Indmaxb = 3;
			maxb	= b3on;
		}

		if (maxa - epsilon < minb || maxb - epsilon < mina)
			return { 0.0, -1 };

		if (maxa - minb < maxb - mina) {
			//bのほうが高い
			if (maxa - minb < dist) {
				if ((Inda0 == Indmaxa || Inda1 == Indmaxa) && (Indb0 == Indminb || Indb1 == Indminb)) {
					dist   = maxa - minb;
					status = 64 + 2 * (i - 8) + 0;

					pt[0] = 0.0;
					pt[1] = 0.0;
					pt[2] = 0.0;
					pt[3] = 0.0;
					pt[4] = 0.0;
					pt[5] = 0.0;
					pt[6] = 0.0;
					pt[7] = 0.0;

					auto [hoge, p] = DistLineLine(*list[Inda0], *list[Inda1], *list[4 + Indb0], *list[4 + Indb1]);
					pt[Inda0]      = p[0];
					pt[Inda1]      = 1.0f - p[0];
					pt[4 + Indb0]  = p[1];
					pt[4 + Indb1]  = 1.0f - p[1];
				}
			}
		} else {
			//aのほうが高い
			if (maxb - mina < dist) {
				if ((Inda0 == Indmina || Inda1 == Indmina) && (Indb0 == Indmaxb || Indb1 == Indmaxb)) {
					dist   = maxb - mina;
					status = 64 + 2 * (i - 8) + 1;

					pt[0] = 0.0;
					pt[1] = 0.0;
					pt[2] = 0.0;
					pt[3] = 0.0;
					pt[4] = 0.0;
					pt[5] = 0.0;
					pt[6] = 0.0;
					pt[7] = 0.0;

					auto [hoge, p] = DistLineLine(*list[Inda0], *list[Inda1], *list[4 + Indb0], *list[4 + Indb1]);
					pt[Inda0]      = p[0];
					pt[Inda1]      = 1.0f - p[0];
					pt[4 + Indb0]  = p[1];
					pt[4 + Indb1]  = 1.0f - p[1];
				}
			}
		}
	}

	return { dist, status, pt[0], pt[1], pt[2], pt[3], pt[4], pt[5], pt[6], pt[7] };
}
