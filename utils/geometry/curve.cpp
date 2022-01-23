#include <cmath>

#include "utils/geometry/curve.hpp"

constexpr float epsilon = 0.0000001;

fvec3 bezier_curve::operator()(const float& t)
{
	fvec3* Array0 = new fvec3[N];
	fvec3* Array1 = new fvec3[N];

	for (uint32_t i = 0; i < N; i++) {
		Array0[i] = (1 - t) * P[i] + t * P[i + 1];
	}

	for (uint32_t i = 2; i < N + 1; i++) {
		for (uint32_t j = 0; j < N + 1 - i; j++) {

			if (i % 2 == 1)
				Array0[j] = (1 - t) * Array1[j] + t * Array1[j + 1];
			else
				Array1[j] = (1 - t) * Array0[j] + t * Array0[j + 1];
		}
	}

	fvec3 ret;
	if (N % 2 == 0)
		ret = Array0[0];
	else
		ret = Array1[0];

	delete[] Array0;
	delete[] Array1;

	return ret;
}

fvec3 b_spline::operator()(const float& t)
{

	uint32_t i = 0;
	if (t <= knots[p] + epsilon)
		i = p;
	else if (t + epsilon >= knots[n + 1])
		i = n;
	else {
		for (uint32_t j = p; j <= n; j++)
			if (t <= knots[j + 1] + epsilon && knots[j] <= t + epsilon)
				i = j;
	}

	fvec3* Array0 = new fvec3[p];
	fvec3* Array1 = new fvec3[p];

	for (uint32_t j = i - p + 1; j <= i; j++) {
		if (std::abs(knots[j + p] - knots[j]) < epsilon)
			Array0[j - (i - p + 1)] = fvec3(0.0);
		else {
			float alpha		= (t - knots[j]) / (knots[j + p] - knots[j]);
			Array0[j - (i - p + 1)] = alpha * P[j] + (1.0 - alpha) * P[j - 1];
		}
	}

	for (uint32_t k = 2; k < p + 1; k++) {
		for (uint32_t j = i - (p - k); j <= i; j++) {
			if (k % 2 == 1) {

				if (std::abs(knots[j + p - k + 1] - knots[j]) < epsilon)
					Array0[j - (i - p + 1)] = fvec3(0.0);
				else {
					float alpha		= (t - knots[j]) / (knots[j + p - k + 1] - knots[j]);
					Array0[j - (i - p + 1)] = alpha * Array1[j - (i - p + 1)] + (1.0 - alpha) * Array1[j - 1 - (i - p + 1)];
				}
			} else {

				if (std::abs(knots[j + p - k + 1] - knots[j]) < epsilon)
					Array1[j - (i - p + 1)] = fvec3(0.0);
				else {
					float alpha		= (t - knots[j]) / (knots[j + p - k + 1] - knots[j]);
					Array1[j - (i - p + 1)] = alpha * Array0[j - (i - p + 1)] + (1.0 - alpha) * Array0[j - 1 - (i - p + 1)];
				}
			}
		}
	}

	fvec3 ret;
	if (p % 2 == 0)
		ret = Array1[p - 1];
	else
		ret = Array0[p - 1];

	delete[] Array0;
	delete[] Array1;

	return ret;
}

b_spline b_spline::diff()
{
	fvec3* Q = new fvec3[n];
	for (uint32_t i = 0; i < n; i++) {
		if (std::abs(knots[i + p + 1] - knots[i + 1]) < epsilon)
			Q[i] = fvec3(0.0);
		else {
			Q[i] = p * (P[i + 1] - P[i]) / (knots[i + p + 1] - knots[i + 1]);
		}
	}

	b_spline diff(p - 1, n - 1, Q, knots + 1);

	delete[] Q;

	return diff;
}

fvec3 piecewise_cubic_hermite::operator()(const float& t)
{
	int32_t i = t * float(n);
	if (i < 0)
		i = 0;
	if (i >= n)
		i = n - 1;

	float alpha = t - i / float(n);
	alpha *= float(n);

	return (1.0 - alpha) * (1.0 - alpha) * (1.0 - alpha) * P[i] + 3.0 * alpha * (1.0 - alpha) * (1.0 - alpha) * (P[i] + V[i] / 3.0) + 3.0 * alpha * alpha * (1.0 - alpha) * (P[i + 1] - V[i + 1] / 3.0) + alpha * alpha * alpha * P[i + 1];
}

fvec3 piecewise_cubic_hermite::derive(const float& t)
{
	int32_t i = t * float(n);
	if (i < 0)
		i = 0;
	if (i >= n)
		i = n - 1;

	float alpha = t - i / float(n);
	alpha *= float(n);

	return 3.0 * ((1.0 - alpha) * (1.0 - alpha) * (V[i] / 3.0) + 2.0 * (1.0 - alpha) * alpha * (-V[i + 1] / 3.0 - V[i] / 3.0 - P[i] + P[i + 1]) + alpha * alpha * (V[i + 1] / 3.0));
}

fvec3 catmull_rom::operator()(const float& t)
{
	int32_t i = t * float(n);
	if (i < 0)
		i = 0;
	if (i >= n)
		i = n - 1;

	float alpha = t - i / float(n);
	alpha *= float(n);

	return (1.0 - alpha) * (1.0 - alpha) * (1.0 - alpha) * P[i] + 3.0 * alpha * (1.0 - alpha) * (1.0 - alpha) * (P[i] + V[i] / 3.0) + 3.0 * alpha * alpha * (1.0 - alpha) * (P[i + 1] - V[i + 1] / 3.0) + alpha * alpha * alpha * P[i + 1];
}

fvec3 catmull_rom::derive(const float& t)
{
	int32_t i = t * float(n);
	if (i < 0)
		i = 0;
	if (i >= n)
		i = n - 1;

	float alpha = t - i / float(n);
	alpha *= float(n);

	return 3.0 * ((1.0 - alpha) * (1.0 - alpha) * (V[i] / 3.0) + 2.0 * (1.0 - alpha) * alpha * (-V[i + 1] / 3.0 - V[i] / 3.0 - P[i] + P[i + 1]) + alpha * alpha * (V[i + 1] / 3.0));
}
