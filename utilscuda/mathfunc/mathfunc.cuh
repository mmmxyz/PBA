#pragma once

#include <cstdint>

class cudavec3 {
    public:
	float x, y, z;
};

class cudamat3 {
    public:
	float m[9];
	//m[0],m[1],m[2];
	//m[3],m[4],m[5];
	//m[6],m[7],m[8];

	__host__ __device__ cudamat3(void)
	{
		for (uint32_t i = 0; i < 9; i++)
			m[i] = 0.0;
	}

	__host__ __device__ cudamat3(const float3& a, const float3& b, const float3& c)
	{
		m[3 * 0 + 0] = a.x;
		m[3 * 1 + 0] = a.y;
		m[3 * 2 + 0] = a.z;

		m[3 * 0 + 1] = b.x;
		m[3 * 1 + 1] = b.y;
		m[3 * 2 + 1] = b.z;

		m[3 * 0 + 2] = c.x;
		m[3 * 1 + 2] = c.y;
		m[3 * 2 + 2] = c.z;
	}

	__host__ __device__ cudamat3(const float (&array)[9])
	{
		for (uint32_t i = 0; i < 9; i++)
			m[i] = 0.0;
	}

	__host__ __device__ float det(void) const
	{

		float temp[3];
		temp[0] = m[1] * m[5] - m[2] * m[4];
		temp[1] = m[2] * m[3] - m[0] * m[5];
		temp[2] = m[0] * m[4] - m[1] * m[3];

		return m[6] * temp[0] + m[7] * temp[1] + m[8] * temp[2];
	}

	__host__ __device__ cudamat3 inverse(void) const
	{
		float det = this->det();
		float data[9];

		data[0] = (m[4] * m[8] - m[5] * m[7]) / det;
		data[1] = -(m[1] * m[8] - m[2] * m[7]) / det;
		data[2] = (m[1] * m[5] - m[2] * m[4]) / det;
		data[3] = -(m[3] * m[8] - m[5] * m[6]) / det;
		data[4] = (m[0] * m[8] - m[2] * m[6]) / det;
		data[5] = -(m[0] * m[5] - m[2] * m[3]) / det;
		data[6] = (m[3] * m[7] - m[4] * m[6]) / det;
		data[7] = -(m[0] * m[7] - m[1] * m[6]) / det;
		data[8] = (m[0] * m[4] - m[1] * m[3]) / det;

		return cudamat3(data);
	}

	__host__ __device__ float sqlength() const
	{
		float hoge = 0.0;
		for (uint32_t i = 0; i < 9; i++)
			hoge += m[i];
		return hoge;
	}

	__host__ __device__ float trace() const
	{
		return m[0] + m[4] + m[8];
	}

	__host__ __device__ cudamat3 transpose() const
	{
		float data[9];
		data[0] = m[0];
		data[1] = m[3];
		data[2] = m[6];
		data[3] = m[1];
		data[4] = m[4];
		data[5] = m[7];
		data[6] = m[2];
		data[7] = m[5];
		data[8] = m[8];

		return cudamat3(data);
	}

	__host__ __device__ static cudamat3 identity()
	{
		cudamat3 temp;
		temp.m[0] = 1.0;
		temp.m[4] = 1.0;
		temp.m[8] = 1.0;
		return temp;
	}
};

__host__ __device__ const cudamat3 operator+(const cudamat3& mat0, const cudamat3& mat1)
{
	float p[9];
	for (uint32_t i = 0; i < 9; i++)
		p[i] = mat0.m[i] + mat1.m[i];
	return cudamat3(p);
}
__host__ __device__ const cudamat3 operator-(const cudamat3& mat0, const cudamat3& mat1)
{
	float p[9];
	for (uint32_t i = 0; i < 9; i++)
		p[i] = mat0.m[i] - mat1.m[i];
	return cudamat3(p);
}
__host__ __device__ const cudamat3 operator-(const cudamat3& mat)
{
	float p[9];
	for (uint32_t i = 0; i < 9; i++)
		p[i] = -mat.m[i];
	return cudamat3(p);
}
__host__ __device__ const cudamat3 operator*(const cudamat3& mat, const float& a)
{
	float p[9];
	for (uint32_t i = 0; i < 9; i++)
		p[i] = a * mat.m[i];
	return cudamat3(p);
}
__host__ __device__ const cudamat3 operator*(const float& a, const cudamat3& mat)
{
	float p[9];
	for (uint32_t i = 0; i < 9; i++)
		p[i] = a * mat.m[i];
	return cudamat3(p);
}
__host__ __device__ const cudamat3 operator*(const cudamat3& mat0, const cudamat3& mat1)
{
	float p[9];
	for (uint32_t i = 0; i < 3; i++) {
		for (uint32_t j = 0; j < 3; j++) {
			p[i * 3 + j] = mat0.m[i * 3 + 0] * mat1.m[0 * 3 + j] + mat0.m[i * 3 + 1] * mat1.m[1 * 3 + j] + mat0.m[i * 3 + 2] * mat1.m[2 * 3 + j];
		}
	}
	return cudamat3(p);
}
