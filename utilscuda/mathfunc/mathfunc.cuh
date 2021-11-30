#pragma once

//TODO forceinlineをつけるべきか

#include <cstdint>

class cudavec3;
class cudamat3;

class cudavec3 {
    public:
	float x, y, z;

	__host__ __device__ cudavec3(const float& x, const float& y, const float& z)
	    : x(x)
	    , y(y)
	    , z(z)
	{
	}
	__host__ __device__ cudavec3(const float3& v)
	    : x(v.x)
	    , y(v.y)
	    , z(v.z)
	{
	}

	__host__ __device__ cudavec3(void)
	    : x(0.0)
	    , y(0.0)
	    , z(0.0)
	{
	}

	__host__ __device__ cudavec3(const float (&array)[3])
	    : x(array[0])
	    , y(array[1])
	    , z(array[2])
	{
	}
	__host__ __device__ float dot(const cudavec3& v) const;
	__host__ __device__ cudavec3 cross(const cudavec3& a) const;
	__host__ __device__ float sqlength() const;
	__host__ __device__ float length() const;
	__host__ __device__ cudavec3 normalize() const;
	__host__ __device__ cudamat3 tensorproduct(const cudavec3& v) const;
	__host__ __device__ cudamat3 skew() const;
	__host__ __device__ cudamat3 rotation() const;
	__host__ __device__ static float STP(const cudavec3& a, const cudavec3& b, const cudavec3& c);
};

__host__ __device__ const cudavec3 operator+(const cudavec3& a, const cudavec3& b);
__host__ __device__ const cudavec3 operator-(const cudavec3& a, const cudavec3& b);
__host__ __device__ const cudavec3 operator*(const float& a, const cudavec3& v);
__host__ __device__ const cudavec3 operator*(const cudavec3& v, const float& a);
__host__ __device__ const cudavec3 operator/(const cudavec3& v, const float& a);

//////////////////////////////////
//////////////////////////////////
//////////////////////////////////

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

	__host__ __device__ float det(void) const;
	__host__ __device__ cudamat3 inverse(void) const;
	__host__ __device__ float sqlength() const;
	__host__ __device__ float trace() const;
	__host__ __device__ cudamat3 transpose() const;
	__host__ __device__ static cudamat3 identity();
};

__host__ __device__ const cudamat3 operator+(const cudamat3& mat0, const cudamat3& mat1);
__host__ __device__ const cudamat3 operator-(const cudamat3& mat0, const cudamat3& mat1);
__host__ __device__ const cudamat3 operator-(const cudamat3& mat);
__host__ __device__ const cudamat3 operator*(const cudamat3& mat, const float& a);
__host__ __device__ const cudamat3 operator*(const float& a, const cudamat3& mat);
__host__ __device__ const cudamat3 operator*(const cudamat3& mat0, const cudamat3& mat1);

//////////////////////////////////
//////////////////////////////////
//////////////////////////////////

class cudaquaternion {
    public:
	float x, y, z, w;

	__host__ __device__ cudaquaternion(
	    const float& x,
	    const float& y,
	    const float& z,
	    const float& w)
	    : x(x)
	    , y(y)
	    , z(z)
	    , w(w)
	{
	}
	__host__ __device__ cudaquaternion(const cudavec3& v, const float& s)
	    : x(v.x)
	    , y(v.y)
	    , z(v.z)
	    , w(s)
	{
	}
	__host__ __device__ cudaquaternion(void)
	    : x(0.0)
	    , y(0.0)
	    , z(0.0)
	    , w(1.0)
	{
	}
	__host__ __device__ cudaquaternion(const cudavec3& v);

	__host__ __device__ cudavec3 getv() const;
	__host__ __device__ cudaquaternion conjugate() const;
	__host__ __device__ float dot(const cudaquaternion& q) const;
	__host__ __device__ float sqlength() const;
	__host__ __device__ float length() const;
	__host__ __device__ cudaquaternion normalize() const;
	__host__ __device__ cudamat3 qtomat() const;
	__host__ __device__ cudaquaternion inverse() const;

	__host__ __device__ static cudaquaternion slerp(const cudaquaternion& q0, const cudaquaternion& q1, const float& t);
};

__host__ __device__ const cudaquaternion operator+(const cudaquaternion& q0, const cudaquaternion& q1);
__host__ __device__ const cudaquaternion operator-(const cudaquaternion& q0, const cudaquaternion& q1);
__host__ __device__ const cudaquaternion operator-(const cudaquaternion& q);
__host__ __device__ const cudaquaternion operator*(const cudaquaternion& q0, const cudaquaternion& q1);
__host__ __device__ const cudaquaternion operator*(const cudaquaternion& q, const float& a);
__host__ __device__ const cudaquaternion operator*(const float& a, const cudaquaternion& q);
__host__ __device__ const cudaquaternion operator/(const cudaquaternion& q, const float& a);

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////↓  実装///////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

__host__ __device__ float cudavec3::dot(const cudavec3& v) const
{
	return x * v.x + y * v.y + z * v.z;
}
__host__ __device__ cudavec3 cudavec3::cross(const cudavec3& a) const
{
	return cudavec3(y * a.z - z * a.y, z * a.x - x * a.z, x * a.y - y * a.x);
}
__host__ __device__ float cudavec3::sqlength() const
{
	return x * x + y * y + z * z;
}
__host__ __device__ float cudavec3::length() const
{
	return sqrtf(this->sqlength());
}
__host__ __device__ cudavec3 cudavec3::normalize() const
{
	float length = this->length();
	if (length < 0.0000001)
		return cudavec3();
	return (*this) / length;
}
__host__ __device__ cudamat3 cudavec3::tensorproduct(const cudavec3& v) const
{
	cudamat3 hoge;
	hoge.m[0 * 3 + 0] = x * v.x;
	hoge.m[0 * 3 + 1] = x * v.y;
	hoge.m[0 * 3 + 2] = x * v.z;
	hoge.m[1 * 3 + 0] = y * v.x;
	hoge.m[1 * 3 + 1] = y * v.y;
	hoge.m[1 * 3 + 2] = y * v.z;
	hoge.m[2 * 3 + 0] = z * v.x;
	hoge.m[2 * 3 + 1] = z * v.y;
	hoge.m[2 * 3 + 2] = z * v.z;
	return hoge;
}
__host__ __device__ cudamat3 cudavec3::skew() const
{
	cudamat3 hoge;
	hoge.m[0 * 3 + 0] = 0.0;
	hoge.m[0 * 3 + 1] = -z;
	hoge.m[0 * 3 + 2] = y;
	hoge.m[1 * 3 + 0] = z;
	hoge.m[1 * 3 + 1] = 0.0;
	hoge.m[1 * 3 + 2] = -x;
	hoge.m[2 * 3 + 0] = -y;
	hoge.m[2 * 3 + 1] = x;
	hoge.m[2 * 3 + 2] = 0.0;
	return hoge;
}
__host__ __device__ cudamat3 cudavec3::rotation() const
{
	float omega  = this->length();
	float omega2 = this->sqlength();
	if (omega2 < 0.0000001)
		return cudamat3::identity();
	float cos = cosf(omega);
	float sin = sinf(omega);

	return cos * cudamat3::identity() + ((1 - cos) / omega2) * this->tensorproduct(*this) + (sin / omega) * this->skew();
}
__host__ __device__ float STP(const cudavec3& a, const cudavec3& b, const cudavec3& c)
{
	return a.dot(b.cross(c));
}

__host__ __device__ const cudavec3 operator+(const cudavec3& a, const cudavec3& b)
{
	return cudavec3(
	    a.x + b.x,
	    a.y + b.y,
	    a.z + b.z);
}
__host__ __device__ const cudavec3 operator-(const cudavec3& a, const cudavec3& b)
{
	return cudavec3(
	    a.x - b.x,
	    a.y - b.y,
	    a.z - b.z);
}
__host__ __device__ const cudavec3 operator*(const float& a, const cudavec3& v)
{
	return cudavec3(
	    a * v.x,
	    a * v.y,
	    a * v.z);
}
__host__ __device__ const cudavec3 operator*(const cudavec3& v, const float& a)
{
	return cudavec3(
	    a * v.x,
	    a * v.y,
	    a * v.z);
}
__host__ __device__ const cudavec3 operator/(const cudavec3& v, const float& a)
{
	return cudavec3(
	    v.x / a,
	    v.y / a,
	    v.z / a);
}

//////////////////////////////////
//////////////////////////////////
//////////////////////////////////

__host__ __device__ float cudamat3::det(void) const
{

	float temp[3];
	temp[0] = m[1] * m[5] - m[2] * m[4];
	temp[1] = m[2] * m[3] - m[0] * m[5];
	temp[2] = m[0] * m[4] - m[1] * m[3];

	return m[6] * temp[0] + m[7] * temp[1] + m[8] * temp[2];
}
__host__ __device__ cudamat3 cudamat3::inverse(void) const
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
__host__ __device__ float cudamat3::sqlength() const
{
	float hoge = 0.0;
	for (uint32_t i = 0; i < 9; i++)
		hoge += m[i];
	return hoge;
}
__host__ __device__ float cudamat3::trace() const
{
	return m[0] + m[4] + m[8];
}
__host__ __device__ cudamat3 cudamat3::transpose() const
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
__host__ __device__ cudamat3 cudamat3::identity()
{
	cudamat3 temp;
	temp.m[0] = 1.0;
	temp.m[4] = 1.0;
	temp.m[8] = 1.0;
	return temp;
}

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

//////////////////////////////////
//////////////////////////////////
//////////////////////////////////

__host__ __device__ cudaquaternion::cudaquaternion(const cudavec3& v)
{
	float Omega = v.length();
	cudavec3 l  = v.normalize();

	float sin = sinf(Omega / 2.0);
	x	  = sin * l.x;
	y	  = sin * l.y;
	z	  = sin * l.z;
	w	  = cosf(Omega / 2.0);
}

__host__ __device__ cudavec3 cudaquaternion::getv() const
{
	return cudavec3(x, y, z);
}
__host__ __device__ cudaquaternion cudaquaternion::conjugate() const
{
	return cudaquaternion(-1.0 * cudavec3(x, y, z), w);
}
__host__ __device__ float cudaquaternion::dot(const cudaquaternion& q) const
{
	return x * q.x + y * q.y + z * q.z + w * q.w;
}
__host__ __device__ float cudaquaternion::sqlength() const
{
	return x * x + y * y + z * z + w * w;
}
__host__ __device__ float cudaquaternion::length() const
{
	return sqrtf(x * x + y * y + z * z + w * w);
}
__host__ __device__ cudaquaternion cudaquaternion::normalize() const
{
	float length = this->length();
	if (length < 0.0000001)
		return cudaquaternion();
	return (*this) / length;
}
__host__ __device__ cudamat3 cudaquaternion::qtomat() const
{
	cudavec3 v = this->getv();
	float w	   = this->w;

	return (w * w - v.sqlength()) * cudamat3::identity() + 2.0 * w * v.skew() + 2.0 * v.tensorproduct(v);
}
__host__ __device__ cudaquaternion cudaquaternion::inverse() const
{
	return cudaquaternion(-x, -y, -z, w);
}

__host__ __device__ cudaquaternion slerp(const cudaquaternion& q0, const cudaquaternion& q1, const float& t)
{
	if (q0.dot(q1) < 0.0) {
		double x = 0.5 * acosf(-q0.dot(q1));
		return (sinf(t * x) / sinf(x)) * -q1 + (sinf(x - t * x) / sinf(x)) * q0;
	} else {
		double x = 0.5 * acosf(q0.dot(q1));
		return (sinf(t * x) / sinf(x)) * q1 + (sinf(x - t * x) / sinf(x)) * q0;
	}
}

__host__ __device__ const cudaquaternion operator+(const cudaquaternion& q0, const cudaquaternion& q1)
{
	return cudaquaternion(
	    q0.x + q1.x,
	    q0.y + q1.y,
	    q0.z + q1.z,
	    q0.w + q1.w);
}
__host__ __device__ const cudaquaternion operator-(const cudaquaternion& q0, const cudaquaternion& q1)
{
	return cudaquaternion(
	    q0.x - q1.x,
	    q0.y - q1.y,
	    q0.z - q1.z,
	    q0.w - q1.w);
}
__host__ __device__ const cudaquaternion operator-(const cudaquaternion& q)
{
	return cudaquaternion(
	    -q.x,
	    -q.y,
	    -q.z,
	    -q.w);
}
__host__ __device__ const cudaquaternion operator*(const cudaquaternion& q0, const cudaquaternion& q1)
{
	float w	   = q0.w * q1.w - (q0.x * q1.x + q0.y * q1.y + q0.z * q1.z);
	cudavec3 v = q0.w * q1.getv() + q1.w * q0.getv() + q0.getv().cross(q1.getv());

	return cudaquaternion(v, w);
}
__host__ __device__ const cudaquaternion operator*(const cudaquaternion& q, const float& a)
{
	return cudaquaternion(
	    a * q.x,
	    a * q.y,
	    a * q.z,
	    a * q.w);
}
__host__ __device__ const cudaquaternion operator*(const float& a, const cudaquaternion& q)
{
	return cudaquaternion(
	    a * q.x,
	    a * q.y,
	    a * q.z,
	    a * q.w);
}
__host__ __device__ const cudaquaternion operator/(const cudaquaternion& q, const float& a)
{
	return cudaquaternion(
	    q.x / a,
	    q.y / a,
	    q.z / a,
	    q.w / a);
}
