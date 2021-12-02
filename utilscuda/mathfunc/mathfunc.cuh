#pragma once

//TODO forceinlineをつけるべきか

//utils/mathfunc/mathfunc.hpp
//float型に特殊化せずに template<class T> vec2<T>::hogehoge みたいにして定義すると当然エラーとなるが
//__host__ __device__をつけて特殊化すると __host__ __device__　コードだけを生成される。なぜかはわからない
//__host__ __host__ __device__をつけて特殊化すると両方のコードだけを生成される。このファイルをincludeしたnvccでコンパイルされるコード(デバイスで動作するコードとカーネルを呼び出すCPUで動作するコード)はこっちの定義を呼び出すが、それ以外の普通のcppのコードはmathfunc.hppの定義の方を呼び出す。
//普通のcppのコードがこちらをincludeしても差し支えないように __CUDA_ARCH__ でガードする。
//gccが __host__ __device__ を見つけると
// error: ‘__host__ __device__’ does not name a type; did you mean ‘__dev_t’?
//とエラーになる。
//ある関数を __host__ __host__ __device__ で異なる定義を別に書くことはできないが、
//特殊化の際に __host__ __host__ __device__ をつければできるということだろうか
//TODOあとでしらべる

#include <cstdint>

#include "utils/mathfunc/mathfunc.hpp"

#if defined(__CUDA_ARCH__)

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////↓  宣言///////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

template <>
__host__ __device__ fvec2::vec2(const float& x, const float& y)
    : x(x)
    , y(y)
{
}
template <>
__host__ __device__ fvec2::vec2(const float& value)
    : x(value)
    , y(value)
{
}
template <>
__host__ __device__ fvec2::vec2(const float (&array)[2])
    : x(array[0])
    , y(array[1])
{
}
template <>
__host__ __device__ fvec2::vec2(void)
    : x(0.0)
    , y(0.0)
{
}

template <>
__host__ __device__ float fvec2::dot(const fvec2& v) const;
template <>
__host__ __device__ float fvec2::cross(const fvec2& v) const;
template <>
__host__ __device__ float fvec2::sqlength() const;
template <>
__host__ __device__ float fvec2::length() const;
template <>
__host__ __device__ fvec2 fvec2::normalize() const;
//template <>
// __host__ __device__ fvec2 fvec2::rot() const;

__host__ __device__ const fvec2 operator+(const fvec2& a, const fvec2& b);
__host__ __device__ const fvec2 operator-(const fvec2& a, const fvec2& b);
__host__ __device__ const fvec2 operator-(const fvec2& v);
__host__ __device__ const fvec2 operator*(const float& a, const fvec2& v);
__host__ __device__ const fvec2 operator*(const fvec2& v, const float& a);
__host__ __device__ const fvec2 operator/(const fvec2& v, const float& a);

//////////////////////////////////
//////////////////////////////////
//////////////////////////////////

template <>
__host__ __device__ fvec3::vec3(const float& x, const float& y, const float& z)
    : x(x)
    , y(y)
    , z(z)
{
}

template <>
__host__ __device__ fvec3::vec3(const float& value)
    : x(value)
    , y(value)
    , z(value)
{
}

template <>
__host__ __device__ fvec3::vec3(void)
    : x(0.0)
    , y(0.0)
    , z(0.0)
{
}

template <>
__host__ __device__ fvec3::vec3(const float (&array)[3])
    : x(array[0])
    , y(array[1])
    , z(array[2])
{
}

template <>
__host__ __device__ float fvec3::dot(const fvec3& v) const;
template <>
__host__ __device__ fvec3 fvec3::cross(const fvec3& a) const;
template <>
__host__ __device__ float fvec3::sqlength() const;
template <>
__host__ __device__ float fvec3::length() const;
template <>
__host__ __device__ fvec3 fvec3::normalize() const;
template <>
__host__ __device__ fmat3 fvec3::tensorproduct(const fvec3& v) const;
template <>
__host__ __device__ fmat3 fvec3::skew() const;
template <>
__host__ __device__ fmat3 fvec3::rotation() const;
template <>
__host__ __device__ float fvec3::STP(const fvec3& a, const fvec3& b, const fvec3& c);

__host__ __device__ const fvec3 operator+(const fvec3& a, const fvec3& b);
__host__ __device__ const fvec3 operator-(const fvec3& a, const fvec3& b);
__host__ __device__ const fvec3 operator-(const fvec3& v);
__host__ __device__ const fvec3 operator*(const float& a, const fvec3& v);
__host__ __device__ const fvec3 operator*(const fvec3& v, const float& a);
__host__ __device__ const fvec3 operator/(const fvec3& v, const float& a);

//////////////////////////////////
//////////////////////////////////
//////////////////////////////////

template <>
__host__ __device__ fvec4::vec4(
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
template <>
__host__ __device__ fvec4::vec4(
    const fvec3& v)
    : x(v.x)
    , y(v.y)
    , z(v.z)
    , w(1.0)
{
}

template <>
__host__ __device__ float fvec4::dot(const fvec4& v) const;
template <>
__host__ __device__ float fvec4::sqlength(void) const;

__host__ __device__ const fvec4 operator+(const fvec4& a, const fvec4& b);
__host__ __device__ const fvec4 operator-(const fvec4& a, const fvec4& b);
__host__ __device__ const fvec4 operator*(const float& a, const fvec4& v);
__host__ __device__ const fvec4 operator*(const fvec4& v, const float& a);

//////////////////////////////////
//////////////////////////////////
//////////////////////////////////

template <>
__host__ __device__ fmat2::mat2(void)
{
	for (uint32_t i = 0; i < 9; i++)
		m[i] = 0.0;
}
template <>
__host__ __device__ fmat2::mat2(const fvec2& a, const fvec2& b)
{
	m[0 * 2 + 0] = a.x;
	m[0 * 2 + 1] = b.x;
	m[1 * 2 + 0] = a.y;
	m[1 * 2 + 1] = b.y;
}
template <>
__host__ __device__ fmat2::mat2(const float (&array)[4])
{

	for (uint32_t i = 0; i < 9; i++)
		m[i] = array[i];
}
template <>
__host__ __device__ fmat2::mat2(const double& omega)
{
	float cos = cosf(omega);
	float sin = sinf(omega);

	m[0] = cos;
	m[1] = -sin;
	m[2] = sin;
	m[3] = cos;
}

template <>
__host__ __device__ float fmat2::det(void) const;
template <>
__host__ __device__ fmat2 fmat2::inverse(void) const;
template <>
__host__ __device__ float fmat2::sqlength(void) const;
template <>
__host__ __device__ float fmat2::trace(void) const;
template <>
__host__ __device__ fmat2 fmat2::transpose(void) const;
template <>
__host__ __device__ fmat2 fmat2::identity();

__host__ __device__ const fmat2 operator+(const fmat2& mat0, const fmat2& mat1);
__host__ __device__ const fmat2 operator-(const fmat2& mat0, const fmat2& mat1);
__host__ __device__ const fmat2 operator-(const fmat2& mat);
__host__ __device__ const fmat2 operator*(const fmat2& mat, const float& a);
__host__ __device__ const fmat2 operator*(const float& a, const fmat2& mat);
__host__ __device__ const fmat2 operator*(const fmat2& mat0, const fmat2& mat1);
__host__ __device__ const fvec2 operator*(const fmat2& mat, const fvec2& vec);
__host__ __device__ const fmat2 operator/(const fmat2& mat, const float& a);

//////////////////////////////////
//////////////////////////////////
//////////////////////////////////

template <>
__host__ __device__ fmat3::mat3(void)
{
	for (uint32_t i = 0; i < 9; i++)
		m[i] = 0.0;
}

template <>
__host__ __device__ fmat3::mat3(const fvec3& a, const fvec3& b, const fvec3& c)
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

template <>
__host__ __device__ fmat3::mat3(const float (&array)[9])
{
	for (uint32_t i = 0; i < 9; i++)
		m[i] = array[i];
}

template <>
__host__ __device__ float fmat3::det(void) const;
template <>
__host__ __device__ fmat3 fmat3::inverse(void) const;
template <>
__host__ __device__ float fmat3::sqlength() const;
template <>
__host__ __device__ float fmat3::trace() const;
template <>
__host__ __device__ fmat3 fmat3::transpose() const;
template <>
__host__ __device__ fmat3 fmat3::indentity();

__host__ __device__ const fmat3 operator+(const fmat3& mat0, const fmat3& mat1);
__host__ __device__ const fmat3 operator-(const fmat3& mat0, const fmat3& mat1);
__host__ __device__ const fmat3 operator-(const fmat3& mat);
__host__ __device__ const fmat3 operator*(const fmat3& mat, const float& a);
__host__ __device__ const fmat3 operator*(const float& a, const fmat3& mat);
__host__ __device__ const fmat3 operator*(const fmat3& mat0, const fmat3& mat1);
__host__ __device__ const fvec3 operator*(const fmat3& mat, const fvec3& vec);
__host__ __device__ const fmat3 operator/(const fmat3& mat, const float& a);

//////////////////////////////////
//////////////////////////////////
//////////////////////////////////

template <>
__host__ __device__ fmat32::mat32(const fvec3& a, const fvec3& b)
{
	m[0] = a.x;
	m[2] = a.y;
	m[4] = a.z;
	m[1] = b.x;
	m[3] = b.y;
	m[5] = b.z;
}
template <>
__host__ __device__ fmat32::mat32(const float (&a)[6])
{
	for (uint32_t i = 0; i < 6; i++)
		m[i] = a[i];
}

template <>
__host__ __device__ float fmat32::sqlength() const;

__host__ __device__ const fmat32 operator+(const fmat32& mat0, const fmat32& mat1);
__host__ __device__ const fmat32 operator-(const fmat32& mat0, const fmat32& mat1);
__host__ __device__ const fmat32 operator*(const fmat32& mat, const float& a);
__host__ __device__ const fmat32 operator*(const float& a, const fmat32& mat);
__host__ __device__ const fmat32 operator*(const fmat32& mat0, const fmat2& mat1);
__host__ __device__ const fvec3 operator*(const fmat32& mat, const fvec2& vec);

//////////////////////////////////
//////////////////////////////////
//////////////////////////////////

template <>
__host__ __device__ fmat4::mat4(const fvec4& a, const fvec4& b, const fvec4& c, const fvec4& d)
{
	m[4 * 0 + 0] = a.x;
	m[4 * 1 + 0] = a.y;
	m[4 * 2 + 0] = a.z;
	m[4 * 3 + 0] = a.w;

	m[4 * 0 + 1] = b.x;
	m[4 * 1 + 1] = b.y;
	m[4 * 2 + 1] = b.z;
	m[4 * 3 + 1] = b.w;

	m[4 * 0 + 2] = c.x;
	m[4 * 1 + 2] = c.y;
	m[4 * 2 + 2] = c.z;
	m[4 * 3 + 2] = c.w;

	m[4 * 0 + 3] = d.x;
	m[4 * 1 + 3] = d.y;
	m[4 * 2 + 3] = d.z;
	m[4 * 3 + 3] = d.w;
}
template <>
__host__ __device__ fmat4::mat4(const float (&a)[16])
{
	for (uint32_t i = 0; i < 16; i++)
		m[i] = a[i];
}
template <>
__host__ __device__ fmat4::mat4(void)
{
	for (uint32_t i = 0; i < 16; i++)
		m[i] = 0.0;
}
template <>
__host__ __device__ fmat4::mat4(const fmat3& mat, const fvec3& vec)
{
	m[0]  = mat.m[0];
	m[1]  = mat.m[1];
	m[2]  = mat.m[2];
	m[3]  = vec.x;
	m[4]  = mat.m[3];
	m[5]  = mat.m[4];
	m[6]  = mat.m[5];
	m[7]  = vec.y;
	m[8]  = mat.m[6];
	m[9]  = mat.m[7];
	m[10] = mat.m[8];
	m[11] = vec.z;
	m[12] = 0.0;
	m[13] = 0.0;
	m[14] = 0.0;
	m[15] = 1.0;
}

template <>
__host__ __device__ fmat4 fmat4::transpose() const;
template <>
__host__ __device__ fmat4 fmat4::indentity();

__host__ __device__ const fmat4 operator*(const fmat4& mat0, const fmat4& mat1);
__host__ __device__ const fvec4 operator*(const fmat4& mat, const fvec4& vec);

//////////////////////////////////
//////////////////////////////////
//////////////////////////////////

template <>
__host__ __device__ fquaternion::quaternion(
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

template <>
__host__ __device__ fquaternion::quaternion(const fvec3& v, const float& s)
    : x(v.x)
    , y(v.y)
    , z(v.z)
    , w(s)
{
}
template <>
__host__ __device__ fquaternion::quaternion(void)
    : x(0.0)
    , y(0.0)
    , z(0.0)
    , w(1.0)
{
}

template <>
__host__ __device__ fquaternion::quaternion(const fvec3& v);

template <>
__host__ __device__ fvec3 fquaternion::getv() const;
template <>
__host__ __device__ fquaternion fquaternion::conjugate() const;
template <>
__host__ __device__ float fquaternion::dot(const fquaternion& q) const;
template <>
__host__ __device__ float fquaternion::sqlength() const;
template <>
__host__ __device__ float fquaternion::length() const;
template <>
__host__ __device__ fquaternion fquaternion::normalize() const;
template <>
__host__ __device__ fmat3 fquaternion::qtomat() const;
template <>
__host__ __device__ fquaternion fquaternion::inverse() const;

template <>
__host__ __device__ fquaternion fquaternion::slerp(const fquaternion& q0, const fquaternion& q1, const float& t);

__host__ __device__ const fquaternion operator+(const fquaternion& q0, const fquaternion& q1);
__host__ __device__ const fquaternion operator-(const fquaternion& q0, const fquaternion& q1);
__host__ __device__ const fquaternion operator-(const fquaternion& q);
__host__ __device__ const fquaternion operator*(const fquaternion& q0, const fquaternion& q1);
__host__ __device__ const fquaternion operator*(const fquaternion& q, const float& a);
__host__ __device__ const fquaternion operator*(const float& a, const fquaternion& q);
__host__ __device__ const fquaternion operator/(const fquaternion& q, const float& a);

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////↓  実装///////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

template <>
__host__ __device__ float fvec2::dot(const fvec2& v) const
{
	return x * v.x + y * v.y;
}
template <>
__host__ __device__ float fvec2::cross(const fvec2& v) const
{
	return x * v.y - y * v.x;
}
template <>
__host__ __device__ float fvec2::sqlength() const
{
	return x * x + y * y;
}
template <>
__host__ __device__ float fvec2::length() const
{
	return sqrtf(this->sqlength());
}
template <>
__host__ __device__ fvec2 fvec2::normalize() const
{
	float length = this->length();
	if (length < 0.0000001)
		return fvec2();
	return (*this) / length;
}

__host__ __device__ const fvec2 operator+(const fvec2& a, const fvec2& b)
{
	return fvec2(
	    a.x + b.x,
	    a.y + b.y);
}
__host__ __device__ const fvec2 operator-(const fvec2& a, const fvec2& b)
{
	return fvec2(
	    a.x - b.x,
	    a.y - b.y);
}
__host__ __device__ const fvec2 operator-(const fvec2& v)
{
	return fvec2(
	    -v.x,
	    -v.y);
}
__host__ __device__ const fvec2 operator*(const float& a, const fvec2& v)
{
	return fvec2(
	    a * v.x,
	    a * v.y);
}
__host__ __device__ const fvec2 operator*(const fvec2& v, const float& a)
{
	return fvec2(
	    a * v.x,
	    a * v.y);
}
__host__ __device__ const fvec2 operator/(const fvec2& v, const float& a)
{
	return fvec2(
	    v.x / a,
	    v.y / a);
}

//////////////////////////////////
//////////////////////////////////
//////////////////////////////////

template <>
__host__ __device__ float fvec3::dot(const fvec3& v) const
{
	return x * v.x + y * v.y + z * v.z;
}
template <>
__host__ __device__ fvec3 fvec3::cross(const fvec3& a) const
{
	return fvec3(y * a.z - z * a.y, z * a.x - x * a.z, x * a.y - y * a.x);
}
template <>
__host__ __device__ float fvec3::sqlength() const
{
	return x * x + y * y + z * z;
}
template <>
__host__ __device__ float fvec3::length() const
{
	return sqrtf(this->sqlength());
}
template <>
__host__ __device__ fvec3 fvec3::normalize() const
{
	float length = this->length();
	if (length < 0.0000001)
		return fvec3();
	return (*this) / length;
}
template <>
__host__ __device__ fmat3 fvec3::tensorproduct(const fvec3& v) const
{
	fmat3 hoge;
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
template <>
__host__ __device__ fmat3 fvec3::skew() const
{
	fmat3 hoge;
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
template <>
__host__ __device__ fmat3 fvec3::rotation() const
{
	float omega  = this->length();
	float omega2 = this->sqlength();
	if (omega2 < 0.0000001)
		return fmat3::indentity();
	float cos = cosf(omega);
	float sin = sinf(omega);

	return cos * fmat3::indentity() + ((1 - cos) / omega2) * this->tensorproduct(*this) + (sin / omega) * this->skew();
}
template <>
__host__ __device__ float fvec3::STP(const fvec3& a, const fvec3& b, const fvec3& c)
{
	return a.dot(b.cross(c));
}

__host__ __device__ const fvec3 operator+(const fvec3& a, const fvec3& b)
{
	return fvec3(
	    a.x + b.x,
	    a.y + b.y,
	    a.z + b.z);
}
__host__ __device__ const fvec3 operator-(const fvec3& a, const fvec3& b)
{
	return fvec3(
	    a.x - b.x,
	    a.y - b.y,
	    a.z - b.z);
}
__host__ __device__ const fvec3 operator-(const fvec3& v)
{
	return fvec3(
	    -v.x,
	    -v.y,
	    -v.z);
}
__host__ __device__ const fvec3 operator*(const float& a, const fvec3& v)
{
	return fvec3(
	    a * v.x,
	    a * v.y,
	    a * v.z);
}
__host__ __device__ const fvec3 operator*(const fvec3& v, const float& a)
{
	return fvec3(
	    a * v.x,
	    a * v.y,
	    a * v.z);
}
__host__ __device__ const fvec3 operator/(const fvec3& v, const float& a)
{
	return fvec3(
	    v.x / a,
	    v.y / a,
	    v.z / a);
}

//////////////////////////////////
//////////////////////////////////
//////////////////////////////////

template <>
__host__ __device__ float fvec4::dot(const fvec4& v) const
{
	return x * v.x + y * v.y + z * v.z + w * v.w;
}
template <>
__host__ __device__ float fvec4::sqlength(void) const
{
	return x * x + y * y + z * z + w * w;
}

__host__ __device__ const fvec4 operator+(const fvec4& a, const fvec4& b)
{
	return fvec4(
	    a.x + b.x,
	    a.y + b.y,
	    a.z + b.z,
	    a.w + b.w);
}
__host__ __device__ const fvec4 operator-(const fvec4& a, const fvec4& b)
{
	return fvec4(
	    a.x - b.x,
	    a.y - b.y,
	    a.z - b.z,
	    a.w - b.w);
}
__host__ __device__ const fvec4 operator*(const float& a, const fvec4& v)
{
	return fvec4(
	    a * v.x,
	    a * v.y,
	    a * v.z,
	    a * v.w);
}
__host__ __device__ const fvec4 operator*(const fvec4& v, const float& a)
{
	return fvec4(
	    a * v.x,
	    a * v.y,
	    a * v.z,
	    a * v.w);
}

//////////////////////////////////
//////////////////////////////////
//////////////////////////////////

template <>
__host__ __device__ float fmat2::det(void) const
{
	return m[0] * m[3] - m[1] * m[2];
}
template <>
__host__ __device__ fmat2 fmat2::inverse(void) const
{
	float det     = this->det();
	float data[4] = { m[3] / det, -m[1] / det, -m[2] / det, m[0] / det };

	return fmat2(data);
}
template <>
__host__ __device__ float fmat2::sqlength(void) const
{
	return m[0] * m[0] + m[1] * m[1] + m[2] * m[2] + m[3] * m[3];
}
template <>
__host__ __device__ float fmat2::trace(void) const
{
	return m[0] + m[3];
}

template <>
__host__ __device__ fmat2 fmat2::transpose(void) const
{
	float data[4];
	data[0] = m[0];
	data[1] = m[2];
	data[2] = m[1];
	data[3] = m[3];
	return fmat2(data);
}
template <>
__host__ __device__ fmat2 fmat2::identity()
{
	fmat2 temp;
	temp.m[0] = 1.0;
	temp.m[3] = 1.0;
	return temp;
}

__host__ __device__ const fmat2 operator+(const fmat2& mat0, const fmat2& mat1)
{
	float p[4];
	for (uint32_t i = 0; i < 4; i++)
		p[i] = mat0.m[i] + mat1.m[i];
	return fmat2(p);
}
__host__ __device__ const fmat2 operator-(const fmat2& mat0, const fmat2& mat1)
{
	float p[4];
	for (uint32_t i = 0; i < 4; i++)
		p[i] = mat0.m[i] - mat1.m[i];
	return fmat2(p);
}
__host__ __device__ const fmat2 operator-(const fmat2& mat)
{
	float p[4];
	for (uint32_t i = 0; i < 4; i++)
		p[i] = -mat.m[i];
	return fmat2(p);
}
__host__ __device__ const fmat2 operator*(const fmat2& mat, const float& a)
{
	float p[4];
	for (uint32_t i = 0; i < 4; i++)
		p[i] = a * mat.m[i];
	return fmat2(p);
}
__host__ __device__ const fmat2 operator*(const float& a, const fmat2& mat)
{
	float p[4];
	for (uint32_t i = 0; i < 4; i++)
		p[i] = a * mat.m[i];
	return fmat2(p);
}
__host__ __device__ const fmat2 operator*(const fmat2& mat0, const fmat2& mat1)
{
	float p[4];
	for (uint32_t i = 0; i < 2; i++) {
		for (uint32_t j = 0; j < 2; j++) {
			p[i * 2 + j] = mat0.m[i * 2 + 0] * mat1.m[0 * 2 + j] + mat0.m[i * 2 + 1] * mat1.m[1 * 2 + j];
		}
	}
	return fmat2(p);
}
__host__ __device__ const fvec2 operator*(const fmat2& mat, const fvec2& vec)
{
	return fvec2(mat.m[0] * vec.x + mat.m[1] * vec.y, mat.m[2] * vec.x + mat.m[3] * vec.y);
}
__host__ __device__ const fmat2 operator/(const fmat2& mat, const float& a)
{
	float p[4];
	for (uint32_t i = 0; i < 4; i++)
		p[i] = mat.m[i] / a;
	return fmat2(p);
}

//////////////////////////////////
//////////////////////////////////
//////////////////////////////////

template <>
__host__ __device__ float fmat3::det(void) const
{

	float temp[3];
	temp[0] = m[1] * m[5] - m[2] * m[4];
	temp[1] = m[2] * m[3] - m[0] * m[5];
	temp[2] = m[0] * m[4] - m[1] * m[3];

	return m[6] * temp[0] + m[7] * temp[1] + m[8] * temp[2];
}
template <>
__host__ __device__ fmat3 fmat3::inverse(void) const
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

	return fmat3(data);
}
template <>
__host__ __device__ float fmat3::sqlength() const
{
	float hoge = 0.0;
	for (uint32_t i = 0; i < 9; i++)
		hoge += m[i];
	return hoge;
}
template <>
__host__ __device__ float fmat3::trace() const
{
	return m[0] + m[4] + m[8];
}
template <>
__host__ __device__ fmat3 fmat3::transpose() const
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

	return fmat3(data);
}

template <>
__host__ __device__ fmat3 fmat3::indentity()
{
	fmat3 temp;
	temp.m[0] = 1.0;
	temp.m[4] = 1.0;
	temp.m[8] = 1.0;
	return temp;
}

__host__ __device__ const fmat3 operator+(const fmat3& mat0, const fmat3& mat1)
{
	float p[9];
	for (uint32_t i = 0; i < 9; i++)
		p[i] = mat0.m[i] + mat1.m[i];
	return fmat3(p);
}
__host__ __device__ const fmat3 operator-(const fmat3& mat0, const fmat3& mat1)
{
	float p[9];
	for (uint32_t i = 0; i < 9; i++)
		p[i] = mat0.m[i] - mat1.m[i];
	return fmat3(p);
}
__host__ __device__ const fmat3 operator-(const fmat3& mat)
{
	float p[9];
	for (uint32_t i = 0; i < 9; i++)
		p[i] = -mat.m[i];
	return fmat3(p);
}
__host__ __device__ const fmat3 operator*(const fmat3& mat, const float& a)
{
	float p[9];
	for (uint32_t i = 0; i < 9; i++)
		p[i] = a * mat.m[i];
	return fmat3(p);
}
__host__ __device__ const fmat3 operator*(const float& a, const fmat3& mat)
{
	float p[9];
	for (uint32_t i = 0; i < 9; i++)
		p[i] = a * mat.m[i];
	return fmat3(p);
}
__host__ __device__ const fmat3 operator*(const fmat3& mat0, const fmat3& mat1)
{
	float p[9];
	for (uint32_t i = 0; i < 3; i++) {
		for (uint32_t j = 0; j < 3; j++) {
			p[i * 3 + j] = mat0.m[i * 3 + 0] * mat1.m[0 * 3 + j] + mat0.m[i * 3 + 1] * mat1.m[1 * 3 + j] + mat0.m[i * 3 + 2] * mat1.m[2 * 3 + j];
		}
	}
	return fmat3(p);
}
__host__ __device__ const fvec3 operator*(const fmat3& mat, const fvec3& a)
{
	return fvec3(mat.m[0 + 0] * a.x + mat.m[0 + 1] * a.y + mat.m[0 + 2] * a.z,
	    mat.m[3 + 0] * a.x + mat.m[3 + 1] * a.y + mat.m[3 + 2] * a.z,
	    mat.m[6 + 0] * a.x + mat.m[6 + 1] * a.y + mat.m[6 + 2] * a.z);
}
__host__ __device__ const fmat3 operator/(const fmat3& mat, const float& a)
{
	float p[9];
	for (uint32_t i = 0; i < 9; i++)
		p[i] = mat.m[i] / a;
	return fmat3(p);
}

//////////////////////////////////
//////////////////////////////////
//////////////////////////////////

template <>
__host__ __device__ float fmat32::sqlength() const
{
	float hoge = 0.0;
	for (uint32_t i = 0; i < 6; i++) {
		hoge += m[i] * m[i];
	}
	return hoge;
}

__host__ __device__ const fmat32 operator+(const fmat32& mat0, const fmat32& mat1)
{
	float p[6];
	for (uint32_t i = 0; i < 6; i++)
		p[i] = mat0.m[i] + mat1.m[i];
	return fmat32(p);
}
__host__ __device__ const fmat32 operator-(const fmat32& mat0, const fmat32& mat1)
{
	float p[6];
	for (uint32_t i = 0; i < 6; i++)
		p[i] = mat0.m[i] - mat1.m[i];
	return fmat32(p);
}
__host__ __device__ const fmat32 operator*(const fmat32& mat, const float& a)
{
	float p[6];
	for (uint32_t i = 0; i < 6; i++)
		p[i] = a * mat.m[i];
	return fmat32(p);
}
__host__ __device__ const fmat32 operator*(const float& a, const fmat32& mat)
{
	float p[6];
	for (uint32_t i = 0; i < 6; i++)
		p[i] = a * mat.m[i];
	return fmat32(p);
}
__host__ __device__ const fmat32 operator*(const fmat32& mat0, const fmat2& mat1)
{
	float p[6];
	for (uint32_t i = 0; i < 3; i++) {
		for (uint32_t j = 0; j < 2; j++) {
			p[i * 2 + j] = mat0.m[2 * i + 0] * mat1.m[2 * 0 + j] + mat0.m[2 * i + 1] * mat1.m[2 * 1 + j];
		}
	}

	return fmat32(p);
}
__host__ __device__ const fvec3 operator*(const fmat32& mat, const fvec2& vec)
{
	return fvec3(
	    mat.m[2 * 0 + 0] * vec.x + mat.m[2 * 0 + 1] * vec.y,
	    mat.m[2 * 1 + 0] * vec.x + mat.m[2 * 1 + 1] * vec.y,
	    mat.m[2 * 2 + 0] * vec.x + mat.m[2 * 2 + 1] * vec.y);
}

//////////////////////////////////
//////////////////////////////////
//////////////////////////////////

template <>
__host__ __device__ fmat4 fmat4::transpose() const
{
	float data[16];

	data[0]	 = m[0];
	data[1]	 = m[4];
	data[2]	 = m[8];
	data[3]	 = m[12];
	data[4]	 = m[1];
	data[5]	 = m[5];
	data[6]	 = m[9];
	data[7]	 = m[13];
	data[8]	 = m[2];
	data[9]	 = m[6];
	data[10] = m[10];
	data[11] = m[14];
	data[12] = m[3];
	data[13] = m[7];
	data[14] = m[11];
	data[15] = m[15];

	return fmat4(data);
}

template <>
__host__ __device__ fmat4 fmat4::indentity()
{
	fmat4 temp;
	temp.m[0]  = 1.0;
	temp.m[5]  = 1.0;
	temp.m[10] = 1.0;
	temp.m[15] = 1.0;
	return temp;
}

__host__ __device__ const fmat4 operator*(const fmat4& mat0, const fmat4& mat1)
{
	float p[16];
	for (uint32_t i = 0; i < 4; i++) {
		for (uint32_t j = 0; j < 4; j++) {
			p[i * 4 + j] = mat0.m[i * 4 + 0] * mat1.m[0 * 4 + j] + mat0.m[i * 4 + 1] * mat1.m[1 * 4 + j] + mat0.m[i * 4 + 2] * mat1.m[2 * 4 + j] + mat0.m[i * 4 + 3] * mat1.m[3 * 4 + j];
		}
	}

	return fmat4(p);
}

__host__ __device__ const fvec4 operator*(const fmat4& mat, const fvec4& vec)
{
	return fvec4(
	    mat.m[4 * 0 + 0] * vec.x + mat.m[4 * 0 + 1] * vec.y + mat.m[4 * 0 + 2] * vec.z + mat.m[4 * 0 + 3] * vec.w,
	    mat.m[4 * 1 + 0] * vec.x + mat.m[4 * 1 + 1] * vec.y + mat.m[4 * 1 + 2] * vec.z + mat.m[4 * 1 + 3] * vec.w,
	    mat.m[4 * 2 + 0] * vec.x + mat.m[4 * 2 + 1] * vec.y + mat.m[4 * 2 + 2] * vec.z + mat.m[4 * 2 + 3] * vec.w,
	    mat.m[4 * 3 + 0] * vec.x + mat.m[4 * 3 + 1] * vec.y + mat.m[4 * 3 + 2] * vec.z + mat.m[4 * 3 + 3] * vec.w);
}

//////////////////////////////////
//////////////////////////////////
//////////////////////////////////

template <>
__host__ __device__ fquaternion::quaternion(const fvec3& v)
{
	float Omega = v.length();
	fvec3 l	    = v.normalize();

	float sin = sinf(Omega / 2.0);
	x	  = sin * l.x;
	y	  = sin * l.y;
	z	  = sin * l.z;
	w	  = cosf(Omega / 2.0);
}

template <>
__host__ __device__ fvec3 fquaternion::getv() const
{
	return fvec3(x, y, z);
}
template <>
__host__ __device__ fquaternion fquaternion::conjugate() const
{
	return fquaternion(-1.0f * fvec3(x, y, z), w);
}
template <>
__host__ __device__ float fquaternion::dot(const fquaternion& q) const
{
	return x * q.x + y * q.y + z * q.z + w * q.w;
}
template <>
__host__ __device__ float fquaternion::sqlength() const
{
	return x * x + y * y + z * z + w * w;
}
template <>
__host__ __device__ float fquaternion::length() const
{
	return sqrtf(x * x + y * y + z * z + w * w);
}
template <>
__host__ __device__ fquaternion fquaternion::normalize() const
{
	float length = this->length();
	if (length < 0.0000001)
		return fquaternion();
	return (*this) / length;
}
template <>
__host__ __device__ fmat3 fquaternion::qtomat() const
{
	fvec3 v = this->getv();
	float w = this->w;

	return (w * w - v.sqlength()) * fmat3::indentity() + 2.0f * w * v.skew() + 2.0f * v.tensorproduct(v);
}
template <>
__host__ __device__ fquaternion fquaternion::inverse() const
{
	return fquaternion(-x, -y, -z, w);
}

template <>
__host__ __device__ fquaternion fquaternion::slerp(const fquaternion& q0, const fquaternion& q1, const float& t)
{
	if (q0.dot(q1) < 0.0) {
		double x = 0.5 * acosf(-q0.dot(q1));
		return (sinf(t * x) / sinf(x)) * -q1 + (sinf(x - t * x) / sinf(x)) * q0;
	} else {
		double x = 0.5 * acosf(q0.dot(q1));
		return (sinf(t * x) / sinf(x)) * q1 + (sinf(x - t * x) / sinf(x)) * q0;
	}
}

__host__ __device__ const fquaternion operator+(const fquaternion& q0, const fquaternion& q1)
{
	return fquaternion(
	    q0.x + q1.x,
	    q0.y + q1.y,
	    q0.z + q1.z,
	    q0.w + q1.w);
}
__host__ __device__ const fquaternion operator-(const fquaternion& q0, const fquaternion& q1)
{
	return fquaternion(
	    q0.x - q1.x,
	    q0.y - q1.y,
	    q0.z - q1.z,
	    q0.w - q1.w);
}
__host__ __device__ const fquaternion operator-(const fquaternion& q)
{
	return fquaternion(
	    -q.x,
	    -q.y,
	    -q.z,
	    -q.w);
}
__host__ __device__ const fquaternion operator*(const fquaternion& q0, const fquaternion& q1)
{
	float w = q0.w * q1.w - (q0.x * q1.x + q0.y * q1.y + q0.z * q1.z);
	fvec3 v = q0.w * q1.getv() + q1.w * q0.getv() + q0.getv().cross(q1.getv());

	return fquaternion(v, w);
}
__host__ __device__ const fquaternion operator*(const fquaternion& q, const float& a)
{
	return fquaternion(
	    a * q.x,
	    a * q.y,
	    a * q.z,
	    a * q.w);
}
__host__ __device__ const fquaternion operator*(const float& a, const fquaternion& q)
{
	return fquaternion(
	    a * q.x,
	    a * q.y,
	    a * q.z,
	    a * q.w);
}
__host__ __device__ const fquaternion operator/(const fquaternion& q, const float& a)
{
	return fquaternion(
	    q.x / a,
	    q.y / a,
	    q.z / a,
	    q.w / a);
}

#endif
