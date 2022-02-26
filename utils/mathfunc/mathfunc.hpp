#pragma once

#include <iostream>
#include <cstdint>

//TODO 宣言と定義を分離する(cudamathと同じ形式にする)

template <class T>
class vec2;
template <class T>
class vec3;
template <class T>
class vec4;
template <class T>
class mat2;
template <class T>
class mat3;
template <class T>
class mat32;
template <class T>
class mat4;
template <class T>
class quaternion;

////vec2

template <class T>
class vec2 {
    public:
	T x, y;

	//constructor
	inline vec2(const T& x, const T& y)
	    : x(x)
	    , y(y)
	{
	}
	inline vec2(const T& value)
	    : x(value)
	    , y(value)
	{
	}
	inline vec2(const T (&a)[2])
	{
		x = a[0];
		y = a[1];
	}
	inline vec2(void)
	    : x(0.0)
	    , y(0.0)
	{
	}

	//menber function
	inline T dot(const vec2<T>& a) const
	{
		return x * a.x + y * a.y;
	}
	inline T cross(const vec2<T>& a) const
	{
		return x * a.y - y * a.x;
	}

	inline T sqlength() const
	{
		return x * x + y * y;
	}

	T length() const;

	vec2<T> normalize() const;

	vec2<T> rot() const;

	template <class U>
	operator vec2<U>() const
	{
		U cx = U(x);
		U cy = U(y);
		return vec2<U>(cx, cy);
	}

	//static function
	inline static T dot(const vec2<T>& a, const vec2<T>& b)
	{
		return a.dot(b);
	}
	inline static T cross(const vec2<T>& a, const vec2<T>& b)
	{
		return a.cross(b);
	}
};

//io operator

template <class T>
std::ostream& operator<<(std::ostream& os, const vec2<T>& vec);

//operator

template <class T>
inline vec2<T> operator+(const vec2<T>& v, const vec2<T>& a)
{
	return vec2<T>(v.x + a.x, v.y + a.y);
}
template <class T>
inline vec2<T> operator-(const vec2<T>& v, const vec2<T>& a)
{
	return vec2<T>(v.x - a.x, v.y - a.y);
}

template <class T>
inline vec2<T> operator-(const vec2<T>& v)
{
	return vec2<T>(-v.x, -v.y);
}

template <class T>
inline vec2<T> operator*(const vec2<T>& v, const vec2<T>& a)
{
	return vec2<T>(v.x * a.x, v.y * a.y);
}

template <class T, class U>
inline vec2<T> operator*(const vec2<T>& v, const U& a)
{
	return vec2<T>(v.x * a, v.y * a);
}
template <class T, class U>
inline vec2<T> operator*(const U& a, const vec2<T>& v)
{
	return vec2<T>(v.x * a, v.y * a);
}

template <class T>
inline vec2<T> operator/(const vec2<T>& v, const vec2<T>& a)
{
	return vec2<T>(v.x / a.x, v.y / a.y);
}

template <class T, class U>
inline vec2<T> operator/(const vec2<T>& v, const U& a)
{
	return vec2<T>(v.x / a, v.y / a);
}
/*
template <class T, class U>
inline vec2<T> operator/(const U& a, const vec2<T>& v)
{
	return vec2<T>(v.x / a, v.y / a);
}
*/

//alias

using fvec2 = vec2<float>;
using dvec2 = vec2<double>;

////vec3

template <class T>
class vec3 {
    public:
	T x, y, z;

	//constructor
	inline vec3(const T& x, const T& y, const T& z)
	    : x(x)
	    , y(y)
	    , z(z)
	{
	}
	inline vec3(const T& value)
	    : x(value)
	    , y(value)
	    , z(value)
	{
	}
	inline vec3(void)
	    : x(0.0)
	    , y(0.0)
	    , z(0.0)
	{
	}
	inline vec3(const T (&a)[3])
	{
		x = a[0];
		y = a[1];
		z = a[2];
	}

	//menber function
	inline T dot(const vec3<T>& a) const
	{
		return x * a.x + y * a.y + z * a.z;
	}
	inline vec3<T> cross(const vec3<T>& a) const
	{
		return vec3(y * a.z - z * a.y, z * a.x - x * a.z, x * a.y - y * a.x);
	}

	inline T sqlength() const
	{
		return x * x + y * y + z * z;
	}

	T length() const;

	vec3<T> normalize() const;

	mat3<T> tensorproduct(const vec3<T>& v) const;

	mat3<T> skew() const;

	mat3<T> rotation() const;

	//cast operator

	template <class U>
	operator vec3<U>() const
	{
		U cx = U(x);
		U cy = U(y);
		U cz = U(z);
		return vec3<U>(cx, cy, cz);
	}

	//static function
	inline static T dot(const vec3<T>& a, const vec3<T>& b)
	{
		return a.x * b.x + a.y * b.y + a.z * b.z;
	}
	inline static vec3<T> cross(const vec3<T>& a, const vec3<T>& b)
	{
		return vec3<T>(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
	}
	inline static T STP(const vec3<T>& a, const vec3<T>& b, const vec3<T>& c)
	{
		return vec3::dot(a, vec3::cross(b, c));
	}
};

//io operator

template <class T>
std::ostream& operator<<(std::ostream& os, const vec3<T>& vec);

//operator

template <class T>
inline vec3<T> operator+(const vec3<T>& v, const vec3<T>& a)
{
	return vec3<T>(v.x + a.x, v.y + a.y, v.z + a.z);
}
template <class T>
inline vec3<T> operator-(const vec3<T>& v, const vec3<T>& a)
{
	return vec3<T>(v.x - a.x, v.y - a.y, v.z - a.z);
}

template <class T>
inline vec3<T> operator-(const vec3<T>& v)
{
	return vec3<T>(-v.x, -v.y, -v.z);
}

template <class T>
inline vec3<T> operator*(const vec3<T>& v, const vec3<T>& a)
{
	return vec3<T>(v.x * a.x, v.y * a.y, v.z * a.z);
}
template <class T, class U>
inline vec3<T> operator*(const vec3<T>& v, const U& a)
{
	return vec3<T>(v.x * a, v.y * a, v.z * a);
}
template <class T, class U>
inline vec3<T> operator*(const U& a, const vec3<T>& v)
{
	return vec3<T>(v.x * a, v.y * a, v.z * a);
}

template <class T>
inline vec3<T> operator/(const vec3<T>& v, const vec3<T>& a)
{
	return vec3<T>(v.x / a.x, v.y / a.y, v.z / a.z);
}
template <class T, class U>
inline vec3<T> operator/(const vec3<T>& v, const U& a)
{
	return vec3<T>(v.x / a, v.y / a, v.z / a);
}
//template <class T, class U>
//inline vec3<T> operator/(const U& a, const vec3<T>& v)
//{
//	return vec3<T>(v.x / a, v.y / a, v.z / a);
//}

//alias

using fvec3 = vec3<float>;
using dvec3 = vec3<float>;

////vec4

template <class T>
class vec4 {
    public:
	T x, y, z, w;

	//constructor
	inline vec4(
	    const T& x,
	    const T& y,
	    const T& z,
	    const T& w)
	    : x(x)
	    , y(y)
	    , z(z)
	    , w(w)
	{
	}

	inline vec4(void)
	    : x(0.0)
	    , y(0.0)
	    , z(0.0)
	    , w(0.0)
	{
	}

	inline vec4(
	    const vec3<T>& a)
	    : x(a.x)
	    , y(a.y)
	    , z(a.z)
	    , w(1.0)
	{
	}

	//member func
	inline T dot(const vec4<T>& a) const
	{
		return x * a.x + y * a.y + z * a.z + w * a.w;
	}

	template <class U>
	operator vec4<U>() const
	{
		U cx = U(x);
		U cy = U(y);
		U cz = U(z);
		U cw = U(w);

		return vec4<U>(cx, cy, cz, cw);
	}

	inline T sqlength() const
	{
		return x * x + y * y + z * z + w * w;
	}
};

//io operator

template <class T>
std::ostream& operator<<(std::ostream& os, const vec4<T>& vec);

//operator

template <class T>
inline const vec4<T> operator+(const vec4<T>& v, const vec4<T>& a)
{
	return vec4<T>(v.x + a.x, v.y + a.y, v.z + a.z, v.w + a.w);
}

template <class T>
inline const vec4<T> operator-(const vec4<T>& v, const vec4<T>& a)
{
	return vec4<T>(v.x - a.x, v.y - a.y, v.z - a.z, v.w - a.w);
}

template <class T, class U>
inline const vec4<T> operator*(const vec4<T>& v, const U& a)
{
	return vec4<T>(v.x * a, v.y * a, v.z * a, v.w * a);
}
template <class T, class U>
inline const vec4<T> operator*(const U& a, const vec4<T>& v)
{
	return vec4<T>(v.x * a, v.y * a, v.z * a, v.w * a);
}

//alias

using fvec4 = vec4<float>;
using dvec4 = vec4<double>;

////mat2

template <class T>
class mat2 {
    public:
	T m[4];
	// m[0],m[1]
	// m[2],m[3]

	//constructor
	inline mat2(const vec2<T>& a, const vec2<T>& b)
	{
		m[0 * 2 + 0] = a.x;
		m[0 * 2 + 1] = b.x;
		m[1 * 2 + 0] = a.y;
		m[1 * 2 + 1] = b.y;
	}

	inline mat2(const T (&a)[4])
	{
		for (uint32_t i = 0; i < 4; i++)
			m[i] = a[i];
	}

	inline mat2(void)
	{
		for (uint32_t i = 0; i < 4; i++)
			m[i] = 0.0;
	}

	mat2(const double& omega);

	//menber func
	inline T det() const
	{
		return m[0] * m[3] - m[1] * m[2];
	}

	inline mat2<T> inverse() const
	{
		T det	  = this->det();
		T data[4] = { m[3] / det, -m[1] / det, -m[2] / det, m[0] / det };

		return mat2<T>(data);
	}
	inline T sqlength() const
	{
		return m[0] * m[0] + m[1] * m[1] + m[2] * m[2] + m[3] * m[3];
	}

	inline T trace() const
	{
		return m[0] + m[3];
	}

	inline mat2<T> transpose() const
	{
		T data[4];
		data[0] = m[0];
		data[1] = m[2];
		data[2] = m[1];
		data[3] = m[3];
		return mat2<T>(data);
	}

	//cast operator

	template <class U>
	operator mat2<U>() const
	{
		U cm[4];
		for (uint32_t i = 0; i < 4; i++)
			cm[i] = U(m[i]);

		return mat2<U>(cm);
	}

	//static func
	inline static mat2<T> identity()
	{
		mat2<T> temp = mat2<T>();
		temp.m[0]    = 1.0;
		temp.m[3]    = 1.0;
		return temp;
	}
};

//io operator

template <class T>
std::ostream& operator<<(std::ostream& os, const mat2<T>& mat);

//operator

template <class T>
const inline vec2<T> operator*(const mat2<T>& mat, const vec2<T>& a)
{
	return vec2<T>(mat.m[0] * a.x + mat.m[1] * a.y, mat.m[2] * a.x + mat.m[3] * a.y);
}

template <class T>
const inline mat2<T> operator*(const mat2<T>& mat, const mat2<T>& a)
{
	T p[4];
	for (uint32_t i = 0; i < 2; i++) {
		for (uint32_t j = 0; j < 2; j++) {
			p[i * 2 + j] = mat.m[i * 2 + 0] * a.m[0 * 2 + j] + mat.m[i * 2 + 1] * a.m[1 * 2 + j];
		}
	}

	return mat2<T>(p);
}

template <class T, class U>
inline mat2<T> operator*(const mat2<T>& v, const U& a)
{
	T data[4];
	data[0] = v.m[0] * a;
	data[1] = v.m[1] * a;
	data[2] = v.m[2] * a;
	data[3] = v.m[3] * a;
	return mat2<T>(data);
}
template <class T, class U>
inline mat2<T> operator*(const U& a, const mat2<T>& v)
{
	T data[4];
	data[0] = v.m[0] * a;
	data[1] = v.m[1] * a;
	data[2] = v.m[2] * a;
	data[3] = v.m[3] * a;
	return mat2<T>(data);
}

template <class T>
inline mat2<T> operator+(const mat2<T>& v, const mat2<T>& a)
{
	T data[4];
	data[0] = v.m[0] + a.m[0];
	data[1] = v.m[1] + a.m[1];
	data[2] = v.m[2] + a.m[2];
	data[3] = v.m[3] + a.m[3];
	return mat2<T>(data);
}
template <class T>
inline mat2<T> operator-(const mat2<T>& v, const mat2<T>& a)
{
	T data[4];
	data[0] = v.m[0] - a.m[0];
	data[1] = v.m[1] - a.m[1];
	data[2] = v.m[2] - a.m[2];
	data[3] = v.m[3] - a.m[3];
	return mat2<T>(data);
}

// alias

using fmat2 = mat2<float>;
using dmat2 = mat2<double>;

////mat3

template <class T>
class mat3 {
    public:
	T m[9];
	// m[0],m[1],m[2]
	// m[3],m[4],m[5]
	// m[6],m[7],m[8]

	//constructor
	inline mat3(const vec3<T>& a, const vec3<T>& b, const vec3<T>& c)
	{
		m[0 * 3 + 0] = a.x;
		m[0 * 3 + 1] = b.x;
		m[0 * 3 + 2] = c.x;

		m[1 * 3 + 0] = a.y;
		m[1 * 3 + 1] = b.y;
		m[1 * 3 + 2] = c.y;

		m[2 * 3 + 0] = a.z;
		m[2 * 3 + 1] = b.z;
		m[2 * 3 + 2] = c.z;
	}

	inline mat3(const T (&a)[9])
	{
		for (uint32_t i = 0; i < 9; i++)
			m[i] = a[i];
	}

	inline mat3(void)
	{
		for (uint32_t i = 0; i < 9; i++)
			m[i] = 0.0;
	}

	//menber func
	inline T det() const
	{
		T temp[3];
		temp[0] = m[1] * m[5] - m[2] * m[4];
		temp[1] = m[2] * m[3] - m[0] * m[5];
		temp[2] = m[0] * m[4] - m[1] * m[3];

		return m[6] * temp[0] + m[7] * temp[1] + m[8] * temp[2];
	}

	inline mat3<T> inverse() const
	{
		T det = this->det();
		T data[9];

		data[0] = (m[4] * m[8] - m[5] * m[7]) / det;
		data[1] = -(m[1] * m[8] - m[2] * m[7]) / det;
		data[2] = (m[1] * m[5] - m[2] * m[4]) / det;
		data[3] = -(m[3] * m[8] - m[5] * m[6]) / det;
		data[4] = (m[0] * m[8] - m[2] * m[6]) / det;
		data[5] = -(m[0] * m[5] - m[2] * m[3]) / det;
		data[6] = (m[3] * m[7] - m[4] * m[6]) / det;
		data[7] = -(m[0] * m[7] - m[1] * m[6]) / det;
		data[8] = (m[0] * m[4] - m[1] * m[3]) / det;

		return mat3<T>(data);
	}

	inline T sqlength() const
	{
		T hoge = 0.0;
		for (uint32_t i = 0; i < 9; i++)
			hoge += m[i] * m[i];

		return hoge;
	}

	inline T trace() const
	{
		return m[0] + m[4] + m[8];
	}

	inline mat3<T> transpose() const
	{
		T data[9];

		data[0] = m[0];
		data[1] = m[3];
		data[2] = m[6];
		data[3] = m[1];
		data[4] = m[4];
		data[5] = m[7];
		data[6] = m[2];
		data[7] = m[5];
		data[8] = m[8];

		return mat3<T>(data);
	}

	vec3<T> mattorot() const;

	//cast operator

	template <class U>
	operator mat3<U>() const
	{
		U cm[9];
		for (uint32_t i = 0; i < 9; i++)
			cm[i] = U(m[i]);

		return mat3<U>(cm);
	}

	//static func
	inline static mat3<T> indentity()
	{
		mat3<T> temp = mat3<T>();
		temp.m[0]    = 1.0;
		temp.m[4]    = 1.0;
		temp.m[8]    = 1.0;
		return temp;
	}
};

//io operator

template <class T>
std::ostream& operator<<(std::ostream& os, const mat3<T>& mat);

//operator

template <class T>
const inline mat3<T> operator+(const mat3<T>& mat0, const mat3<T>& mat1)
{
	T p[9];
	for (uint32_t i = 0; i < 9; i++)
		p[i] = mat0.m[i] + mat1.m[i];
	return mat3<T>(p);
}

template <class T>
const inline mat3<T> operator-(const mat3<T>& mat0, const mat3<T>& mat1)
{
	T p[9];
	for (uint32_t i = 0; i < 9; i++)
		p[i] = mat0.m[i] - mat1.m[i];
	return mat3<T>(p);
}

template <class T>
const inline mat3<T> operator-(const mat3<T>& mat)
{
	T p[9];
	for (uint32_t i = 0; i < 9; i++)
		p[i] = -mat.m[i];
	return mat3<T>(p);
}

template <class T, class U>
const inline mat3<T> operator*(const mat3<T>& mat, const U& a)
{
	T p[9];
	for (uint32_t i = 0; i < 9; i++)
		p[i] = mat.m[i] * a;
	return mat3<T>(p);
}

template <class T, class U>
const inline mat3<T> operator*(const U& a, const mat3<T>& mat)
{
	T p[9];
	for (uint32_t i = 0; i < 9; i++)
		p[i] = mat.m[i] * a;
	return mat3<T>(p);
}

template <class T>
const inline vec3<T> operator*(const mat3<T>& mat, const vec3<T>& a)
{
	return vec3<T>(mat.m[0 + 0] * a.x + mat.m[0 + 1] * a.y + mat.m[0 + 2] * a.z,
	    mat.m[3 + 0] * a.x + mat.m[3 + 1] * a.y + mat.m[3 + 2] * a.z,
	    mat.m[6 + 0] * a.x + mat.m[6 + 1] * a.y + mat.m[6 + 2] * a.z);
}

template <class T>
const inline mat3<T> operator*(const mat3<T>& mat, const mat3<T>& a)
{
	T p[9];
	for (uint32_t i = 0; i < 3; i++) {
		for (uint32_t j = 0; j < 3; j++) {
			p[i * 3 + j] = mat.m[i * 3 + 0] * a.m[0 * 3 + j] + mat.m[i * 3 + 1] * a.m[1 * 3 + j] + mat.m[i * 3 + 2] * a.m[2 * 3 + j];
		}
	}

	return mat3<T>(p);
}

// alias

using fmat3 = mat3<float>;
using dmat3 = mat3<float>;

////mat3x2

template <class T>
class mat32 {

    public:
	T m[6];
	//m[0],m[1]
	//m[2],m[3]
	//m[4],m[5]

	//constructor
	inline mat32(const vec3<T>& a, const vec3<T>& b)
	{
		m[0] = a.x;
		m[2] = a.y;
		m[4] = a.z;
		m[1] = b.x;
		m[3] = b.y;
		m[5] = b.z;
	}
	inline mat32(const T (&a)[6])
	{
		for (uint32_t i = 0; i < 6; i++)
			m[i] = a[i];
	}

	//memberfunc
	inline T sqlength() const
	{
		T hoge = 0.0;
		for (uint32_t i = 0; i < 6; i++) {
			hoge += m[i] * m[i];
		}
		return hoge;
	}

	//cast operator
	template <class U>
	operator mat32<U>() const
	{
		U cm[6];
		for (uint32_t i = 0; i < 6; i++) {
			cm[i] = U(m[i]);
		}
		return mat32<U>(cm);
	}
};

template <class T>
std::ostream& operator<<(std::ostream& os, const mat32<T>& mat);

//operator

template <class T>
const inline mat32<T> operator+(const mat32<T>& mat0, const mat32<T>& mat1)
{
	T p[6];
	for (uint32_t i = 0; i < 6; i++)
		p[i] = mat0.m[i] + mat1.m[i];
	return mat32<T>(p);
}

template <class T, class U>
const inline mat32<T> operator*(const mat32<T>& mat, const U& a)
{
	T p[6];
	for (uint32_t i = 0; i < 6; i++)
		p[i] = mat.m[i] * a;
	return mat32<T>(p);
}

template <class T, class U>
const inline mat32<T> operator*(const U& a, const mat32<T>& mat)
{
	T p[6];
	for (uint32_t i = 0; i < 6; i++)
		p[i] = mat.m[i] * a;
	return mat32<T>(p);
}

template <class T>
const inline mat32<T> operator*(const mat32<T>& mat, const mat2<T>& a)
{
	T p[6];
	for (uint32_t i = 0; i < 3; i++) {
		for (uint32_t j = 0; j < 2; j++) {
			p[i * 2 + j] = mat.m[2 * i + 0] * a.m[2 * 0 + j] + mat.m[2 * i + 1] * a.m[2 * 1 + j];
		}
	}

	return mat32<T>(p);
}

template <class T>
const inline vec3<T> operator*(const mat32<T>& mat, const vec2<T>& a)
{
	return vec3<T>(
	    mat.m[2 * 0 + 0] * a.x + mat.m[2 * 0 + 1] * a.y,
	    mat.m[2 * 1 + 0] * a.x + mat.m[2 * 1 + 1] * a.y,
	    mat.m[2 * 2 + 0] * a.x + mat.m[2 * 2 + 1] * a.y);
}

using fmat32 = mat32<float>;
using dmat32 = mat32<double>;

////mat4

template <class T>
class mat4 {
    public:
	T m[16];
	// m[ 0],m[ 1],m[ 2],m[ 3]
	// m[ 4],m[ 5],m[ 6],m[ 7]
	// m[ 8],m[ 9],m[10],m[11]
	// m[12],m[13],m[14],m[15]

	//constructor
	inline mat4(const vec4<T>& a, const vec4<T>& b, const vec4<T>& c, const vec4<T>& d)
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
	inline mat4(const T (&a)[16])
	{
		for (uint32_t i = 0; i < 16; i++)
			m[i] = a[i];
	}

	inline mat4(void)
	{
		for (uint32_t i = 0; i < 16; i++)
			m[i] = 0.0;
	}

	template <class U>
	mat4(const mat3<U>& mat, const vec3<U>& vec)
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

	//menber func

	inline mat4<T> transpose() const
	{
		T data[16];

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

		return mat4<T>(data);
	}

	//cast operator

	template <class U>
	operator mat4<U>() const
	{
		U cm[16];
		for (uint32_t i = 0; i < 16; i++)
			cm[i] = U(m[i]);

		return mat4<U>(cm);
	}

	//static func
	inline static mat4<T> indentity()
	{
		mat4<T> temp = mat4<T>();
		temp.m[0]    = 1.0;
		temp.m[5]    = 1.0;
		temp.m[10]   = 1.0;
		temp.m[15]   = 1.0;
		return temp;
	}
};

//io operator

template <class T>
std::ostream& operator<<(std::ostream& os, const mat4<T>& mat);

//operator

template <class T>
const inline mat4<T> operator*(const mat4<T>& mat, const mat4<T>& a)
{
	T p[16];
	for (uint32_t i = 0; i < 4; i++) {
		for (uint32_t j = 0; j < 4; j++) {
			p[i * 4 + j] = mat.m[i * 4 + 0] * a.m[0 * 4 + j] + mat.m[i * 4 + 1] * a.m[1 * 4 + j] + mat.m[i * 4 + 2] * a.m[2 * 4 + j] + mat.m[i * 4 + 3] * a.m[3 * 4 + j];
		}
	}

	return mat4<T>(p);
}

template <class T>
const inline vec4<T> operator*(const mat4<T>& mat, const vec4<T>& a)
{
	return vec4<T>(
	    mat.m[4 * 0 + 0] * a.x + mat.m[4 * 0 + 1] * a.y + mat.m[4 * 0 + 2] * a.z + mat.m[4 * 0 + 3] * a.w,
	    mat.m[4 * 1 + 0] * a.x + mat.m[4 * 1 + 1] * a.y + mat.m[4 * 1 + 2] * a.z + mat.m[4 * 1 + 3] * a.w,
	    mat.m[4 * 2 + 0] * a.x + mat.m[4 * 2 + 1] * a.y + mat.m[4 * 2 + 2] * a.z + mat.m[4 * 2 + 3] * a.w,
	    mat.m[4 * 3 + 0] * a.x + mat.m[4 * 3 + 1] * a.y + mat.m[4 * 3 + 2] * a.z + mat.m[4 * 3 + 3] * a.w);
}

// alias

using fmat4 = mat4<float>;
using dmat4 = mat4<float>;

//others

template <class T>
mat3<T> vec3<T>::tensorproduct(const vec3<T>& v) const
{
	mat3<T> hoge;
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

template <class T>
mat3<T> vec3<T>::skew() const
{
	mat3<T> hoge;
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

////quaternion

template <class T>
class quaternion {
    public:
	T x, y, z, w;
	//scalar w
	//vector x,y,z

	inline quaternion(
	    const T& x,
	    const T& y,
	    const T& z,
	    const T& w)
	    : x(x)
	    , y(y)
	    , z(z)
	    , w(w)
	{
	}

	inline quaternion(const vec3<T>& v, const float& s)
	    : x(v.x)
	    , y(v.y)
	    , z(v.z)
	    , w(s)
	{
	}

	quaternion(const vec3<T>& v);

	inline quaternion(void)
	    : x(0.0)
	    , y(0.0)
	    , z(0.0)
	    , w(1.0)
	{
	}

	//menber function

	inline vec3<T> getv() const
	{
		return vec3<T>(x, y, z);
	}

	inline quaternion<T> conjugate() const
	{
		return quaternion<T>(-1.0 * vec3<T>(x, y, z), w);
	}

	inline T dot(const quaternion<T>& q) const
	{
		return x * q.x + y * q.y + z * q.z + w * q.w;
	}

	inline T sqlength() const
	{
		return x * x + y * y + z * z + w * w;
	}

	T length() const;

	quaternion<T> normalize() const;

	mat3<T> qtomat() const;

	inline quaternion<T> inverse() const
	{
		return quaternion<T>(-x, -y, -z, w);
	}

	quaternion<T> slerp(const quaternion<T>& q0, const quaternion<T>& q1, const float& t);
};

//io operator

template <class T>
std::ostream& operator<<(std::ostream& os, const quaternion<T>& q);

//operator

template <class T>
inline quaternion<T> operator+(const quaternion<T>& q0, const quaternion<T>& q1)
{
	return quaternion<T>(
	    q0.x + q1.x,
	    q0.y + q1.y,
	    q0.z + q1.z,
	    q0.w + q1.w);
}

template <class T>
inline quaternion<T> operator-(const quaternion<T>& q0, const quaternion<T>& q1)
{

	return quaternion<T>(
	    q0.x - q1.x,
	    q0.y - q1.y,
	    q0.z - q1.z,
	    q0.w - q1.w);
}

template <class T>
inline quaternion<T> operator-(const quaternion<T>& q)
{

	return quaternion<T>(
	    -q.x,
	    -q.y,
	    -q.z,
	    -q.w);
}

template <class T, class U>
inline quaternion<T> operator*(const quaternion<T>& q, const U& a)
{

	return quaternion<T>(
	    a * q.x,
	    a * q.y,
	    a * q.z,
	    a * q.w);
}

template <class T, class U>
inline quaternion<T> operator*(const U& a, const quaternion<T>& q)
{

	return quaternion<T>(
	    a * q.x,
	    a * q.y,
	    a * q.z,
	    a * q.w);
}

template <class T, class U>
inline quaternion<T> operator/(const quaternion<T>& q, const U& a)
{

	return quaternion<T>(
	    q.x / a,
	    q.y / a,
	    q.z / a,
	    q.w / a);
}

template <class T>
inline quaternion<T> operator*(const quaternion<T>& q0, const quaternion<T>& q1)
{

	T w	  = q0.w * q1.w - (q0.x * q1.x + q0.y * q1.y + q0.z * q1.z);
	vec3<T> v = q0.w * q1.getv() + q1.w * q0.getv() + q0.getv().cross(q1.getv());

	return quaternion<T>(v, w);
}

//alias

using fquaternion = quaternion<float>;
using dquaternion = quaternion<double>;
