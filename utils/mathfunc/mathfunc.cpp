#include "utils/mathfunc/mathfunc.hpp"
#include <iostream>
#include <cstdint>
#include <cmath>

//vec2

//翻訳単位の外で使われる型について実体化する。
template class vec2<float>;
template class vec2<double>;

//menber function

template float vec2<float>::length() const;
template double vec2<double>::length() const;

template vec2<float> vec2<float>::normalize() const;
template vec2<double> vec2<double>::normalize() const;

template vec2<float> vec2<float>::rot() const;
template vec2<double> vec2<double>::rot() const;

template <class T>
T vec2<T>::length() const
{
	return sqrt(this->sqlength());
}

template <class T>
vec2<T> vec2<T>::normalize() const
{
	double length = this->length();
	if (length < 0.000001)
		return vec2<T>(0.0, 0.0);
	return (*this) / length;
}

template <class T>
vec2<T> vec2<T>::rot() const
{
	return vec2<T>(-1.0 * (this->y), this->x);
}

//standard io

template std::ostream& operator<<(std::ostream& os, const vec2<float>& vec);
template std::ostream& operator<<(std::ostream& os, const vec2<double>& vec);

template <class T>
std::ostream& operator<<(std::ostream& os, const vec2<T>& vec)
{
	os << vec.x << " " << vec.y;
	return os;
}

////vec3

template class vec3<float>;
template class vec3<double>;

//menber function

template float vec3<float>::length() const;
template double vec3<double>::length() const;

template vec3<float> vec3<float>::normalize() const;
template vec3<double> vec3<double>::normalize() const;

template <class T>
T vec3<T>::length() const
{
	return sqrt(this->sqlength());
}

template <class T>
vec3<T> vec3<T>::normalize() const
{
	double length = this->length();
	if (length < 0.0000001)
		return vec3<T>(0.0, 0.0, 0.0);
	return (*this) / length;
}

//standard io

template std::ostream& operator<<(std::ostream& os, const vec3<float>& vec);
template std::ostream& operator<<(std::ostream& os, const vec3<double>& vec);

template <class T>
std::ostream& operator<<(std::ostream& os, const vec3<T>& vec)
{
	os << vec.x << " " << vec.y << " " << vec.z;
	return os;
}

template mat3<float> vec3<float>::rotation() const;
template mat3<double> vec3<double>::rotation() const;

template <class T>
mat3<T> vec3<T>::rotation() const
{
	T omega	 = this->length();
	T omega2 = this->sqlength();
	if (omega2 < 0.0000001)
		return mat3<T>::indentity();
	T cos = std::cos(omega);
	T sin = std::sin(omega);

	return cos * mat3<T>::indentity() + ((1 - cos) / omega2) * this->tensorproduct(*this) + (sin / omega) * this->skew();
}

////vec4

template class vec4<float>;
template class vec4<double>;

//standard io

template std::ostream& operator<<(std::ostream& os, const vec4<float>& vec);
template std::ostream& operator<<(std::ostream& os, const vec4<double>& vec);

template <class T>
std::ostream& operator<<(std::ostream& os, const vec4<T>& vec)
{
	os << vec.x << " " << vec.y << " " << vec.z << " " << vec.w;
	return os;
}

////mat2

template class mat2<float>;
template class mat2<double>;

template mat2<float>::mat2(const double& omega);
template mat2<double>::mat2(const double& omega);

template <class T>
mat2<T>::mat2(const double& omega)
{
	float cos = std::cos(omega);
	float sin = std::sin(omega);

	m[0] = cos;
	m[1] = -sin;
	m[2] = sin;
	m[3] = cos;
}

//standard io

template std::ostream& operator<<(std::ostream& os, const mat2<float>& mat);
template std::ostream& operator<<(std::ostream& os, const mat2<double>& mat);

template <class T>
std::ostream& operator<<(std::ostream& os, const mat2<T>& mat)
{
	for (uint32_t i = 0; i < 2; i++) {
		for (uint32_t j = 0; j < 2; j++) {
			os << mat.m[i * 2 + j] << " ";
		}
		os << std::endl;
	}
	return os;
}

////mat3

template class mat3<float>;
template class mat3<double>;

//standard io

template std::ostream& operator<<(std::ostream& os, const mat3<float>& mat);
template std::ostream& operator<<(std::ostream& os, const mat3<double>& mat);

template <class T>
std::ostream& operator<<(std::ostream& os, const mat3<T>& mat)
{
	for (uint32_t i = 0; i < 3; i++) {
		for (uint32_t j = 0; j < 3; j++) {
			os << mat.m[i * 3 + j] << " ";
		}
		os << std::endl;
	}
	return os;
}

template vec3<float> mat3<float>::mattorot() const;
template vec3<double> mat3<double>::mattorot() const;

template <class T>
vec3<T> mat3<T>::mattorot() const
{
	T omega = std::acos((this->trace() - 1) / 2.0);
	vec3<T> l;
	l.x = this->m[7] - this->m[5];
	l.y = this->m[2] - this->m[6];
	l.y = this->m[3] - this->m[1];

	l = l.normalize();
	return omega * l;
}

//static function

////mat32

template class mat32<float>;
template class mat32<double>;

template std::ostream& operator<<(std::ostream& os, const mat32<float>& mat);
template std::ostream& operator<<(std::ostream& os, const mat32<double>& mat);

template <class T>
std::ostream& operator<<(std::ostream& os, const mat32<T>& mat)
{
	for (uint32_t i = 0; i < 3; i++) {
		for (uint32_t j = 0; j < 2; j++) {
			os << mat.m[2 * i + j] << " ";
		}
		os << std::endl;
	}
	return os;
}

////mat4

template class mat4<float>;
template class mat4<double>;

//standard io

template std::ostream& operator<<(std::ostream& os, const mat4<float>& mat);
template std::ostream& operator<<(std::ostream& os, const mat4<double>& mat);

template <class T>
std::ostream& operator<<(std::ostream& os, const mat4<T>& mat)
{
	for (uint32_t i = 0; i < 4; i++) {
		for (uint32_t j = 0; j < 4; j++) {
			os << mat.m[i * 4 + j] << " ";
		}
		os << std::endl;
	}
	return os;
}

//static function

////quaternion

template class quaternion<float>;
template class quaternion<double>;

//constructor

template quaternion<float>::quaternion(const vec3<float>& v);
template quaternion<double>::quaternion(const vec3<double>& v);

template <class T>
quaternion<T>::quaternion(const vec3<T>& v)
{

	if (v.sqlength() < 0.00000000001) {
		x = 0.0;
		y = 0.0;
		z = 0.0;
		w = 1.0;
	} else {
		T O	  = v.length();
		vec3<T> l = v.normalize();

		T sin = std::sin(O / 2.0);
		x     = sin * l.x;
		y     = sin * l.y;
		z     = sin * l.z;
		w     = std::cos(O / 2.0);
	}
}

//menber function

template float quaternion<float>::length() const;
template double quaternion<double>::length() const;

template quaternion<float> quaternion<float>::normalize() const;
template quaternion<double> quaternion<double>::normalize() const;

template mat3<float> quaternion<float>::qtomat() const;
template mat3<double> quaternion<double>::qtomat() const;

template <class T>
T quaternion<T>::length() const
{
	return std::sqrt(this->sqlength());
}

template <class T>
quaternion<T> quaternion<T>::normalize() const
{
	double length = this->length();
	if (length < 0.00000000001)
		return quaternion<T>(0.0, 0.0, 0.0, 1.0);
	return (*this) / (length);
}

template <class T>
mat3<T> quaternion<T>::qtomat() const
{
	vec3<T> v = this->getv();
	T w	  = this->w;

	return (w * w - v.sqlength()) * mat3<T>::indentity() + 2.0 * w * v.skew() + 2.0 * v.tensorproduct(v);
}

//standard io

template std::ostream& operator<<(std::ostream& os, const quaternion<float>& q);
template std::ostream& operator<<(std::ostream& os, const quaternion<double>& q);

template <class T>
std::ostream& operator<<(std::ostream& os, const quaternion<T>& q)
{
	os << "( " << q.x << " , " << q.y << " , " << q.z << " , " << q.w << " ) " << std::endl;

	return os;
}

//

template quaternion<float> quaternion<float>::slerp(const quaternion<float>& q0, const quaternion<float>& q1, const float& t);
template quaternion<double> quaternion<double>::slerp(const quaternion<double>& q0, const quaternion<double>& q1, const float& t);

template <class T>
quaternion<T> quaternion<T>::slerp(const quaternion<T>& q, const float& t) const
{
	return quaternion<T>::slerp(*this, q, t);
}

template <class T>
quaternion<T> quaternion<T>::slerp(const quaternion<T>& q0, const quaternion<T>& q1, const float& t)
{
	if (q0.dot(q1) < 0.0) {
		double x = 0.5 * std::acos(-q0.dot(q1));
		if (std::abs(x) < 0.0000001)
			return q0;
		return (std::sin(t * x) / std::sin(x)) * -q1 + (std::sin(x - t * x) / std::sin(x)) * q0;
	} else {
		double x = 0.5 * std::acos(q0.dot(q1));
		if (std::abs(x) < 0.0000001)
			return q0;
		return (std::sin(t * x) / std::sin(x)) * q1 + (std::sin(x - t * x) / std::sin(x)) * q0;
	}
}
