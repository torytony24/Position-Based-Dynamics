#ifndef __VECTOR3_H__
#define __VECTOR3_H__

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>
#include "Util.h"

using namespace std;

//==============================================================================
//	Vector3 class declaration
//!		Class to store a three dimensional vector.
//!		Expected uses include point in R3, RGB color, etc.
template <typename T> class Vector3 {
public:
	////////////////////////////////////////////////////////////////////////////
	// Constructors
	inline Vector3();
	inline Vector3(T d);
	inline Vector3(T d0, T d1, T d2);

	inline Vector3(const Vector3 &d);
	inline Vector3(const T    *d);
	// d should point to a T[3] that will be copied

////////////////////////////////////////////////////////////////////////////
// Index operators
	inline T &operator[](unsigned int i);
	inline const T&  operator[](unsigned int i) const;

	inline T &operator()(unsigned int i);
	inline T  operator()(unsigned int i) const;

	////////////////////////////////////////////////////////////////////////////
	// Assignment and set
	inline Vector3 &set(T d);
	inline Vector3 &set(T d0, T d1, T d2);
	inline Vector3 &set(const Vector3 &d);
	inline Vector3 &set(const T *d);

	inline Vector3 &operator=(T d);
	inline Vector3 &operator=(const Vector3 &d);
	inline Vector3 &operator=(const T *d);

	////////////////////////////////////////////////////////////////////////////
	// Getter methods
	inline const Vector3 &zero();
	inline const Vector3 &half();
	inline const Vector3 &one();

	////////////////////////////////////////////////////////////////////////////
	// Comparison operators
	inline int operator==(const Vector3 &d) const;
	inline int operator!=(const Vector3 &d) const;

	inline int operator==(T d) const;
	inline int operator!=(T d) const;

	////////////////////////////////////////////////////////////////////////////
	// In place arithmetic
	inline Vector3 &operator+=(T d);
	inline Vector3 &operator-=(T d);
	inline Vector3 &operator*=(T d);
	inline Vector3 &operator/=(T d);

	inline Vector3 &operator+=(const Vector3 &d);
	inline Vector3 &operator-=(const Vector3 &d);
	inline Vector3 &operator*=(const Vector3 &d);
	inline Vector3 &operator/=(const Vector3 &d);
	// Component-wise operations

	inline Vector3 &maxSet(const Vector3 &d);
	inline Vector3 &minSet(const Vector3 &d);

	inline T        norm();

	inline void		polar(T& d, T& x_rot, T& y_rot) const;

	//inline friend Vector3<T> cross( const Vector3<T> &a, const Vector3<T> &b ) {
	//	return Vector3<T>(a.y * b.z - b.y * a.z, a.z * b.x - b.z * a.x, a.x * b.y - b.x * a.y);
	//}


	////////////////////////////////////////////////////////////////////////////
	// Static methods
	inline static unsigned int cycleAxis(unsigned int axis, int direction);


	/////////////////////////////////////////////////////////////////////////
	// general methods
	inline void print();

	// Public data members
	T x, y, z;
};

template <typename T> inline Vector3<T>::Vector3() {
	x = y = z = 0;
}

template <typename T> inline Vector3<T>::Vector3(T d) {
	x = y = z = d;
}

template <typename T> inline Vector3<T>::Vector3(T d0, T d1, T d2) {
	x = d0; y = d1; z = d2;
}

template <typename T> inline Vector3<T>::Vector3(const Vector3 &d) {
	x = d.x; y = d.y; z = d.z;
}

template <typename T> inline Vector3<T>::Vector3(const T *d) {
	x = d[0]; y = d[1]; z = d[2];
}

template <typename T> inline T &Vector3<T>::operator[](unsigned int i) {
	return (&x)[i];
}

template <typename T> inline T &Vector3<T>::operator()(unsigned int i) {
	return (&x)[i];
}

template <typename T> inline const T&  Vector3<T>::operator[](unsigned int i) const {
	return (&x)[i];
}
  
template <typename T> inline T  Vector3<T>::operator()(unsigned int i) const {
	return (&x)[i];
}

template <typename T> inline Vector3<T> &Vector3<T>::set(T d) {
	x = y = z = d;
	return (*this);
}
  
template <typename T> inline Vector3<T> &Vector3<T>::set(T d0, T d1, T d2) {
	x = d0; y = d1; z = d2;
	return (*this);
}

template <typename T> inline Vector3<T> &Vector3<T>::set(const Vector3 &d) {
	x = d.x; y = d.y; z = d.z;
	return (*this);
}
  
template <typename T> inline Vector3<T> &Vector3<T>::set(const T *d) {
	x = d[0]; y = d[1]; z = d[2];
	return (*this);
}

template <typename T> inline Vector3<T> &Vector3<T>::operator=(T d) {
	x = y = z = d;
	return (*this);
}
 
template <typename T> inline Vector3<T> &Vector3<T>::operator=(const Vector3 &d) {
	x = d.x; y = d.y; z = d.z;
	return (*this);
}

template <typename T> inline Vector3<T> &Vector3<T>::operator=(const T *d) {
	x = d[0]; y = d[1]; z = d[2];
	return (*this);
}

// Getter methods
template <typename T> inline const Vector3<T> &Vector3<T>::zero() {
	const Vector3<T> d(0, 0, 0);
	return d;
}

template <typename T> inline const Vector3<T> &Vector3<T>::half() {
	const Vector3<T> d(0.5f, 0.5f, 0.5f);
	return d;
}

template <typename T> inline const Vector3<T> &Vector3<T>::one() {
	const Vector3<T> d(1, 1, 1);
	return d;
}
	  
// Comparison operators
template <typename T> inline int Vector3<T>::operator==(const Vector3 &d) const {
	return ((x == d.x) && (y == d.y) && (z == d.z));
}

template <typename T> inline int Vector3<T>::operator!=(const Vector3 &d) const {
	return ((x != d.x) || (y != d.y) || (z != d.z));
}
  
template <typename T> inline int Vector3<T>::operator==(T d) const {
	return ((x == d) && (y == d) && (z == d));
}

template <typename T> inline int Vector3<T>::operator!=(T d) const {
	return ((x != d) || (y != d) || (z != d));
}

// In place arithmetic
template <typename T> inline Vector3<T> &Vector3<T>::operator+=(T d) {
	x += d; y += d; z += d;
	return (*this);
}
    
template <typename T> inline Vector3<T> &Vector3<T>::operator-=(T d) {
	x -= d; y -= d; z -= d;
	return (*this);
}

template <typename T> inline Vector3<T> &Vector3<T>::operator*=(T d) {
	x *= d; y *= d; z *= d;
	return (*this);
}

template <typename T> inline Vector3<T> &Vector3<T>::operator/=(T d) {
	x /= d; y /= d; z /= d;
	return (*this);
}

// In place arithmetic; component wise operations
template <typename T> inline Vector3<T> &Vector3<T>::operator+=(const Vector3 &d) {
	x += d.x; y += d.y; z += d.z;
	return (*this);
}

template <typename T> inline Vector3<T> &Vector3<T>::operator-=(const Vector3 &d) {
	x -= d.x; y -= d.y; z -= d.z;
	return (*this);
}

template <typename T> inline Vector3<T> &Vector3<T>::operator*=(const Vector3 &d) {
	x *= d.x; y *= d.y; z *= d.z;
	return (*this);
}

template <typename T> inline Vector3<T> &Vector3<T>::operator/=(const Vector3 &d) {
	x /= d.x; y /= d.y; z /= d.z;
	return (*this);
}

// Setter as maximum and minimum
//! Setter as component-wise maximum.
template <typename T> inline Vector3<T> &Vector3<T>::maxSet(const Vector3 &d) {
	if (d.x > x) x = d.x;
	if (d.y > y) y = d.y;
	if (d.z > z) z = d.z;
	return (*this);
}

//! Setter as component-wise minimum.
template <typename T> inline Vector3<T> &Vector3<T>::minSet(const Vector3 &d) {
	if (d.x < x) x = d.x;
	if (d.y < y) y = d.y;
	if (d.z < z) z = d.z;
	return (*this);
}

template <typename T> inline T Vector3<T>::norm() {
	return sqrt(x*x + y * y + z * z);
}

//-------------------------------------------------------------------

template <typename T> Vector3<T> operator-(const Vector3<T> &d) {
	return Vector3<T>(-d.x, -d.y, -d.z);
}

//-------------------------------------------------------------------

template <typename T> inline Vector3<T> operator+(const Vector3<T> &a, const Vector3<T> &b) {
	return Vector3<T>(a.x + b.x, a.y + b.y, a.z + b.z);
}

template <typename T> inline Vector3<T> operator-(const Vector3<T> &a, const Vector3<T> &b) {
	return Vector3<T>(a.x - b.x, a.y - b.y, a.z - b.z);
}

template <typename T> inline T operator*(const Vector3<T> &a, const Vector3<T> &b) {
	return (T)(a.x * b.x + a.y * b.y + a.z * b.z);
}

template <typename T> inline Vector3<T> operator/(const Vector3<T> &a, const Vector3<T> &b) {
	return Vector3<T>(a.x / b.x, a.y / b.y, a.z / b.z);
}

//-------------------------------------------------------------------

template <typename T> inline Vector3<T> operator+(const Vector3<T> &a, T b) {
	return Vector3<T>(a.x + b, a.y + b, a.z + b);
}

template <typename T> inline Vector3<T> operator-(const Vector3<T> &a, T b) {
	return Vector3<T>(a.x - b, a.y - b, a.z - b);
}

template <typename T> inline Vector3<T> operator*(const Vector3<T> &a, T b) {
	return Vector3<T>(a.x * b, a.y * b, a.z * b);
}

template <typename T> inline Vector3<T> operator/(const Vector3<T> &a, T b) {
	return Vector3<T>(a.x / b, a.y / b, a.z / b);
}

//-------------------------------------------------------------------

template <typename T> inline Vector3<T> operator+(T a, const Vector3<T> &b) {
	return Vector3<T>(a + b.x, a + b.y, a + b.z);
}

template <typename T> inline Vector3<T> operator-(T a, const Vector3<T> &b) {
	return Vector3<T>(a - b.x, a - b.y, a - b.z);
}

template <typename T> inline Vector3<T> operator*(T a, const Vector3<T> &b) {
	return Vector3<T>(a * b.x, a * b.y, a * b.z);
}

template <typename T> inline Vector3<T> operator/(T a, const Vector3<T> &b) {
	return Vector3<T>(a / b.x, a / b.y, a / b.z);
}

//-------------------------------------------------------------------

template <typename T> inline bool operator<(const Vector3<T> &a, T b) {
	return (a.x < b && a.y < b && a.z < b);
}

template <typename T> inline bool operator>(const Vector3<T> &a, T b) {
	return (a.x > b && a.y > b && a.z > b);
}

//-------------------------------------------------------------------

template <typename T> std::ostream &operator<<(std::ostream &strm, const Vector3<T> &v) {
	strm << "[" << v.x << ", " << v.y << ", " << v.z << "]";
	return strm;
}

template <typename T> inline void Vector3<T>::polar(T& d, T& x_rot, T& y_rot) const {
	d = sqrt(x * x + y * y + z * z);
	float d2 = sqrt(x * x + z * z);
	if (!(fabs(d2) < 0.000001)) {
		x_rot = atan(-y / d2);
		y_rot = asin(x / d2);
		if (z < 0) y_rot = PI - y_rot;
	}
	else {
		x_rot = (y > 0) ? -PI / 2 : PI / 2;
		y_rot = 0;
	}
}

//-------------------------------------------------------------------

template <typename T> inline T l1Norm(const Vector3<T> &a) {
	return (((a.x > 0) ? a.x : -a.x) + ((a.y > 0) ? a.y : -a.y) + ((a.z > 0) ? a.z : -a.z));
}
  
template <typename T> inline T l2Norm(const Vector3<T> &a) {
	return (T)sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
}
  
template <typename T> inline T lInfNorm(const Vector3<T> &a) {
	return maximum(abs(a));
}

template <typename T> inline T mag(const Vector3<T> &a) {
	return (T)sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
}

template <typename T> inline T sqrMag(const Vector3<T> &a) {
	return (a.x * a.x + a.y * a.y + a.z * a.z);
}

template <typename T> inline void normalize(Vector3<T> &a) {
	T m = mag(a);
	if (m != 0) a /= m;
}

//-------------------------------------------------------------------

template <typename T> inline unsigned int dominantAxis(const Vector3<T> &v) {
	T x, y, z;
	if (v.x > 0) x = v.x; else x = -v.x;
	if (v.y > 0) y = v.y; else y = -v.y;
	if (v.z > 0) z = v.z; else z = -v.z;
	return (x > y) ? ((x > z) ? 0 : 2) : ((y > z) ? 1 : 2);
}

template <typename T> inline unsigned int subinantAxis(const Vector3<T> &v) {
	T x, y, z;
	if (v.x > 0) x = v.x; else x = -v.x;
	if (v.y > 0) y = v.y; else y = -v.y;
	if (v.z > 0) z = v.z; else z = -v.z;
	return (x < y) ? ((x < z) ? 0 : 2) : ((y < z) ? 1 : 2);
}

template <typename T> inline unsigned int midinantAxis(const Vector3<T> &v) {
	T x, y, z;
	if (v.x > 0) x = v.x; else x = -v.x;
	if (v.y > 0) y = v.y; else y = -v.y;
	if (v.z > 0) z = v.z; else z = -v.z;
	 unsigned int d = (x > y) ? ((x > z) ? 0 : 2) : ((y > z) ? 1 : 2);
	 unsigned int s = (x < y) ? ((x < z) ? 0 : 2) : ((y < z) ? 1 : 2);
	 unsigned int m;
	if (d == 0) {
		if (s != 1) m = 1;
		else m = 2;
	}
	else if (d == 1) {
		if (s != 0) m = 0;
		else m = 2;
	}
	else if (d == 2) {
		if (s != 0) m = 0;
		else m = 1;
	}
	return m;
}

//-------------------------------------------------------------------

template <typename T> inline T dot(const Vector3<T> &a, const Vector3<T> &b) {
	return (a.x * b.x + a.y * b.y + a.z * b.z);
}

//! Cross-product of two vectors.
template <typename T> inline Vector3<T> cross(const Vector3<T> &a, const Vector3<T> &b) {
	return Vector3<T>(a.y * b.z - b.y * a.z, a.z * b.x - b.z * a.x, a.x * b.y - b.x * a.y);
}

template <typename T> inline T dist(const Vector3<T> &a, const Vector3<T> &b) {
	return sqrt((a.x - b.x)*(a.x - b.x) + (a.y - b.y)*(a.y - b.y) + (a.z - b.z)*(a.z - b.z));
}

template <typename T> inline T angle(const Vector3<T> &a, const Vector3<T> &b) {
	return ACOS(a*b / (mag(a)*mag(b)));
}


//-------------------------------------------------------------------
//! Volume of parallelepipe generated by three vectors.
template <typename T> inline T box(const Vector3<T> &a, const Vector3<T> &b, const Vector3<T> &c) {
	return dot(cross(a, b), c);
}

//-------------------------------------------------------------------

template <typename T> inline Vector3<T> abs(const Vector3<T> &a) {
	return Vector3<T>(((a.x > 0) ? a.x : -a.x), ((a.y > 0) ? a.y : -a.y), ((a.z > 0) ? a.z : -a.z));
}

template <typename T> inline T sum(const Vector3<T> &a) {
	return a.x + a.y + a.z;
}

//-------------------------------------------------------------------

template <typename T> inline T maximum(const Vector3<T> &a) {
	return ((a.x > a.y) ? ((a.x > a.z) ? a.x : a.z) : (a.y > a.z) ? a.y : a.z);
}

template <typename T> inline T minimum(const Vector3<T> &a) {
	return ((a.x < a.y) ? ((a.x < a.z) ? a.x : a.z) : (a.y < a.z) ? a.y : a.z);
}

//-------------------------------------------------------------------

template <typename T> inline Vector3<T> maximum(const Vector3<T> &a, const Vector3<T> &b) {
	return Vector3<T>((a.x > b.x) ? a.x : b.x, (a.y > b.y) ? a.y : b.y, (a.z > b.z) ? a.z : b.z);
}
 
template <typename T> inline Vector3<T> minimum(const Vector3<T> &a, const Vector3<T> &b) {
	return Vector3<T>((a.x < b.x) ? a.x : b.x, (a.y < b.y) ? a.y : b.y, (a.z < b.z) ? a.z : b.z);
}
//-------------------------------------------------------------------

template <typename T>
inline Vector3<T> normalized(const Vector3<T>& v) {
	Vector3<T> ret(v);
	normalize(ret);
	return ret;
}

template <typename T>
inline Vector3<T> tangential(const Vector3<T>& v) {
	return normalized(
		cross((abs(v.x) < 1. ? Vector3<T>(1, 0, 0) : Vector3<T>(0, 1, 0)), v)
	);
}

template <typename T> inline void Vector3<T>::print() {
	cout << "[ " << x << " , " << y << " , " << z << " ]" << endl;
}

template <typename T> inline T dot(const std::vector<Vector3<T> > &a, const std::vector<Vector3<T> > &b) {
	const int N = (int)a.size();

	T tmp = (T)0.0;

#pragma omp parallel for reduction(+:tmp)
	for (int i = 0; i < N; i++) {
		tmp += (a[i] * b[i]);
	}
	return tmp;
}

template <typename T> inline void scale(const T s, const std::vector<Vector3<T> > &v, std::vector<Vector3<T> > &r) {
	const int N = (int)v.size();
#pragma omp parallel for
	for (int i = 0; i < N; i++) {
		r[i] = s * v[i];
	}
}

template <typename T> inline void scale_add(const T s, const std::vector<Vector3<T> > &v, std::vector<Vector3<T> > &r) {
	const int N = (int)v.size();
#pragma omp parallel for
	for (int i = 0; i < N; i++) {
		r[i] += s * v[i];
	}
}

#endif // __VECTOR3_H__

