#ifndef __VECTOR2_H__
#define __VECTOR2_H__

#include <cmath>
#include <cstdlib>
#include <iostream>

//============================================================================
//	Vector2 class definition
//!		Class to store a two dimensional vector. Expected uses include point or 
//!		velocity in R2, UVs in texture, etc.
template <typename T> class Vector2 {
public:
	//////////////////////////////////////////////////////////////////////////
	// Constructors
	inline Vector2();									//!< Default constructor.
	inline Vector2(T d);								//!< Constructor. 
														//!< @param d - default value.
	inline Vector2(T d0, T d1);						//!< Constructor.
														//!< @param d0 - 1st element.
														//!< @param d1 - 2nd element.

	inline Vector2(const Vector2 &d);					//!< Constructor.
														//!< @param d - default constant vector.
	inline Vector2(const T *d);						//!< Constructor
														//!< @param d - pointer of data Array to be assinged as default constant value.

	//////////////////////////////////////////////////////////////////////////
	// Index operators
	inline T &get(unsigned int i);					//!< Getter by reference.
														//!< @param i - index for i_th component.
	inline T  get(unsigned int i) const;				//!< Getter by value.
														//!< @param i - index for i_th component.

	inline T &operator[](unsigned int i);				//!< Getter by reference.
														//!< @param i - index for i_th component.
	inline T  operator[](unsigned int i) const;		//!< Getter by value.
														//!< @param i - index for i_th component.

	inline T &operator()(unsigned int i);				//!< Getter by reference.
														//!< @param i - index for i_th component.
	inline T  operator()(unsigned int i) const;		//!< Getter by value.
														//!< @param i - index for i_th component.

	//////////////////////////////////////////////////////////////////////////
	// Assignment and set
	inline Vector2 &set(T d);							//!< Global setter.
														//!< @param d - Constant value to be assigned.
	inline Vector2 &set(T d0, T d1);					//!< Component-wise setter.
														//!< @param d0 - Value to be set to 1st component.
														//!< @param d1 - Value to be set to 2nd component.

	inline Vector2 &set(const Vector2 &d);			//!< Setter sourced from Vector2.
														//!< @param d - Vector2 to be assinged.
	inline Vector2 &set(const T *d);					//!< Setter sourced from data array.
														//!< @param d - pointer of data Array to be assinged.
														// d should point to a T[2] that will be copied

	inline Vector2 &operator=(T d);					//!< Assignment operation.
														//!< @param d - value to be assigned.
	inline Vector2 &operator=(const Vector2 &d);		//!< Assignment operation.
														//!< @param d - constant Vector2 to be assigned.
	inline Vector2 &operator=(const T *d);			//!< Assignment operation.
														//!< @param d - constant pointer of data array to be assigned.

	//////////////////////////////////////////////////////////////////////////
	// Comparison operators
	inline int isEqual(const Vector2 &d) const;		//!< Equality comparison function.
														//!< @param d - constant Vector2 to be compared.
	inline int isNotEqual(const Vector2 &d) const;	//!< Inequality comparison function. 
														//!< @param d - constant Vector2 to be compared.
	inline int operator==(const Vector2 &d) const;	//!< Equality comparison operator.
														//!< @param d - constant Vector2 to be compared.
	inline int operator!=(const Vector2 &d) const;	//!< Inequality comparison operator.
														//!< @param d - constant Vector2 to be compared.

	inline int operator==(T d) const;					//!< Equality comparison operator.
														//!< @param d - value to be compared.
	inline int operator!=(T d) const;					//!< Inequality comparison operator.
														//!< @param d - value to be compared.

	//////////////////////////////////////////////////////////////////////////
	// In place arithmetic
	inline Vector2 &add(T d);							//!< Global addition function.
														//!< @param d - value to be added.
	inline Vector2 &subtract(T d);					//!< Global subtraction function.
														//!< @param d - value to be subtractd.
	inline Vector2 &multiply(T d);					//!< Global multiplication function.
														//!< @param d - value to be multiplicated.
	inline Vector2 &divide(T d);						//!< Global division function.
														//!< @param d - value to be divided.

	inline Vector2 &operator+=(T d);					//!< Global addition operator.
														//!< @param d - value to be added.
	inline Vector2 &operator-=(T d);					//!< Global subtraction operator.
														//!< @param d - value to be subtracted.
	inline Vector2 &operator*=(T d);					//!< Global multiplication operator.
														//!< @param d - value to be multiplicated.
	inline Vector2 &operator/=(T d);					//!< Global division operator.
														//!< @param d - value to be divided.

	//////////////////////////////////////////////////////////////////////////
	// In place arithmetic; componentwise operations
	inline Vector2 &operator+=(const Vector2 &d);		//!< Component-wise addition operator.
														//!< @param d - constant Vector2 to be added.
	inline Vector2 &operator-=(const Vector2 &d);		//!< Component-wise subtraction operator.
														//!< @param d - constant Vector2 to be subtracted.
	inline Vector2 &operator*=(const Vector2 &d);		//!< Component-wise multiplication operator.
														//!< @param d - constant Vector2 to be multiplied.
	inline Vector2 &operator/=(const Vector2 &d);		//!< Component-wise division operator.
														//!< @param d - constant Vector2 to be divided.

	//////////////////////////////////////////////////////////////////////////
	// Setter as maximum and minimum
	inline Vector2 &maxSet(const Vector2 &d);			//! Component-wise maximum setter fuction.
														//!< @param d - constant Vector2.
														// Sets val[i] = minimum(val[i],d[i])
	inline Vector2 &minSet(const Vector2 &d);			//! Component-wise minimum setter fuction.
														//!< @param d - constant Vector2.
														// Sets val[i] = maximum(val[i],d[i])

	inline T norm();

	//////////////////////////////////////////////////////////////////////////
	// Static methods
	inline static unsigned int cycleAxis(unsigned int axis, int direction);

	//////////////////////////////////////////////////////////////////////////
	// Define Components
	enum Index { x = 0, y = 1, u = 0, v = 1 };

private:
	T val[2];
};
	
template <typename T> inline T &Vector2<T>::get(unsigned int i) {
	DB_CHECK(i < 2);
	return val[i];
}

template <typename T> inline T  Vector2<T>::get(unsigned int i) const {
	DB_CHECK(i < 2);
	return val[i];
}

template <typename T> inline T &Vector2<T>::operator[](unsigned int i) {
	DB_CHECK(i < 2);
	return val[i];
}

template <typename T> inline T &Vector2<T>::operator()(unsigned int i) {
	DB_CHECK(i < 2);
	return val[i];
}

template <typename T> inline T  Vector2<T>::operator[](unsigned int i) const {
	DB_CHECK(i < 2);
	return val[i];
}
  
template <typename T> inline T  Vector2<T>::operator()(unsigned int i) const {
	DB_CHECK(i < 2);
	return val[i];
}

template <typename T> inline Vector2<T>::Vector2() {
	val[0] = val[1] = (T)0.0;
}

template <typename T> inline Vector2<T>::Vector2(T d) {
	val[0] = val[1] = d;
}

template <typename T> inline Vector2<T>::Vector2(T d0, T d1) {
	val[0] = d0;
	val[1] = d1;
}

template <typename T> inline Vector2<T>::Vector2(const Vector2 &d) {
	val[0] = d[0];
	val[1] = d[1];
}

template <typename T> inline Vector2<T>::Vector2(const T *d) {
	val[0] = d[0];
	val[1] = d[1];
}

template <typename T> inline Vector2<T> &Vector2<T>::set(T d) {
	val[0] = d;
	val[1] = d;
	return (*this);
}
  
template <typename T> inline Vector2<T> &Vector2<T>::set(T d0, T d1) {
	val[0] = d0;
	val[1] = d1;
	return (*this);
}

template <typename T> inline Vector2<T> &Vector2<T>::set(const Vector2 &d) {
	val[0] = d[0];
	val[1] = d[1];
	return (*this);
}
  
template <typename T> inline Vector2<T> &Vector2<T>::set(const T *d) {
	val[0] = d[0];
	val[1] = d[1];
	return (*this);
}
 
template <typename T> inline Vector2<T> &Vector2<T>::operator=(T d) {
	return set(d);
}
 
template <typename T> inline Vector2<T> &Vector2<T>::operator=(const Vector2 &d) {
	return set(d);
}

template <typename T> inline Vector2<T> &Vector2<T>::operator=(const T *d) {
	return set(d);
}

template <typename T> inline int Vector2<T>::isEqual(const Vector2 &d) const {
	return ((val[0] == d[0]) && (val[1] == d[1]));
}

template <typename T> inline int Vector2<T>::isNotEqual(const Vector2 &d) const {
	return ((val[0] != d[0]) || (val[1] != d[1]));
}
 
template <typename T> inline int Vector2<T>::operator==(const Vector2 &d) const {
	return ((val[0] == d[0]) && (val[1] == d[1]));
}

template <typename T> inline int Vector2<T>::operator!=(const Vector2 &d) const {
	return ((val[0] != d[0]) || (val[1] != d[1]));
}
  
template <typename T> inline int Vector2<T>::operator==(T d) const {
	return ((val[0] == d) && (val[1] == d));
}

template <typename T> inline int Vector2<T>::operator!=(T d) const {
	return ((val[0] != d) || (val[1] != d));
}

template <typename T> inline Vector2<T> &Vector2<T>::add(T d) {
	val[0] += d;
	val[1] += d;
	return (*this);
}
    
template <typename T> inline Vector2<T> &Vector2<T>::subtract(T d) {
	val[0] -= d;
	val[1] -= d;
	return (*this);
}

template <typename T> inline Vector2<T> &Vector2<T>::multiply(T d) {
	val[0] *= d;
	val[1] *= d;
	return (*this);
}

template <typename T> inline Vector2<T> &Vector2<T>::divide(T d) {
	val[0] /= d;
	val[1] /= d;
	return (*this);
}

template <typename T> inline Vector2<T> &Vector2<T>::operator+= (T d) {
	val[0] += d;
	val[1] += d;
	return (*this);
}
    
template <typename T> inline Vector2<T> &Vector2<T>::operator-=(T d) {
	val[0] -= d;
	val[1] -= d;
	return (*this);
}

template <typename T> inline Vector2<T> &Vector2<T>::operator*=(T d) {
	val[0] *= d;
	val[1] *= d;
	return (*this);
}

template <typename T> inline Vector2<T> &Vector2<T>::operator/=(T d) {
	val[0] /= d;
	val[1] /= d;
	return (*this);
}

template <typename T> inline Vector2<T> &Vector2<T>::operator+=(const Vector2 &d) {
	val[0] += d[0];
	val[1] += d[1];
	return (*this);
}

template <typename T> inline Vector2<T> &Vector2<T>::operator-=(const Vector2 &d) {
	val[0] -= d[0];
	val[1] -= d[1];
	return (*this);
}

template <typename T> inline Vector2<T> &Vector2<T>::operator*=(const Vector2 &d) {
	val[0] *= d[0];
	val[1] *= d[1];
	return (*this);
}

template <typename T> inline Vector2<T> &Vector2<T>::operator/=(const Vector2 &d) {
	val[0] /= d[0];
	val[1] /= d[1];
	return (*this);
}

template <typename T> inline Vector2<T> &Vector2<T>::maxSet(const Vector2 &d) {
	if (d[0] > val[0]) val[0] = d[0];
	if (d[1] > val[1]) val[1] = d[1];
	return (*this);
}

template <typename T> inline Vector2<T> &Vector2<T>::minSet(const Vector2 &d) {
	if (d[0] < val[0]) val[0] = d[0];
	if (d[1] < val[1]) val[1] = d[1];
	return (*this);
}

template <typename T> inline unsigned int Vector2<T>::cycleAxis(unsigned int axis, int direction) {
	switch (axis + direction) {
	case 0: case 2: case 4: return 0;
	case 1: case 3: case 5: return 1;
	default: return (axis + direction) % 2;
	}
}

template <typename T> inline T Vector2<T>::norm() {
	return (T)sqrt(val[0] * val[0] + val[1] * val[1]);
}

//! Multiply -1 to all component.
//! @param a - Vector2 to be modified.
template <typename T> inline Vector2<T> neg(const Vector2<T> &a) {
	return Vector2<T>(-a[0], -a[1]);
}

//! Multiply -1 to all component.
//! @param a - Vector2 to be modified.
template <typename T> inline Vector2<T> operator-(const Vector2<T> &a) {
	return Vector2<T>(-a[0], -a[1]);
}

//! Component-wise summation of two vectors.
//! @param a - Vector2 to be added.
//! @param b - Vector2 to added.
template <typename T> inline Vector2<T> sum(const Vector2<T> &a, const Vector2<T> &b) {
	return Vector2<T>(a[0] + b[0], a[1] + b[1]);
}

//! Component-wise subtraction of two vectors.
//! @param a - Vector2 to be subtracted.
//! @param b - Vector2 to subtract.
template <typename T> inline Vector2<T> dif(const Vector2<T> &a, const Vector2<T> &b) {
	return Vector2<T>(a[0] - b[0], a[1] - b[1]);
}

//! Component-wise multiplication of two vectors.
//! @param a - Vector2 to be multiplicated.
//! @param b - Vector2 to multiply.
template <typename T> inline Vector2<T> mul(const Vector2<T> &a, const Vector2<T> &b) {
	return Vector2<T>(a[0] * b[0], a[1] * b[1]);
}

//! Component-wise division of two vectors.
//! @param a - Vector2 to be divided.
//! @param b - Vector2 to divide.
template <typename T> inline Vector2<T> div(const Vector2<T> &a, const Vector2<T> &b) {
	return Vector2<T>(a[0] / b[0], a[1] / b[1]);
}

//! Component-wise summation of two vectors.
//! @param a - Vector2 to be added.
//! @param b - Vector2 to added.
template <typename T> inline Vector2<T> operator+(const Vector2<T> &a, const Vector2<T> &b) {
	return Vector2<T>(a[0] + b[0], a[1] + b[1]);
}

//! Component-wise subtraction of two vectors.
//! @param a - Vector2 to be subtracted.
//! @param b - Vector2 to subtract.
template <typename T> inline Vector2<T> operator-(const Vector2<T> &a, const Vector2<T> &b) {
	return Vector2<T>(a[0] - b[0], a[1] - b[1]);
}

//! Component-wise multiplication of two vectors.
//! @param a - Vector2 to be multiplicated.
//! @param b - Vector2 to multiply.
template <typename T> inline T operator*(const Vector2<T> &a, const Vector2<T> &b) {
	return (a[0] * b[0]) + (a[1] * b[1]);
}

//! Component-wise division of two vectors.
//! @param a - Vector2 to be divided.
//! @param b - Vector2 to divide.
template <typename T> inline Vector2<T> operator/(const Vector2<T> &a, const Vector2<T> &b) {
	return Vector2<T>(a[0] / b[0], a[1] / b[1]);
}

//! Component-wise summation of vector and scalar.
//! @param a - Vector2 to be added.
//! @param b - scalar to added.
template <typename T> inline Vector2<T> sum(const Vector2<T> &a, T b) {
	return Vector2<T>(a[0] + b, a[1] + b);
}

//! Component-wise subtraction of vector and scalar.
//! @param a - Vector2 to be added.
//! @param b - scalar to added.
template <typename T> inline Vector2<T> dif(const Vector2<T> &a, T b) {
	return Vector2<T>(a[0] - b, a[1] - b);
}

//! Component-wise multiplication of vector and scalar.
//! @param a - Vector2 to be added.
//! @param b - scalar to added.
template <typename T> inline Vector2<T> mul(const Vector2<T> &a, T b) {
	return Vector2<T>(a[0] * b, a[1] * b);
}

//! Component-wise division of vector and scalar.
//! @param a - Vector2 to be added.
//! @param b - scalar to added.
template <typename T> inline Vector2<T> div(const Vector2<T> &a, T b) {
	return Vector2<T>(a[0] / b, a[1] / b);
}

//! Component-wise summation of vector and scalar.
//! @param a - Vector2 to be added.
//! @param b - scalar to added.
template <typename T> inline Vector2<T> operator+(const Vector2<T> &a, T b) {
	return Vector2<T>(a[0] + b, a[1] + b);
}

//! Component-wise subtraction of vector and scalar.
//! @param a - Vector2 to be added.
//! @param b - scalar to added.
template <typename T> inline Vector2<T> operator-(const Vector2<T> &a, T b) {
	return Vector2<T>(a[0] - b, a[1] - b);
}

//! Component-wise multiplication of vector and scalar.
//! @param a - Vector2 to be added.
//! @param b - scalar to added.
template <typename T> inline Vector2<T> operator*(const Vector2<T> &a, T b) {
	return Vector2<T>(a[0] * b, a[1] * b);
}

//! Component-wise division of vector and scalar.
//! @param a - Vector2 to be added.
//! @param b - scalar to added.
template <typename T> inline Vector2<T> operator/(const Vector2<T> &a, T b) {
	return Vector2<T>(a[0] / b, a[1] / b);
}

//! Component-wise summation of scalar and vector.
//! @param a - scalar to be added.
//! @param b - Vector2 to added.
template <typename T> inline Vector2<T> sum(T a, const Vector2<T> &b) {
	return Vector2<T>(a + b[0], a + b[1]);
}

//! Component-wise subtraction of scalar and vector.
//! @param a - scalar to be added.
//! @param b - Vector2 to added.
template <typename T> inline Vector2<T> dif(T a, const Vector2<T> &b) {
	return Vector2<T>(a - b[0], a - b[1]);
}

//! Component-wise multiplication of scalar and vector.
//! @param a - scalar to be added.
//! @param b - Vector2 to added.
template <typename T> inline Vector2<T> mul(T a, const Vector2<T> &b) {
	return Vector2<T>(a * b[0], a * b[1]);
}

//! Component-wise division of scalar and vector.
//! @param a - scalar to be added.
//! @param b - Vector2 to added.
template <typename T> inline Vector2<T> div(T a, const Vector2<T> &b) {
	return Vector2<T>(a / b[0], a / b[1]);
}

//! Component-wise summation of scalar and vector.
//! @param a - scalar to be added.
//! @param b - Vector2 to added.
template <typename T> inline Vector2<T> operator+(T a, const Vector2<T> &b) {
	return Vector2<T>(a + b[0], a + b[1]);
}

//! Component-wise subtraction of scalar and vector.
//! @param a - scalar to be added.
//! @param b - Vector2 to added.
template <typename T> inline Vector2<T> operator-(T a, const Vector2<T> &b) {
	return Vector2<T>(a - b[0], a - b[1]);
}

//! Component-wise multiplication of scalar and vector.
//! @param a - scalar to be added.
//! @param b - Vector2 to added.
template <typename T> inline Vector2<T> operator*(T a, const Vector2<T> &b) {
	return Vector2<T>(a * b[0], a * b[1]);
}

//! Component-wise division of scalar and vector.
//! @param a - scalar to be added.
//! @param b - Vector2 to added.
template <typename T> inline Vector2<T> operator/(T a, const Vector2<T> &b) {
	return Vector2<T>(a / b[0], a / b[1]);
}

//! Operator overloading for data output to stream.
//! @param strm - Output stream.
//! @param v - Vector2 to be output.
template <typename T> std::ostream &operator<<(std::ostream &strm, const Vector2<T> &v) {
	strm << "[" << v[0] << ", " << v[1] << "]";
	return strm;
}

// Computes the l1, l2 or lInfinity norm of a
//! L1 norm of Vector2.
//! @param a - Vector2 to be calculated.
template <typename T> inline T l1Norm(const Vector2<T> &a) {
	return (((a[0] > 0) ? a[0] : -a[0]) + ((a[1] > 0) ? a[1] : -a[1]));
}
 
//! L2 norm of Vector2.
//! @param a - Vector2 to be calculated.
template <typename T> inline T l2Norm(const Vector2<T> &a) {
	return (T)sqrt(a[0] * a[0] + a[1] * a[1]);
}

//! L-infinity norm of Vector2.
//! @param a - Vector2 to be calculated.
template <typename T> inline T lInfNorm(const Vector2<T> &a) {
	return maximum(abs(a));
}

// mag is the l2Norm or magnitude of the vector
//! Magnitude of Vector2. It is equivalant to the L2 norm.
//! @param a - Vector2 to be calculated.
template <typename T> inline T mag(const Vector2<T> &a) {
	return (T)sqrt(abs(a[0] * a[0] + a[1] * a[1]));
}

//! Square of magnitude of Vector2. Its fast to compute.
//! @param a - Vector2 to be calculated.
template <typename T> inline T sqrMag(const Vector2<T> &a) {
	return (a[0] * a[0] + a[1] * a[1]);
}

// Sets a = a/mag(a)
//! Normalization. The vector is set to unit size.
//! @param a - Vector2 to be calculated.
template <typename T> inline void normalize(Vector2<T> &a) {
	T m = mag(a);
	if (m != 0) a /= m;
}

//! The index of the component with the largest value. If multiple components have the same value, then the results are not unique.
//! @param v - Vector2 to be referred.
template <typename T> inline unsigned int dominantAxis(const Vector2<T> &v) {
	T x, y;
	if (v[0] > 0) x = v[0]; else x = -v[0];
	if (v[1] > 0) y = v[1]; else y = -v[1];
	return (x > y) ? 0 : 1;
}

//! The index of the component with the smallest value. If multiple components have the same value, then the results are not unique.
//! @param v - Vector2 to be referred.
template <typename T> inline unsigned int subinantAxis(const Vector2<T> &v) {
	T x, y;
	if (v[0] > 0) x = v[0]; else x = -v[0];
	if (v[1] > 0) y = v[1]; else y = -v[1];
	return (x < y) ? 0 : 1;
}

//! Dot product of two vectors.
template <typename T> inline T dot(const Vector2<T> &a, const Vector2<T> &b) {
	return (a[0] * b[0] + a[1] * b[1]);
}

//! Cross product of two vectors.
template <typename T> inline T cross(const Vector2<T> &a, const Vector2<T> &b) {
	return (a[0] * b[1] - b[0] * a[1]);
}

//! Absolute value of a vector.
template <typename T> inline Vector2<T> abs(const Vector2<T> &a) {
	return  Vector2<T>(((a[0] > 0) ? a[0] : -a[0]), ((a[1] > 0) ? a[1] : -a[1]));
}

//! Summation of two components of a vector.
template <typename T> inline T sum(const Vector2<T> &a) {
	return a[0] + a[1];
}

//! Value of the component with the largest value. If multiple components have the same value, then the results are not unique.
template <typename T> inline T maximum(const Vector2<T> &a) {
	return ((a[0] > a[1]) ? a[0] : a[1]);
}

//! Value of the component with the largest value. If multiple components have the same value, then the results are not unique.
template <typename T> inline T minimum(const Vector2<T> &a) {
	return ((a[0] < a[1]) ? a[0] : a[1]);
}

//! A new vector by taking the max component in the x,y direction from a and b.  ie: r[i] = maximum(a[i],b[i]). Note: signed values are used.
template <typename T> inline Vector2<T> maximum(const Vector2<T> &a, const Vector2<T> &b) {
	return Vector2<T>((a[0] > b[0]) ? a[0] : b[0], (a[1] > b[1]) ? a[1] : b[1]);
}

//! A new vector by taking the min component in the x,y direction from a and b.
template <typename T> inline Vector2<T> minimum(const Vector2<T> &a, const Vector2<T> &b) {
	return Vector2<T>((a[0] < b[0]) ? a[0] : b[0], (a[1] < b[1]) ? a[1] : b[1]);
}

#endif // __VECTOR2_H__
