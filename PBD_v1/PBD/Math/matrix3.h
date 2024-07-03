#pragma once

#include <cstring>
#include <vector>

#include "vector.h"
#include "Util.h"

//==============================================================================
//	Matrix3 class declaration.
//!	3x3 matrix class.
template <typename T> class Matrix3 {
public:
	////////////////////////////////////////////////////////////////////////////
	// Constructors
	Matrix3();
	Matrix3(const Matrix3<T> &d);
	Matrix3(T m00, T m10, T m20, T m01, T m11, T m21, T m02, T m12, T m22); // column-major
	Matrix3(T m[9]);
	Matrix3(Vector3<T> col0, Vector3<T> col1, Vector3<T> col2);

	////////////////////////////////////////////////////////////////////////////
	// Index operators
	inline const T& operator()(unsigned int i, unsigned int j) const;

	////////////////////////////////////////////////////////////////////////////
	// Setter and Assignment 
	inline Matrix3 &set(T d);
	inline Matrix3 &set(T m00, T m10, T m20, T m01, T m11, T m21, T m02, T m12, T m22);
	inline Matrix3 &set(const Matrix3 &d);
	inline Matrix3 &set(const T *d);
	inline Matrix3 &operator=(T d);
	inline Matrix3 &operator=(const Matrix3 &d);
	inline Matrix3 &operator=(const T *d);
	inline Matrix3 &zero();

	//////////////////////////////////////////////////////////////////////////
	// In place arithmetic
	inline Matrix3 &operator+=(const Matrix3 &rhs);
	inline Matrix3 &operator-=(const Matrix3 &rhs);
	inline Matrix3 &operator*=(const Matrix3 &rhs);
	inline Matrix3 &operator*=(T d);
	inline Matrix3 &operator/=(T d);

	////////////////////////////////////////////////////////////////////////////
	// Linear algebra methods
	inline T determinant();
	inline T inverse(Matrix3<T>& invM);
	inline Matrix3 &transpose();
	inline Matrix3<T> transp() const;


	/////////////////////////////////////////////////////////////////////////
	// general methods
	inline bool isAlmostZero();
	inline void Identity();

public:
	T data[9]; // row-major
};


template <typename T> inline Matrix3<T>::Matrix3() {
	memset(data, 0, 9 * sizeof(T));
	data[0] = 1.0;
	data[4] = 1.0;
	data[8] = 1.0;
}

template <typename T> inline Matrix3<T>::Matrix3(const Matrix3<T> &d) {
	memcpy(data, d.data, 9 * sizeof(T));
}

template <typename T> inline Matrix3<T>::Matrix3(T m00, T m10, T m20, T m01, T m11, T m21, T m02, T m12, T m22) {
	data[0] = m00; // column 1 
	data[3] = m10;
	data[6] = m20;

	data[1] = m01; // column 2 
	data[4] = m11;
	data[7] = m21;

	data[2] = m02; // column 3 
	data[5] = m12;
	data[8] = m22;
}

template <typename T> inline Matrix3<T>::Matrix3(T d[9]) {
	memcpy(data, d, 9 * sizeof(T));
}

template <typename T> inline Matrix3<T>::Matrix3(Vector3<T> col0, Vector3<T> col1, Vector3<T> col2) {
	data[0] = col0[0];
	data[3] = col0[1];
	data[6] = col0[2];

	data[1] = col1[0];
	data[4] = col1[1];
	data[7] = col1[2];

	data[2] = col2[0];
	data[5] = col2[1];
	data[8] = col2[2];
}

template <typename T> inline Matrix3<T> &Matrix3<T>::set(T d) {
	for (int i = 0; i < 9; ++i) data[i] = d;
	return (*this);
}

template <typename T> inline Matrix3<T> &Matrix3<T>::set(T m00, T m10, T m20, T m01, T m11, T m21, T m02, T m12, T m22) {
	data[0] = m00; // column 1 
	data[3] = m10;
	data[6] = m20;

	data[1] = m01; // column 2 
	data[4] = m11;
	data[7] = m21;

	data[2] = m02; // column 3 
	data[5] = m12;
	data[8] = m22;
	return (*this);
}

template <typename T> inline Matrix3<T> &Matrix3<T>::set(const Matrix3<T> &d) {
	memcpy(data, d.data, 9 * sizeof(T));
	return (*this);
}

template <typename T> inline Matrix3<T> &Matrix3<T>::set(const T *d) {
	memcpy(data, d, 9 * sizeof(T));
	return (*this);
}

template <typename T> inline Matrix3<T>& Matrix3<T>::zero() {
	memset(data, 0, 9 * sizeof(T));
	return (*this);
}

template <typename T> inline Matrix3<T> &Matrix3<T>::operator=(T d) {
	return set(d);
}

template <typename T> inline Matrix3<T> &Matrix3<T>::operator=(const Matrix3<T> &d) {
	return set(d);
}

template <typename T> inline Matrix3<T> &Matrix3<T>::operator=(const T *d) {
	return set(d);
}

template <typename T> inline Matrix3<T> &Matrix3<T>::operator+=(const Matrix3<T> &rhs) {
	data[0] += rhs.data[0];
	data[1] += rhs.data[1];
	data[2] += rhs.data[2];
	data[3] += rhs.data[3];
	data[4] += rhs.data[4];
	data[5] += rhs.data[5];
	data[6] += rhs.data[6];
	data[7] += rhs.data[7];
	data[8] += rhs.data[8];

	return (*this);
}

template <typename T> inline Matrix3<T> &Matrix3<T>::operator-=(const Matrix3<T> &rhs) {
	data[0] -= rhs.data[0];
	data[1] -= rhs.data[1];
	data[2] -= rhs.data[2];
	data[3] -= rhs.data[3];
	data[4] -= rhs.data[4];
	data[5] -= rhs.data[5];
	data[6] -= rhs.data[6];
	data[7] -= rhs.data[7];
	data[8] -= rhs.data[8];

	return (*this);
}

template <typename T> inline Matrix3<T> &Matrix3<T>::operator*=(const Matrix3<T> &rhs) {
	Matrix3<T> lhs(data);

	data[0] = lhs.data[0] * rhs.data[0] + lhs.data[1] * rhs.data[3] + lhs.data[2] * rhs.data[6];
	data[3] = lhs.data[3] * rhs.data[0] + lhs.data[4] * rhs.data[3] + lhs.data[5] * rhs.data[6];
	data[6] = lhs.data[6] * rhs.data[0] + lhs.data[7] * rhs.data[3] + lhs.data[8] * rhs.data[6];
	data[1] = lhs.data[0] * rhs.data[1] + lhs.data[1] * rhs.data[4] + lhs.data[2] * rhs.data[7];
	data[4] = lhs.data[3] * rhs.data[1] + lhs.data[4] * rhs.data[4] + lhs.data[5] * rhs.data[7];
	data[7] = lhs.data[6] * rhs.data[1] + lhs.data[7] * rhs.data[4] + lhs.data[8] * rhs.data[7];
	data[2] = lhs.data[0] * rhs.data[2] + lhs.data[1] * rhs.data[5] + lhs.data[2] * rhs.data[8];
	data[5] = lhs.data[3] * rhs.data[2] + lhs.data[4] * rhs.data[5] + lhs.data[5] * rhs.data[8];
	data[8] = lhs.data[6] * rhs.data[2] + lhs.data[7] * rhs.data[5] + lhs.data[8] * rhs.data[8];

	return (*this);
}

template <typename T> inline Matrix3<T> &Matrix3<T>::operator*=(T d) {
	data[0] *= d;
	data[1] *= d;
	data[2] *= d;
	data[3] *= d;
	data[4] *= d;
	data[5] *= d;
	data[6] *= d;
	data[7] *= d;
	data[8] *= d;

	return (*this);
}

template <typename T> inline Matrix3<T> &Matrix3<T>::operator/=(T d) {
	data[0] /= d;
	data[1] /= d;
	data[2] /= d;
	data[3] /= d;
	data[4] /= d;
	data[5] /= d;
	data[6] /= d;
	data[7] /= d;
	data[8] /= d;

	return (*this);
}

template <typename T> const inline T&	Matrix3<T>::operator()(unsigned int i, unsigned int j) const {
	return data[i + 3 * j];
}

template <typename T> inline void Matrix3<T>::Identity() {
	memset(data, 0, 9 * sizeof(T));
	data[0] = data[4] = data[8] = (T)1.0;
}

template <typename T> inline bool Matrix3<T>::isAlmostZero() {
	bool bZero = true;
	for (int i = 0; i < 9; ++i) {
		if (abs(data[i]) > EPS) bZero = false;
	}
	return bZero;
}

template <typename T> inline T Matrix3<T>::determinant() {
	return data[0] * (data[4] * data[8] - data[7] * data[5])
		- data[1] * (data[3] * data[8] - data[5] * data[6])
		+ data[2] * (data[3] * data[7] - data[4] * data[6]);
}

template <typename T> inline T Matrix3<T>::inverse(Matrix3<T>& invM) {
	T det = determinant();
	if (/*det<1.e-9*/det == 0) {
		//cout << "deteminant is zero" << endl;
		memset(&(invM.data[0]), 0, sizeof(T) * 9);
		return det;
	}

	invM.data[0] = (data[4] * data[8] - data[5] * data[7]) / det;
	invM.data[3] = (data[5] * data[6] - data[3] * data[8]) / det;
	invM.data[6] = (data[3] * data[7] - data[4] * data[6]) / det;
	invM.data[1] = (data[2] * data[7] - data[1] * data[8]) / det;
	invM.data[4] = (data[0] * data[8] - data[2] * data[6]) / det;
	invM.data[7] = (data[1] * data[6] - data[0] * data[7]) / det;
	invM.data[2] = (data[1] * data[5] - data[2] * data[4]) / det;
	invM.data[5] = (data[2] * data[3] - data[0] * data[5]) / det;
	invM.data[8] = (data[0] * data[4] - data[1] * data[3]) / det;
	return det;
}

template <typename T> inline Matrix3<T> &Matrix3<T>::transpose() {
	//for (int i=0; i<3 ; i++)
	//	for (int j=i+1; j<3 ; j++)
	//		swap(data[i+j*3],data[i*3+j]);

	std::swap(data[3], data[1]);
	std::swap(data[6], data[2]);
	std::swap(data[5], data[7]);
	return (*this);
}

template <typename T> inline Matrix3<T> Matrix3<T>::transp() const {
	Matrix3<T> mat;
	mat.data[0] = data[0];
	mat.data[1] = data[3];
	mat.data[2] = data[6];
	mat.data[3] = data[1];
	mat.data[4] = data[4];
	mat.data[5] = data[7];
	mat.data[6] = data[2];
	mat.data[7] = data[5];
	mat.data[8] = data[8];
	return mat;
}

template <typename T> inline Matrix3<T> operator+(const Matrix3<T>& lhs, const Matrix3<T>& rhs) {
	return Matrix3<T>(
		lhs.data[0] + rhs.data[0],
		lhs.data[3] + rhs.data[3],
		lhs.data[6] + rhs.data[6],
		lhs.data[1] + rhs.data[1],
		lhs.data[4] + rhs.data[4],
		lhs.data[7] + rhs.data[7],
		lhs.data[2] + rhs.data[2],
		lhs.data[5] + rhs.data[5],
		lhs.data[8] + rhs.data[8]);
}

template <typename T> inline Matrix3<T> operator-(const Matrix3<T>& lhs, const Matrix3<T>& rhs) {
	return Matrix3<T>(
		lhs.data[0] - rhs.data[0],
		lhs.data[3] - rhs.data[3],
		lhs.data[6] - rhs.data[6],
		lhs.data[1] - rhs.data[1],
		lhs.data[4] - rhs.data[4],
		lhs.data[7] - rhs.data[7],
		lhs.data[2] - rhs.data[2],
		lhs.data[5] - rhs.data[5],
		lhs.data[8] - rhs.data[8]);
}

template <typename T> inline Matrix3<T> operator*(const T a, const Matrix3<T>& mat) {
	return Matrix3<T>(
		a*mat.data[0],
		a*mat.data[3],
		a*mat.data[6],
		a*mat.data[1],
		a*mat.data[4],
		a*mat.data[7],
		a*mat.data[2],
		a*mat.data[5],
		a*mat.data[8]);
}

template <typename T> inline Matrix3<T> operator*(const Matrix3<T>& mat, const T a) {
	return Matrix3<T>(
		a*mat.data[0],
		a*mat.data[3],
		a*mat.data[6],
		a*mat.data[1],
		a*mat.data[4],
		a*mat.data[7],
		a*mat.data[2],
		a*mat.data[5],
		a*mat.data[8]);
}

/*template <typename T> inline Vector3<T> operator&( const Matrix3<T>& mat, const Vector3<T>& vec) {
	double xy[2];
	xy[0] = vec[0];
	xy[1] = vec[1];

	double zw[2];
	zw[0] = vec[2];
	zw[1] = 0.0;

	double xx_xy[2];
	xx_xy[0] = mat.data[0];
	xx_xy[1] = mat.data[3];

	double xz_yx[2];
	xz_yx[0] = mat.data[6];
	xz_yx[1] = mat.data[1];

	double yy_yz[2];
	yy_yz[0] = mat.data[4];
	yy_yz[1] = mat.data[7];

	double zx_zy[2];
	zx_zy[0] = mat.data[2];
	zx_zy[1] = mat.data[5];

	double zz_ww[2];
	zz_ww[0] = mat.data[8];
	zz_ww[1] = 0;


	__m128d t1 = _mm_mul_pd( _mm_unpacklo_pd( *reinterpret_cast<__m128d*>(xy), *reinterpret_cast<__m128d*>(xy)),
								_mm_shuffle_pd( *reinterpret_cast<__m128d*>(xx_xy), *reinterpret_cast<__m128d*>(xz_yx), 2) );

	__m128d t2 = _mm_mul_pd( _mm_unpackhi_pd( *reinterpret_cast<__m128d*>(xy), *reinterpret_cast<__m128d*>(xy)),
								_mm_shuffle_pd( *reinterpret_cast<__m128d*>(xx_xy), *reinterpret_cast<__m128d*>(yy_yz), 1) );

	__m128d t3 = _mm_mul_pd( _mm_unpacklo_pd( *reinterpret_cast<__m128d*>(zw), *reinterpret_cast<__m128d*>(zw)),
								_mm_shuffle_pd( *reinterpret_cast<__m128d*>(xz_yx), *reinterpret_cast<__m128d*>(yy_yz), 2) );

	__m128d t4 = _mm_mul_pd( *reinterpret_cast<__m128d*>(zx_zy), *reinterpret_cast<__m128d*>(xy) );
	__m128d t5 = _mm_mul_pd( *reinterpret_cast<__m128d*>(zz_ww), *reinterpret_cast<__m128d*>(zw) );
	__m128d t6 = _mm_unpackhi_pd( t4, t4);;

	__m128d r1 = _mm_add_pd( _mm_add_pd(t1, t2), t3);
	__m128d r2 = _mm_add_sd( _mm_add_sd(t4, t6), t5);

	*reinterpret_cast<__m128d*>(xy) = r1;
	*reinterpret_cast<__m128d*>(zw) = r2;

	return Vector3<T>  (xy[0],xy[1],zw[0] );
}*/

template <typename T> inline Vector3<T> operator*(const Matrix3<T>& mat, const Vector3<T>& vec) {
	return Vector3<T>
		(mat.data[0] * vec[0] + mat.data[1] * vec[1] + mat.data[2] * vec[2],
		mat.data[3] * vec[0] + mat.data[4] * vec[1] + mat.data[5] * vec[2],
		mat.data[6] * vec[0] + mat.data[7] * vec[1] + mat.data[8] * vec[2]);
}

template <typename T> inline void mul_add(const Matrix3<T>& mat, const Vector3<T>& vec, Vector3<T>& ret) {
	ret[0] += mat.data[0] * vec[0] + mat.data[1] * vec[1] + mat.data[2] * vec[2];
	ret[1] += mat.data[3] * vec[0] + mat.data[4] * vec[1] + mat.data[5] * vec[2];
	ret[2] += mat.data[6] * vec[0] + mat.data[7] * vec[1] + mat.data[8] * vec[2];
}


template <typename T> inline Matrix3<T> orthoProjector(const Vector3<T>& vec) {
	Matrix3<T> orthoP(vec[0] * vec[0], vec[0] * vec[1], vec[0] * vec[2], vec[0] * vec[1], vec[1] * vec[1], vec[1] * vec[2], vec[0] * vec[2], vec[1] * vec[2], vec[2] * vec[2]);
	T norm = (T)(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
	return norm != 0.0 ? orthoP * (1.0 / norm) : orthoP;
	//return orthoP ;
}

template <typename T> inline Matrix3<T> outerProduct(const Vector3<T>& vec1, const Vector3<T>& vec2) {
	Matrix3<T> outerP(vec1[0] * vec2[0], vec1[1] * vec2[0], vec1[2] * vec2[0], vec1[0] * vec2[1], vec1[1] * vec2[1], vec1[2] * vec2[1], vec1[0] * vec2[2], vec1[1] * vec2[2], vec1[2] * vec2[2]);
	return outerP;
}

//------------------------------------------------------------------------------------
template <typename T> inline void mul(const std::vector<Matrix3<T> > &m, std::vector<Vector3<T> > &r) {
	const int N = (int)r.size();
#pragma omp parallel for
	for (int i = 0; i < N; i++) {
		r[i] = m[i] * r[i];
	}
}

template <typename T> inline void mul(const std::vector<Matrix3<T> > &m, const std::vector<Vector3<T> > &target, std::vector<Vector3<T> > &r) {
	const int N = (int)target.size();
#pragma omp parallel for
	for (int i = 0; i < N; i++) {
		r[i] = m[i] * target[i];
	}
}

template <typename T> inline void mul(const std::vector<Matrix3<T> > &m1, const std::vector<Matrix3<T> > &m2, const std::vector<Vector3<T> > &target, std::vector<Vector3<T> > &r) {
	const int N = (int)target.size();
#pragma omp parallel for
	for (int i = 0; i < N; i++) {
		r[i] = m1[i] * (m2[i] * target[i]);
	}
}

template <typename T> inline void scale_add_mul(const T s, const std::vector<Vector3<T> > &target, const std::vector<Vector3<T> > &a, const std::vector<Matrix3<T> >& m, std::vector<Vector3<T> > &r) {
	const int N = (int)target.size();
#pragma omp parallel for
	for (int i = 0; i < N; i++) {
		r[i] = m[i] * (s * target[i] + a[i]);
	}
}