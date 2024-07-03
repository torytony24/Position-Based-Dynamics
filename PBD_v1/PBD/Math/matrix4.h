#ifndef __MATRIX4_H__
#define __MATRIX4_H__

#include "Vector.h"

template <typename T> class Matrix4 {
public:
	////////////////////////////////////////////////////////////////////////////
	// Constructors
	Matrix4();
	Matrix4(const Matrix4<T> &d);
	Matrix4(T d00, T d10, T d20, T d30,   // col 0
		T d01, T d11, T d21, T d31,   // col 1
		T d02, T d12, T d22, T d32,   // col 2
		T d03, T d13, T d23, T d33); // col 3
	//Matrix4( T d[16] );
	//Matrix4( Vector4<T> col0, Vector4<T> col1, Vector4<T> col2, Vector4<T> col3 );	

////////////////////////////////////////////////////////////////////////////
	// Index operators
	inline T& operator()(unsigned int i, unsigned int j);

	inline			T&	operator[](unsigned int i);
	inline const	T&  operator[](unsigned int i) const;

	//////////////////////////////////////////////////////////////////////////
	// Setter and Assignment 
	inline Matrix4 &set(T d);
	inline Matrix4 &set(const Matrix4 &d);
	inline Matrix4 &set(const T *d);
	inline Matrix4 &set(T d[16]);
	inline Matrix4 &operator=(T d);
	inline Matrix4 &operator=(const Matrix4 &d);
	inline Matrix4 &operator=(const T *d);

	//////////////////////////////////////////////////////////////////////////
	// Getter and Assignment 
	inline float* getFloat() {
		return (float*)data;
	}

	//////////////////////////////////////////////////////////////////////////
	// In place arithmetic
	//inline void				operator+=( T d );
	//inline void				operator-=( T d );
	//inline void				operator*=( T d );
	//inline void				operator/=( T d );

	//inline void				operator+=( const Matrix4 &d );
	//inline void				operator-=( const Matrix4 &d );
	inline void				operator*=(const Matrix4 &d);
	//inline void				operator/=( const Matrix4 &d );

	inline Matrix4			operator*(const Matrix4 &rhs) const;
	inline Vector3<T>		operator*(const Vector3<T> &rhs) const;
	inline Vector4<T>		operator*(const Vector4<T> &rhs) const;

	inline Matrix4			operator+(const Matrix4 &rhs) const;
	inline Matrix4			operator-(const Matrix4 &rhs) const;
	inline Matrix4			operator*(T d) const;
	inline friend Matrix4	operator*(T d, const Matrix4 &mat) {
		return mat * d;
	}

	////////////////////////////////////////////////////////////////////////////
	// Linear algebra methods
	inline Matrix4 inverse();
	inline Matrix4 inverseSE3();
	inline Matrix4 transpose();
	//inline T       determinant();


	////////////////////////////////////////////////////////////////////////////
	// OpenGL methods
	inline void loadIdentity();
	inline void setTranslation(const Vector3<T> &translation);
	inline void setScale(const Vector3<T> &scale);
	inline void setRotation(const Vector3<T> &axis, const T angle);
	inline void translate(const Vector3<T> &translation);
	inline void scale(const Vector3<T> &scale);
	inline void rotate(const Vector3<T> &axis, const T angle);

	/////////////////////////////////////////////////////////////////////////
	// general methods
	inline void print();

public:
	T data[16];
};

template <typename T> inline Matrix4<T> operator*(T d, const Matrix4<T>& mat) {
	return mat * d;
}

template <typename T> inline Matrix4<T>::Matrix4() {
	memset(data, 0, 16 * sizeof(T));
	data[0] = 1.0f;
	data[5] = 1.0f;
	data[10] = 1.0f;
	data[15] = 1.0f;
}

template <typename T> inline Matrix4<T>::Matrix4(const Matrix4 &d) {
	memcpy(data, d.data, 16 * sizeof(T));
}

template <typename T> inline Matrix4<T>::Matrix4(
	T d0, T d1, T d2, T d3,			// col 0
	T d4, T d5, T d6, T d7,			// col 1
	T d8, T d9, T d10, T d11,		// col 2
	T d12, T d13, T d14, T d15)		// col 3
{
	data[0] = d0;
	data[1] = d1;
	data[2] = d2;
	data[3] = d3;
	data[4] = d4;
	data[5] = d5;
	data[6] = d6;
	data[7] = d7;
	data[8] = d8;
	data[9] = d9;
	data[10] = d10;
	data[11] = d11;
	data[12] = d12;
	data[13] = d13;
	data[14] = d14;
	data[15] = d15;
}


template <typename T> inline T& Matrix4<T>::operator()(unsigned int i, unsigned int j) {
	return data[i + 4 * j];
}

template <typename T> inline T &Matrix4<T>::operator[](unsigned int i) {
	return data[i];
}

template <typename T> inline const T&  Matrix4<T>::operator[](unsigned int i) const {
	return data[i];
}

template <typename T> inline Matrix4<T> &Matrix4<T>::set(T d) {
	for (int i = 0; i < 16; ++i) data[i] = d;
	return (*this);
}

template <typename T> inline Matrix4<T> &Matrix4<T>::set(const Matrix4 &d) {
	memcpy(data, d.data, 16 * sizeof(T));
	return (*this);
}

template <typename T> inline Matrix4<T> &Matrix4<T>::set( const T *d ) {
	memcpy( data, d, 16*sizeof(T) );
	return (*this);
}
 
template <typename T> inline Matrix4<T> &Matrix4<T>::set(T d[16]) {
	data[0] = d[0];
	data[1] = d[1];
	data[2] = d[2];
	data[3] = d[3];
	data[4] = d[4];
	data[5] = d[5];
	data[6] = d[6];
	data[7] = d[7];
	data[8] = d[8];
	data[9] = d[9];
	data[10] = d[10];
	data[11] = d[11];
	data[12] = d[12];
	data[13] = d[13];
	data[14] = d[14];
	data[15] = d[15];

	return (*this);
}

template <typename T> inline Matrix4<T> &Matrix4<T>::operator=(T d) {
	return set(d);
}

template <typename T> inline Matrix4<T> &Matrix4<T>::operator=(const Matrix4 &d) {
	return set(d);
}

template <typename T> inline Matrix4<T> &Matrix4<T>::operator=(const T *d) {
	return set(d);
}

template <typename T> inline void Matrix4<T>::loadIdentity() {
	memset(data, 0, 16 * sizeof(T));
	data[0] = 1.0f;
	data[5] = 1.0f;
	data[10] = 1.0f;
	data[15] = 1.0f;
}

template <typename T> inline void Matrix4<T>::setTranslation(const Vector3<T> &translation) {
	loadIdentity();

	data[12] = translation[0];
	data[13] = translation[1];
	data[14] = translation[2];
}

template <typename T> inline void Matrix4<T>::setScale(const Vector3<T> &scale) {
	loadIdentity();

	data[0] = scale[0];
	data[5] = scale[1];
	data[10] = scale[2];
}

template <typename T> inline void Matrix4<T>::setRotation(const Vector3<T> &axis, const T angle) {
	loadIdentity();

	Vector3<T> u(axis);
	normalize(u);

	T sinAngle = (T)sin(angle);
	T cosAngle = (T)cos(angle);
	T oneMinusCosAngle = (T)(1.0 - cosAngle);

	data[0] = (u[0])*(u[0]) + cosAngle * (1 - (u[0])*(u[0]));
	data[4] = (u[0])*(u[1])*(oneMinusCosAngle)-sinAngle * u[2];
	data[8] = (u[0])*(u[2])*(oneMinusCosAngle)+sinAngle * u[1];

	data[1] = (u[0])*(u[1])*(oneMinusCosAngle)+sinAngle * u[2];
	data[5] = (u[1])*(u[1]) + cosAngle * (1 - (u[1])*(u[1]));
	data[9] = (u[1])*(u[2])*(oneMinusCosAngle)-sinAngle * u[0];

	data[2] = (u[0])*(u[2])*(oneMinusCosAngle)-sinAngle * u[1];
	data[6] = (u[1])*(u[2])*(oneMinusCosAngle)+sinAngle * u[0];
	data[10] = (u[2])*(u[2]) + cosAngle * (1 - (u[2])*(u[2]));
}

template <typename T> inline void Matrix4<T>::translate(const Vector3<T> &translation) {
	data[12] += translation[0];
	data[13] += translation[1];
	data[14] += translation[2];
}

template <typename T> inline void Matrix4<T>::scale(const Vector3<T> &scale) {
	data[0] *= scale[0];
	data[5] *= scale[1];
	data[10] *= scale[2];
}

template <typename T> inline void Matrix4<T>::rotate(const Vector3<T> &axis, const T angle) {
	Matrix4<T> r;

	Vector3<T> u(axis);
	normalize(u);

	T sinAngle = (T)sin(angle);
	T cosAngle = (T)cos(angle);
	T oneMinusCosAngle = (T)(1.0 - cosAngle);

	r.data[0] = (u[0])*(u[0]) + cosAngle * (1 - (u[0])*(u[0]));
	r.data[4] = (u[0])*(u[1])*(oneMinusCosAngle)-sinAngle * u[2];
	r.data[8] = (u[0])*(u[2])*(oneMinusCosAngle)+sinAngle * u[1];

	r.data[1] = (u[0])*(u[1])*(oneMinusCosAngle)+sinAngle * u[2];
	r.data[5] = (u[1])*(u[1]) + cosAngle * (1 - (u[1])*(u[1]));
	r.data[9] = (u[1])*(u[2])*(oneMinusCosAngle)-sinAngle * u[0];

	r.data[2] = (u[0])*(u[2])*(oneMinusCosAngle)-sinAngle * u[1];
	r.data[6] = (u[1])*(u[2])*(oneMinusCosAngle)+sinAngle * u[0];
	r.data[10] = (u[2])*(u[2]) + cosAngle * (1 - (u[2])*(u[2]));

	(*this) = r * (*this);
}

template <typename T> inline Matrix4<T> Matrix4<T>::inverse() {
	Matrix4<T> rtn;

	float matr[4][4], ident[4][4];
	int i, j, k, l, ll;
	int icol = 0, irow = 0;
	int indxc[4], indxr[4], ipiv[4];
	float big, dum, pivinv;

	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			matr[i][j] = (float)this->data[i][j];
			ident[i][j] = 0.0f;
		}
		ident[i][i] = 1.0f;
	}
	// Gauss-Jordan elimination with full pivoting.  Yes, folks, a 
	// GL Matrix4 is inverted like any other, since the identity is		
	// still the identity.

	// by numerical recipe in C second edition, pg 39

	for (j = 0; j <= 3; j++) ipiv[j] = 0;
	for (i = 0; i <= 3; i++) {
		big = 0.0;
		for (j = 0; j <= 3; j++) {
			if (ipiv[j] != 1) {
				for (k = 0; k <= 3; k++) {
					if (ipiv[k] == 0) {
						if (fabs(matr[j][k]) >= big) {
							big = (float)fabs(matr[j][k]);
							irow = j;
							icol = k;
						}
					}
					else if (ipiv[k] > 1) {
						;
					}
				}
			}
		}
		++(ipiv[icol]);
		if (irow != icol) {
			for (l = 0; l <= 3; l++) swap(matr[irow][l], matr[icol][l]);
			for (l = 0; l <= 3; l++) swap(ident[irow][l], ident[icol][l]);
		}
		indxr[i] = irow;
		indxc[i] = icol;
		if (matr[icol][icol] == 0.0f) { ; }
		pivinv = 1.0f / matr[icol][icol];
		matr[icol][icol] = 1.0f;
		for (l = 0; l <= 3; l++) matr[icol][l] *= pivinv;
		for (l = 0; l <= 3; l++) ident[icol][l] *= pivinv;
		for (ll = 0; ll <= 3; ll++) {
			if (ll != icol) {
				dum = matr[ll][icol];
				matr[ll][icol] = 0.0f;
				for (l = 0; l <= 3; l++) matr[ll][l] -= matr[icol][l] * dum;
				for (l = 0; l <= 3; l++) ident[ll][l] -= ident[icol][l] * dum;
			}
		}
	}
	for (l = 3; l >= 0; l--) {
		if (indxr[l] != indxc[l]) {
			for (k = 0; k <= 3; k++) {
				swap(matr[k][indxr[l]], matr[k][indxc[l]]);
			}
		}
	}
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			rtn.data[i][j] = matr[i][j];
		}
	}

	return (rtn);
}

template <typename T> inline Matrix4<T> Matrix4<T>::transpose() {
	for (int i = 0; i < 4; i++)
		for (int j = i + 1; j < 4; j++)
			swap(data[i + j * 4], data[i * 4 + j]);
	return (*this);
}

template <typename T> inline Matrix4<T> Matrix4<T>::inverseSE3() {
	Matrix4<T> inv;
	T det, oneOverDet;

	/* Compute 3*3 upper left matrix portion determinant */
	det = data[0] * ((data[5] * data[10]) - (data[6] * data[9]))
		+ data[1] * ((data[6] * data[8]) - (data[4] * data[10]))
		+ data[2] * ((data[4] * data[9]) - (data[5] * data[8]));

	oneOverDet = 1.0f / det;

	/* Compute affine part inverse */
	inv.data[0] = ((data[5] * data[10]) - (data[6] * data[9])) * oneOverDet;
	inv.data[1] = ((data[2] * data[9]) - (data[1] * data[10]))* oneOverDet;
	inv.data[2] = ((data[1] * data[6]) - (data[2] * data[5])) * oneOverDet;
	inv.data[3] = 0.0f;

	inv.data[4] = ((data[6] * data[8]) - (data[4] * data[10]))* oneOverDet;
	inv.data[5] = ((data[0] * data[10]) - (data[2] * data[8])) * oneOverDet;
	inv.data[6] = ((data[2] * data[4]) - (data[0] * data[6])) * oneOverDet;
	inv.data[7] = 0.0f;

	inv.data[8] = ((data[4] * data[9]) - (data[5] * data[8])) * oneOverDet;
	inv.data[9] = ((data[1] * data[8]) - (data[0] * data[9])) * oneOverDet;
	inv.data[10] = ((data[0] * data[5]) - (data[1] * data[4])) * oneOverDet;
	inv.data[11] = 0.0f;

	/* Compute translation part inverse*/
	inv.data[12] = -((data[12] * inv.data[0]) + (data[13] * inv.data[4]) + (data[14] * inv.data[8]));
	inv.data[13] = -((data[12] * inv.data[1]) + (data[13] * inv.data[5]) + (data[14] * inv.data[9]));
	inv.data[14] = -((data[12] * inv.data[2]) + (data[13] * inv.data[6]) + (data[14] * inv.data[10]));
	inv.data[15] = 1.0f;

	return inv;
}

template <typename T> inline void Matrix4<T>::operator*=(const Matrix4<T> &d) {
	(*this) = (*this)*d;
}

template <typename T> inline Matrix4<T> Matrix4<T>::operator*(const Matrix4<T> &rhs) const {
	return Matrix4<T>(
		data[0] * rhs.data[0] + data[4] * rhs.data[1] + data[8] * rhs.data[2] + data[12] * rhs.data[3],
		data[1] * rhs.data[0] + data[5] * rhs.data[1] + data[9] * rhs.data[2] + data[13] * rhs.data[3],
		data[2] * rhs.data[0] + data[6] * rhs.data[1] + data[10] * rhs.data[2] + data[14] * rhs.data[3],
		data[3] * rhs.data[0] + data[7] * rhs.data[1] + data[11] * rhs.data[2] + data[15] * rhs.data[3],

		data[0] * rhs.data[4] + data[4] * rhs.data[5] + data[8] * rhs.data[6] + data[12] * rhs.data[7],
		data[1] * rhs.data[4] + data[5] * rhs.data[5] + data[9] * rhs.data[6] + data[13] * rhs.data[7],
		data[2] * rhs.data[4] + data[6] * rhs.data[5] + data[10] * rhs.data[6] + data[14] * rhs.data[7],
		data[3] * rhs.data[4] + data[7] * rhs.data[5] + data[11] * rhs.data[6] + data[15] * rhs.data[7],

		data[0] * rhs.data[8] + data[4] * rhs.data[9] + data[8] * rhs.data[10] + data[12] * rhs.data[11],
		data[1] * rhs.data[8] + data[5] * rhs.data[9] + data[9] * rhs.data[10] + data[13] * rhs.data[11],
		data[2] * rhs.data[8] + data[6] * rhs.data[9] + data[10] * rhs.data[10] + data[14] * rhs.data[11],
		data[3] * rhs.data[8] + data[7] * rhs.data[9] + data[11] * rhs.data[10] + data[15] * rhs.data[11],

		data[0] * rhs.data[12] + data[4] * rhs.data[13] + data[8] * rhs.data[14] + data[12] * rhs.data[15],
		data[1] * rhs.data[12] + data[5] * rhs.data[13] + data[9] * rhs.data[14] + data[13] * rhs.data[15],
		data[2] * rhs.data[12] + data[6] * rhs.data[13] + data[10] * rhs.data[14] + data[14] * rhs.data[15],
		data[3] * rhs.data[12] + data[7] * rhs.data[13] + data[11] * rhs.data[14] + data[15] * rhs.data[15]);
}

template <typename T> inline Vector3<T> Matrix4<T>::operator*(const Vector3<T> &rhs) const {
	return Vector3<T>(
		data[0] * rhs[0] + data[4] * rhs[1] + data[8] * rhs[2] + data[12],
		data[1] * rhs[0] + data[5] * rhs[1] + data[9] * rhs[2] + data[13],
		data[2] * rhs[0] + data[6] * rhs[1] + data[10] * rhs[2] + data[14]);
}

template <typename T> inline Vector4<T> Matrix4<T>::operator*(const Vector4<T> &rhs) const {
	return Vector4<T>(
		data[0] * rhs[0] + data[4] * rhs[1] + data[8] * rhs[2] + data[12] * rhs[3],
		data[1] * rhs[0] + data[5] * rhs[1] + data[9] * rhs[2] + data[13] * rhs[3],
		data[2] * rhs[0] + data[6] * rhs[1] + data[10] * rhs[2] + data[14] * rhs[3],
		data[3] * rhs[0] + data[7] * rhs[1] + data[11] * rhs[2] + data[15] * rhs[3]);
}

template <typename T> std::ostream &operator<<(std::ostream &strm, const Matrix4<T> &d) {
	strm << "[" << d(0, 0) << ", " << d(0, 1) << ", " << d(0, 2) << ", " << d.data[12] << "]\n";
	strm << "|" << d(1, 0) << ", " << d(1, 1) << ", " << d(1, 2) << ", " << d.data[13] << "|\n";
	strm << "|" << d(2, 0) << ", " << d(2, 1) << ", " << d(2, 2) << ", " << d.data[14] << "|\n";
	strm << "[" << d(3, 0) << ", " << d(3, 1) << ", " << d(3, 2) << ", " << d.data[15] << "]\n";
	return strm;
}

template <typename T> inline Matrix4<T> Matrix4<T>::operator+(const Matrix4<T> &rhs) const {
	return Matrix4<T>(
		data[0] + rhs.data[0], data[1] + rhs.data[1], data[2] + rhs.data[2], data[3] + rhs.data[3],
		data[4] + rhs.data[4], data[5] + rhs.data[5], data[6] + rhs.data[6], data[7] + rhs.data[7],
		data[8] + rhs.data[8], data[9] + rhs.data[9], data[10] + rhs.data[10], data[11] + rhs.data[11],
		data[12] + rhs.data[12], data[13] + rhs.data[13], data[14] + rhs.data[14], data[15] + rhs.data[15]);
}

template <typename T> inline Matrix4<T> Matrix4<T>::operator-(const Matrix4<T> &rhs) const {
	return Matrix4<T>(
		data[0] - rhs.data[0], data[1] - rhs.data[1], data[2] - rhs.data[2], data[3] - rhs.data[3],
		data[4] - rhs.data[4], data[5] - rhs.data[5], data[6] - rhs.data[6], data[7] - rhs.data[7],
		data[8] - rhs.data[8], data[9] - rhs.data[9], data[10] - rhs.data[10], data[11] - rhs.data[11],
		data[12] - rhs.data[12], data[13] - rhs.data[13], data[14] - rhs.data[14], data[15] - rhs.data[15]);
}

template <typename T> inline Matrix4<T> Matrix4<T>::operator*(T d) const {
	return Matrix4<T>(
		data[0] * d, data[1] * d, data[2] * d, data[3] * d,
		data[4] * d, data[5] * d, data[6] * d, data[7] * d,
		data[8] * d, data[9] * d, data[10] * d, data[11] * d,
		data[12] * d, data[13] * d, data[14] * d, data[15] * d);
}

template <typename T> inline void Matrix4<T>::print() {
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++)
			cout << data[i + j * 4] << " ";
		cout << endl;
	}
}

#endif // __MATRIX4_H__
