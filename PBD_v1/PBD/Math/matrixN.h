#ifndef _MATRIXN_H_
#define _MATRIXN_H_

#include "VectorN.h"
#include <stdlib.h>
#include <iostream>
#include <fstream>

using namespace std;

#ifdef DEBUG
#define DB_CHECK( C ) { if ( ! (C) ) { abort(); } }
#else
#define DB_CHECK( C ) { }
#endif

//----------------------------------------------------------------------------
//	MatrixN class definition.
//!	NxM matrix class.
template <typename T> class MatrixN {
public:
	//------------------------------------------------------------------------
	// Constructor/Destructor
	inline MatrixN() { nRows = nCols = 0;	val = NULL; }
	inline MatrixN(const int n, const int m);
	inline MatrixN(const int n, const int m, const T **a);
	inline MatrixN(const MatrixN<T> &src);
	inline MatrixN(const int n, const int m, const VectorN<T>* v);
	inline ~MatrixN() { if (nCols > 0 || nRows > 0) { delete[] val; nRows = nCols = 0; val = NULL; } }

	//------------------------------------------------------------------------
	// Size
	inline int dimRow()	const { return nRows; }
	inline int dimCol()	const { return nCols; }

	//------------------------------------------------------------------------
	// Simple Setter
	inline int i2i(const int row, const int col) const { return (row * nCols + col); }

	inline void zero();
	inline void identity();
	inline void transpose();

	inline int isZero() const;
	inline int isIdentity() const;

	//----------------------------------------------------------------------------
	// Index Operator
	inline T &operator() (unsigned int i, unsigned int j) { return val[i*nCols + j]; }

	//------------------------------------------------------------------------
	// Setter methods
	inline MatrixN &setDim(const int n, const int m);
	inline MatrixN &set(const int i, const int j, const T d);
	inline MatrixN &set(const MatrixN<T> &src);
	inline MatrixN &set(const T **a);
	inline MatrixN &operator=(const MatrixN &src);
	inline MatrixN &operator=(const T **a);
	inline void		setBlockMatrix3by3(int iBlock, int jBlock, const Matrix3<T> mat);

	//------------------------------------------------------------------------
	// in-place Arithmetic Operator
	inline MatrixN		&operator+=(const MatrixN& a);
	inline MatrixN		&operator-=(const MatrixN& a);
	inline MatrixN		&operator*=(const T d);
	inline MatrixN		&operator/=(const T d);
	inline VectorN<T>	 operator*(const VectorN<T>& a);
	inline MatrixN		 operator*(const T d);

	//------------------------------------------------------------------------
	// Row operation
	inline MatrixN &scaleRow(const int row, const T c);
	inline MatrixN &addRow(const int row0, const int row1, const T c);
	inline MatrixN &interchangeRow(const int row0, const int row1);

	//------------------------------------------------------------------------
	// Getter methods
	inline int isSquareMatrix() const { return (nRows == nCols); }
	inline T get(const int i, const int j) const { return val[i2i(i, j)]; }
	inline T get(const int i) const { return val[i]; }
	T *getMatrixPointer() const { return val; }

	//------------------------------------------------------------------------
	// IO
	void printMatrix() const;
	void printMatrix(char *filename) const;

public:
	int nRows, nCols;
	T *val;
};
						   

// Constructor/Destructor
template <typename T> inline MatrixN<T>::MatrixN(const int n, const int m) {
	nRows = n;
	nCols = m;

	int size = nRows * nCols;

	//cout << nRows <<"      " << nCols << "   " << size << endl;
	val = new T[size];
	for (int i = 0; i < size; i++)
		val[i] = (T)0.0;
}

template <typename T> inline MatrixN<T>::MatrixN(const int n, const int m, const T **a) {
	nRows = n;
	nCols = m;

	val = new T[nRows * nCols];
	for (int i = 0; i < nRows; i++)
		for (int j = 0; j < nCols; j++)
			val[i2i(i, j)] = a[i][j];
}

template <typename T> inline MatrixN<T>::MatrixN(const MatrixN<T> &src) {
	nRows = src.nRows;
	nCols = src.nCols;

	int size = nRows * nCols;
	val = new T[size];
	for (int i = 0; i < size; i++)
		val[i] = src.val[i];
}

template <typename T> inline MatrixN<T>::MatrixN(const int n, const int m, const VectorN<T>* v) {
	//DB_CHECK(n == v[0].dim());

	nRows = n;
	nCols = m;

	val = new T[nRows * nCols];
	for (int i = 0; i < nRows; i++)
		for (int j = 0; j < nCols; j++)
			val[i2i(i, j)] = v[j][i];
}

// Simple Setter
template <typename T> inline void MatrixN<T>::zero() {
	int size = nRows * nCols;
	for (int i = 0; i < size; i++)
		val[i] = 0;
}

template <typename T> inline void MatrixN<T>::identity() {
	//DB_CHECK(nRows == nCols);
	zero();
	for (int i = 0; i < nRows; i++)
		val[i2i(i, i)] = (T)1;
}

template <typename T> inline void MatrixN<T>::transpose() {
	MatrixN<T> tempM(*this);

	for (int i = 0; i < nRows; i++) {
		for (int j = 0; j < nCols; j++) {
			val[i2i(i, j)] = tempM.val[i2i(j, i)];
		}
	}
}

template <typename T> inline int MatrixN<T>::isZero() const {
	int size = nRows * nCols;
	for (int i = 0; i < size; i++)
		if (!IsAlmostZero(val[i]))
			return false;
	return true;
}

template <typename T> inline int MatrixN<T>::isIdentity() const {
	//DB_CHECK(nRows == nCols);

	for (int i = 0; i < nRows; i++) {
		for (int j = 0; j < nCols; j++) {
			if (i == j && !IsAlmostZero(val[i2i(i, j)] - 1))
				return false;
			if (i != j && !IsAlmostZero(val[i2i(i, j)]))
				return false;
		}
	}
	return true;
}
	
// Setter methods
template <typename T> inline MatrixN<T> &MatrixN<T>::setDim(const int n, const int m) {
	if (val != NULL)
		delete[]val;

	nRows = n;
	nCols = m;

	int size = nRows * nCols;
	val = new T[size];
	for (int i = 0; i < size; i++)
		val[i] = (T)0.0;

	return (*this);
}

template <typename T> inline MatrixN<T> &MatrixN<T>::set(const int i, const int j, T d) {
	//DB_CHECK( ( i < nRows && j < nCols ) );
	val[i2i(i, j)] = d;

	return (*this);
}

template <typename T> inline MatrixN<T> &MatrixN<T>::set(const MatrixN<T> &src) {
	//DB_CHECK( ( nRows == src.nRows && nCols == src.nCols  ) );

	int size = nRows * nCols;
	for (int i = 0; i < size; i++)
		val[i] = src.val[i];

	return (*this);
}

template <typename T> inline MatrixN<T> &MatrixN<T>::set(const T **a) {
	val = new T[nRows * nCols];
	for (int i = 0; i < nRows; i++)
		for (int j = 0; j < nCols; j++)
			val[i2i(i, j)] = a[i][j];

	return (*this);
}

template <typename T> inline MatrixN<T> &MatrixN<T>::operator=(const MatrixN &src) {
	set(src);
	return (*this);
}

template <typename T> inline MatrixN<T> &MatrixN<T>::operator=(const T **a) {
	set(a);
	return (*this);
}

template <typename T> inline void		MatrixN<T>::setBlockMatrix3by3(int iBlock, int jBlock, const Matrix3<T> mat) {
	val[i2i(3 * iBlock + 0, 3 * jBlock + 0)] = mat(0, 0);
	val[i2i(3 * iBlock + 1, 3 * jBlock + 0)] = mat(1, 0);
	val[i2i(3 * iBlock + 2, 3 * jBlock + 0)] = mat(2, 0);

	val[i2i(3 * iBlock + 0, 3 * jBlock + 1)] = mat(0, 1);
	val[i2i(3 * iBlock + 1, 3 * jBlock + 1)] = mat(1, 1);
	val[i2i(3 * iBlock + 2, 3 * jBlock + 1)] = mat(2, 1);

	val[i2i(3 * iBlock + 0, 3 * jBlock + 2)] = mat(0, 2);
	val[i2i(3 * iBlock + 1, 3 * jBlock + 2)] = mat(1, 2);
	val[i2i(3 * iBlock + 2, 3 * jBlock + 2)] = mat(2, 2);
}

//-------------------------------------------------------------------
// in-place Arithmetic Operator
template <typename T> inline MatrixN<T> &MatrixN<T>::operator+=(const MatrixN& a) {
	//DB_CHECK( ( nRows == a.nRows && nCols == a.nCols ) );

	int size = nRows * nCols;
	for (int i = 0; i < size; i++)
		val[i] += a.val[i];

	return (*this);
}

template <typename T> inline MatrixN<T> &MatrixN<T>::operator-=(const MatrixN& a) {
	//DB_CHECK( ( nRows == a.nRows && nCols == a.nCols ) );

	int size = nRows * nCols;
	for (int i = 0; i < size; i++)
		val[i] -= a.val[i];

	return (*this);
}

template <typename T> inline MatrixN<T> &MatrixN<T>::operator*=(const T a) {
	//DB_CHECK( ( nRows == a.nRows && nCols == a.nCols ) );

	int size = nRows * nCols;
	for (int i = 0; i < size; i++)
		val[i] *= a;

	return (*this);
}

template <typename T> inline MatrixN<T> &MatrixN<T>::operator/=(const T a) {
	//DB_CHECK( ( nRows == a.nRows && nCols == a.nCols ) );

	int size = nRows * nCols;
	for (int i = 0; i < size; i++)
		val[i] /= a;

	return (*this);
}

template <typename T> inline VectorN<T> MatrixN<T>::operator*(const VectorN<T>& a) {
	//DB_CHECK( ( nRows == a.nRows && nCols == a.nCols ) );

	VectorN<T> b(nRows);

	for (int i = 0; i < nCols; i++) {
		for (int j = 0; j < nRows; j++) {
			b[j] += (val[i2i(j, i)] * a[i]);
		}
	}

	return b;
}

template <typename T> inline MatrixN<T> MatrixN<T>::operator*(const T d) {
	MatrixN<T> A(nRows, nCols);
	int size = nRows * nCols;
	for (int i = 0; i < size; i++) {
		A.val[i] = val[i] * d;
	}
	return A;
}

//-------------------------------------------------------------------
// Row Operation
template <typename T> inline MatrixN<T> &MatrixN<T>::scaleRow(const int row, const T c) {
	for (int i = 0; i < nCols; i++)
		val[i2i(row, i)] *= c;
	return (*this);
}

template <typename T> inline MatrixN<T> &MatrixN<T>::addRow(const int row0, const int row1, const T c) {
	for (int i = 0; i < nCols; i++)
		val[i2i(row0, i)] += val[i2i(row1, i)] * c;
	return (*this);
}
	
template <typename T> inline MatrixN<T> &MatrixN<T>::interchangeRow(const int row0, const int row1) {
	T temp;
	for (int i = 0; i < nCols; i++) {
		temp = val[i2i(row0, i)];
		val[i2i(row0, i)] = val[i2i(row1, i)];
		val[i2i(row1, i)] = temp;
	}
	return (*this);
}

//-------------------------------------------------------------------
// IO
template <typename T> inline void MatrixN<T>::printMatrix() const {
	for (int i = 0; i < nRows; i++) {
		for (int j = 0; j < nCols; j++) {
			int tmp1 = (int)(val[i2i(i, j)] * 1000);
			float tmp2 = tmp1 / 1000.0f;

			cout << tmp2 << " ";

			//cout << val[i2i(i,j)] << " ";
			//printf("%lf ",val[i2i(i,j)]);	
		}
		cout << endl;
	}
	cout << endl;
}

template <typename T> inline void MatrixN<T>::printMatrix(char *filename) const {
	FILE *fp = fopen(filename, "w");

	//fprintf(fp, "[ ");
	for (int i = 0; i < nRows; i++) {
		for (int j = 0; j < nCols; j++) {
			if (val[i2i(i, j)] >= 0 && val[i2i(i, j)] < 10)
				fprintf(fp, "%0.3lf ", val[i2i(i, j)]);
			else
				fprintf(fp, "%0.2lf ", val[i2i(i, j)]);

		}
		if (i != nRows - 1)
			fprintf(fp, " ;\n");
		else {}
		//fprintf(fp, " ]\n");
	}
	fclose(fp);
}

//-------------------------------------------------------------------
// Inline implementation of functions related to MatrixN
template <typename T> inline void neg(MatrixN<T> &a, const MatrixN<T> &b) {
	//DB_CHECK( ( a.dimRow() == b.dimRow() && a.dimCol() == b.dimCol() ) );

	int n = a.dimRow();
	int m = a.dimCol();

	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++)
			a.set(i, j, -b.get(i, j));
}

template <typename T> inline void add(MatrixN<T> &a, const MatrixN<T> &b, const MatrixN<T> &c) {
	//DB_CHECK( ( a.dimRow() == b.dimRow() && a.dimCol() == b.dimCol() )
	//	   && ( b.dimRow() == c.dimRow() && b.dimCol() == c.dimCol() ));

	int n = a.dimRow();
	int m = a.dimCol();

	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++)
			a.set(i, j, b.get(i, j) + c.get(i, j));
}

template <typename T> inline void dif(MatrixN<T> &a, const MatrixN<T> &b, const MatrixN<T> &c) {
	//DB_CHECK( ( a.dimRow() == b.dimRow() && a.dimCol() == b.dimCol() )
	//	   && ( b.dimRow() == c.dimRow() && b.dimCol() == c.dimCol() ));

	int n = a.dimRow();
	int m = a.dimCol();

	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++)
			a.set(i, j, b.get(i, j) - c.get(i, j));
}

template <typename T> inline void mul(MatrixN<T> &a, const MatrixN<T> &b, const T &c) {
	//DB_CHECK( ( a.dimRow() == b.dimRow() && a.dimCol() == b.dimCol() ) );
	int n = a.dimRow();
	int m = a.dimCol();

	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++)
			a.set(i, j, b.get(i, j) * c);
}

template <typename T> inline void mul(MatrixN<T> &a, const T &c, const MatrixN<T> &b) {
	mul(a, b, c);
}

template <typename T> inline void mul(VectorN<T> &a, const MatrixN<T> &b, const VectorN<T> &c) {
	//DB_CHECK( ( a.dim() == b.dimRow() && c.dim() == b.dimCol() ) );

	int n = b.dimRow();
	int m = b.dimCol();

	for (int i = 0; i < n; i++) {
		T val = (T)0;
		for (int j = 0; j < m; j++) {
			val += b.get(i, j) * c.get(j);
		}
		a.set(i, val);
	}
}

template <typename T> inline void mul(MatrixN<T> &a, const MatrixN<T> &b, const MatrixN<T> &c) {
	//DB_CHECK( ( a.dimRow() == b.dimRow() && a.dimCol() == c.dimCol() )
	//	   && ( b.dimCol() == c.dimRow() ) );

	int n = a.dimRow();
	int m = a.dimCol();
	int l = b.dimCol();

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			T val = (T)0;

			for (int k = 0; k < l; k++)
				val += (b.get(i, k) * c.get(k, j));
			a.set(i, j, val);
		}
	}
}

template <typename T> inline void transpose(MatrixN<T> &trans, const MatrixN<T> &a) {
	//DB_CHECK( (trans.dimRow() == a.dimCol() && trans.dimCol() == a.dimRow()) );

	int n = a.dimRow();
	int m = a.dimCol();

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			trans.set(j, i, a.get(i, j));
		}
	}
}

template <typename T> std::ostream &operator<<(std::ostream &strm, MatrixN<T> &a) {
	int NRow = a.dimRow();
	int NCol = a.dimCol();
	for (int i = 0; i < NRow; i++) {
		strm << "[ ";
		for (int j = 0; j < NCol - 1; j++) {
			strm << a(i, j) << ", ";
		}
		strm << a(i, NCol - 1) << " ]" << endl;;
	}
	return strm;
}

template <typename T> inline void directMethod(MatrixN<T> A, VectorN<T> b, VectorN<T> &x) {
	//DB_CHECK( (A.dimCol() == A.dimRow()) && (A.dimRow()==b.dim()) && (A.dimRow()==x.dim());
	int N = A.dimRow();
	float e;

	// Gauss-Elimination
	for (int k = 0; k < N - 1; k++) { // k는 현재 소거시키고 있는 Column의 index (0 ~ N-2)
		for (int i = k + 1; i < N; i++) { // i는 현재 소거중인 (k열의) i번째 행의 index (k+1 ~ N-1)

			//A의 i행 i-1열을 소거하자. 근데 i-1행, i-1열의 성분이 0이면 곤란하지
			//if(A(i-1,i-1)==0){
			//	assert(FALSE);
			//}

			e = A(i, k) / A(k, k);
			for (int j = 0; j < N; j++) {
				// i번째 행에서 k행에 e를 곱해서 뺀다. j는 (현재 소거 중인 Row의) Column의 index (i-1 ~ N-1)
				A(i, j) = A(i, j) - e * A(k, j);
			}
			b[i] -= e * b[k];
		}
	}

	//Back-Substitution
	x[N - 1] = b[N - 1] / A(N - 1, N - 1);
	for (int i = N - 2; i >= 0; i--) {
		e = 0.0f;
		for (int j = i + 1; j <= N - 1; j++) {
			e += A(i, j)*x[j];
		}
		x[i] = (b[i] - e) / A(i, i);
	}
}
#endif