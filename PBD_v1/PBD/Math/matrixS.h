#ifndef _MATRIXS_H_
#define _MATRIXS_H_

#include "VectorN.h"
#include "matrixN.h"
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

//============================================================================
// MatrixS class declaration (Sparse matrix MxN, Compressed Sparse Row(CSR))
//----------------------------------------------------------------------------
// [1 0 0 0] 
// [2 3 0 0] 
// [0 4 5 0]
// [6 0 7 8] 
// int ri[5] = {0,1,3,5,8};
// int ci[8] = {0, 0,1, 1,2, 0,2,3};
// T  val[8] = {1, 2,3, 4,5, 6,7,8};


template <typename T> class MatrixS {
public:
	inline MatrixS();
	inline MatrixS(int n, int nonZeroN);
	inline MatrixS(const MatrixS<T>& A);
	inline MatrixS(const MatrixN<T>& A);
	inline ~MatrixS();

	//------------------------------------------------------------------------
	// Setter methods
	inline void setMatrixN(const MatrixN<T>& A);
	inline void setSparsity(const MatrixS<T>& A);
	template<typename U>
	inline void setSparsity(const MatrixS<U>& A);
	inline void setDim(int n, int nnz);
	inline void zero();
	inline void setValue(int i, int j, T newVal);
	inline void setValue(int i, bool isRow, T newVal);
	inline void setValue(const MatrixS<T> A);
	inline void setValueInc(int i, int j, T IncVal);
	inline void setValueMul(int i, int j, T MulVal);
	inline void setBlockMatrix3by3(int iBlock, int jBlock, const Matrix3<T> mat);


	//------------------------------------------------------------------------
	// Getter methods
	inline int dim() const { return n; }
	inline int nnz() const { return numNonZeros; };

	inline T& get(int i, int j) const;

	template <typename S> friend inline void mul(MatrixS<S> &A, VectorN<S> &x, VectorN<S> &b);
	template <typename S> friend inline void mul(const MatrixS< Matrix3<S> > &A, const std::vector< Vector3<S> > &x, std::vector< Vector3<S> > &b);
	template <typename S> friend std::ostream &operator<<(std::ostream &strm, MatrixN<S> &m);

	template <typename S> friend int cg(MatrixS<S> &A, VectorN<S> &x, VectorN<S> &b, int maxItrNums, double tol);
	template <typename S> friend int iccg(MatrixS<S> &A, VectorN<S> &x, VectorN<S> &b, int maxItrNums, double tol);

	void printMatrix() const;
	void printMatrix(char *filename) const;

public:
	int		n;				// number of row(column)
	int		numNonZeros;	// number of non zeros

	int		*ri;			// first nonzero index
	int		*ci;
	T		*val;
};

template <typename T> inline MatrixS<T>::MatrixS() {
	n = 0;
	numNonZeros = 0;
	ri = NULL;
	ci = NULL;
	val = NULL;
}

template <typename T> inline MatrixS<T>::~MatrixS() {
	if (n > 0) {
		delete[] ri;
		n = 0;
		ri = NULL;
	}
	if (numNonZeros > 0) {
		delete[] ci;
		delete[] val;
		numNonZeros = 0;
		ci = NULL;
		ri = NULL;
	}
}

template <typename T> inline MatrixS<T>::MatrixS(int __n, int nnz) {
	n = __n;
	numNonZeros = nnz;

	ri = new int[n + 1];
	val = new T[numNonZeros];
	ci = new int[numNonZeros];

	for (int i = 0; i < n + 1; i++)
		ri[i] = 0;
	for (int i = 0; i < nnz; i++) {
		val[i] = (T)0;
		ci[i] = 0;
	}
}

template <typename T> inline MatrixS<T>::MatrixS(const MatrixS<T>& A) {
	n = A.n;
	numNonZeros = A.numNonZeros;

	ri = new int[n + 1];
	val = new T[numNonZeros];
	ci = new int[numNonZeros];

	for (int i = 0; i < n + 1; i++)
		ri[i] = A.ri[i];
	for (int i = 0; i < numNonZeros; i++) {
		val[i] = A.val[i];
		ci[i] = A.ci[i];
	}
}

template <typename T> inline MatrixS<T>::MatrixS(const MatrixN<T>& A) {
	int row, col, size, nnz;

	row = A.dimRow();
	col = A.dimCol();
	size = row * col;
	n = row;
	nnz = 0;
	for (int i = 0; i < size; i++) {
		if (A.val[i] > DC_EPS || A.val[i] < (-1.0f)*DC_EPS) {
			nnz++;
		}
	}
	numNonZeros = nnz;

	ri = new int[n + 1];
	ci = new int[numNonZeros];
	val = new T[numNonZeros];

	nnz = 0;
	for (row = 0; row < n; row++) {
		ri[row] = nnz;
		for (col = 0; col < n; col++) {
			if (A.val[A.i2i(row, col)] > DC_EPS || A.val[A.i2i(row, col)] < (-1.0f)*DC_EPS) {
				ci[nnz] = col;
				val[nnz] = A.val[A.i2i(row, col)];
				nnz++;
			}
		}
	}
	ri[row] = nnz;
}

template <typename T> inline void MatrixS<T>::setMatrixN(const MatrixN<T>& A) {
	int row, col, size, nnz;

	row = A.dimRow();
	col = A.dimCol();
	size = row * col;
	n = row;
	nnz = 0;
	for (int i = 0; i < size; i++) {
		if (A.val[i] > DC_EPS || A.val[i] < (-1.0f)*DC_EPS) {
			nnz++;
		}
	}
	numNonZeros = nnz;

	ri = new int[n + 1];
	ci = new int[numNonZeros];
	val = new T[numNonZeros];

	nnz = 0;
	for (row = 0; row < n; row++) {
		ri[row] = nnz;
		for (col = 0; col < n; col++) {
			if (A.val[A.i2i(row, col)] > DC_EPS || A.val[A.i2i(row, col)] < (-1.0f)*DC_EPS) {
				ci[nnz] = col;
				val[nnz] = A.val[A.i2i(row, col)];
				nnz++;
			}
		}
	}
	ri[row] = nnz;
}

template <typename T> inline void MatrixS<T>::setDim(int __n, int nnz) {
	n = __n;
	numNonZeros = nnz;

	ri = new int[n + 1];
	val = new T[numNonZeros];
	ci = new int[numNonZeros];

	for (int i = 0; i < n + 1; i++)
		ri[i] = 0;
	for (int i = 0; i < nnz; i++) {
		val[i] = (T)0;
		ci[i] = 0;
	}
}

template <typename T> inline void MatrixS<T>::setSparsity(const MatrixS<T>& A) {
	assert((n == A.n) && (numNonZeros == A.numNonZeros));

	for (int i = 0; i < n + 1; i++)
		ri[i] = A.ri[i];
	for (int i = 0; i < numNonZeros; i++) {
		ci[i] = A.ci[i];
	}
}

template <typename T> inline T&	MatrixS<T>::get(int i, int j) const {
	// check if row-i, column-j is non-zero-element 
	bool bNNZ = false;
	int k;
	for (k = ri[i]; k < ri[i + 1]; k++) {
		if (ci[k] == j) {
			bNNZ = true;
			break;
		}
	}

	return val[k];
}

template <typename T> template <typename U>
inline void MatrixS<T>::setSparsity(const MatrixS<U>& A) {
	assert((n == A.n) && (numNonZeros == A.numNonZeros));

	for (int i = 0; i < n + 1; i++)
		ri[i] = A.ri[i];
	for (int i = 0; i < numNonZeros; i++) {
		ci[i] = A.ci[i];
	}
}

template <typename T> inline void MatrixS<T>::zero() {
	for (int i = 0; i < numNonZeros; i++) {
		val[i] = (T)0;
	}
}

template <typename T> inline void MatrixS<T>::setValue(const MatrixS<T> A) {
	assert((n == A.n) && (numNonZeros == A.numNonZeros));

	for (int i = 0; i < n + 1; i++)
		ri[i] = A.ri[i];
	for (int i = 0; i < numNonZeros; i++) {
		val[i] = A.val[i];
		ci[i] = A.ci[i];
	}
}

template <typename T> inline void MatrixS<T>::setValue(int i, int j, T newVal) {
	//// i행 j열이 non zero element 인지 판단
	bool bNNZ = false;
	int k;
	for (k = ri[i]; k < ri[i + 1]; k++) {
		if (ci[k] == j) {
			bNNZ = true;
			break;
		}
	}
	assert(bNNZ);

	val[k] = newVal;
}

template <typename T> inline void MatrixS<T>::setValue(int i, bool isRow, T newVal) {
	if (isRow) { // i번째 행(Row)의 값을 모두 newVal로 변경
		for (int k = ri[i]; k < ri[i + 1]; k++) {
			val[k] = newVal;
		}
	}
	else {     // i번째 열(Col)의 값을 모두 newVal로 변경
		for (int k = 0; k < numNonZeros; k++) {
			if (ci[k] == i)
				val[k] = newVal;
		}
	}
}

template <typename T> inline void MatrixS<T>::setBlockMatrix3by3(int i, int j, const Matrix3<T> mat) {
	// i번째 j번째 Block에 3by3 Matrix값으로 변경
	setValueInc(3 * i + 0, 3 * j + 0, mat(0, 0));
	setValueInc(3 * i + 1, 3 * j + 0, mat(1, 0));
	setValueInc(3 * i + 2, 3 * j + 0, mat(2, 0));
	setValueInc(3 * i + 0, 3 * j + 1, mat(0, 1));
	setValueInc(3 * i + 1, 3 * j + 1, mat(1, 1));
	setValueInc(3 * i + 2, 3 * j + 1, mat(2, 1));
	setValueInc(3 * i + 0, 3 * j + 2, mat(0, 2));
	setValueInc(3 * i + 1, 3 * j + 2, mat(1, 2));
	setValueInc(3 * i + 2, 3 * j + 2, mat(2, 2));
}

template <typename T> inline void MatrixS<T>::setValueInc(int i, int j, T IncVal) {
	// i행 j열이 non zero element 인지 판단
	bool bNNZ = false;
	int k;
	for (k = ri[i]; k < ri[i + 1]; k++) {
		if (ci[k] == j) {
			bNNZ = true;
			break;
		}
	}
	assert(bNNZ);

	val[k] += IncVal;
}

template <typename T> inline void MatrixS<T>::setValueMul(int i, int j, T MulVal) {
	// i행 j열이 non zero element 인지 판단
	bool bNNZ = false;
	int k;
	for (k = ri[i]; k < ri[i + 1]; k++) {
		if (ci[k] == j) {
			bNNZ = true;
			break;
		}
	}
	assert(bNNZ);

	val[k] *= MulVal;
}

template <typename T> inline void mul(MatrixS<T> &A, VectorN<T> &x, VectorN<T> &b) {
	int ri, idxCol;
	int nnz = A.nnz();

	ri = idxCol = 0;
	b = 0.f;

	for (int i = 0; i < A.n; i++) {
		b[i] = (T)0;
		for (int j = A.ri[i]; j < A.ri[i + 1]; j++) {
			b[i] += A.val[j] * x[A.ci[j]];
		}
	}
}

template <typename T> inline void mul(const MatrixS< Matrix3<T> > &A, const std::vector< Vector3<T> > &x, std::vector< Vector3<T> > &b) {
	assert(A.n == (int)x.size() && x.size() == b.size());
	const int& N = A.dim();

#pragma omp parallel for 
	for (int i = 0; i < N; ++i) {
		const int& end = A.ri[i + 1];
		Vector3d bi;
		for (int j = A.ri[i]; j < end; ++j) {
			mul_add(A.val[j], x[A.ci[j]], bi);
		}
		b[i] = bi;
	}
}

//-------------------------------------------------------------------
// IO
template <typename T> inline void MatrixS<T>::printMatrix() const {
	cout << "# of row : " << n << endl;
	cout << "# of nnz : " << numNonZeros << endl;
	cout << "[ ";
	for (int i = 0; i < n + 1; i++) {
		cout << ri[i] << "    ";
	}
	cout << "] " << endl;
	cout << "[ ";
	for (int i = 0; i < numNonZeros; i++) {
		cout << ci[i] << "    ";
	}
	cout << "] " << endl;
	cout << "[ ";
	for (int i = 0; i < numNonZeros; i++) {
		cout << val[i] << "    ";
	}
	cout << "] " << endl;
}

template <typename T> inline void MatrixS<T>::printMatrix(char *filename) const {
	FILE *fp = fopen(filename, "w");

	fprintf(fp, "[ ");
	for (int i = 0; i < n + 1; i++) {
		fprintf(fp, "%d  ", ri[i]);
	}
	fprintf(fp, " ] \n");

	fprintf(fp, "[ ");
	for (int i = 0; i < numNonZeros; i++) {
		fprintf(fp, "%d  ", ci[i]);
	}
	fprintf(fp, " ] \n");

	fprintf(fp, "[ ");
	for (int i = 0; i < numNonZeros; i++) {
		fprintf(fp, "%0.4lf  ", val[i]);
	}
	fprintf(fp, " ] \n");

	fclose(fp);
}

template <typename T> inline void mul_dif(const MatrixS< Matrix3<T> > &A, const std::vector< Vector3<T> > &x, const std::vector< Vector3<T> > &b, std::vector< Vector3<T> > &r) {
	assert(A.n == (int)x.size() && x.size() == b.size());
	const int& N = A.dim();

#pragma omp parallel for 
	for (int i = 0; i < N; ++i) {
		const int& end = A.ri[i + 1];
		Vector3d ti;
		for (int j = A.ri[i]; j < end; ++j) {
			mul_add(A.val[j], x[A.ci[j]], ti);
		}
		r[i] = b[i] - ti;
	}
}

#endif