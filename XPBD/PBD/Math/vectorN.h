#ifndef __VECTORN_H__
#define __VECTORN_H__


#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>

//============================================================================
//	3.1. VectorN definition
//	Class to store a N dimensional vector. 
template <typename T> class VectorN {
public:
	//////////////////////////////////////////////////////////////////////////
	// Constructors
	inline VectorN();						//!< Default constructor
	inline VectorN(int N);				//!< constructor.
	inline VectorN(int N, T d);			//!< Constructor.
											//!< @param d - default value.
	inline VectorN(VectorN &d);			//!< Constructor.
											//!< @param d - default constant vector.
	inline VectorN(int N, const T *d);	//!< Constructor.
											//!< @param N - size of vector.
											//!< @param d - pointer of data Array to be assinged as default constant value.
											// d should point to a T[N] that will be copied

	~VectorN() { if (n > 0) { delete[] val; n = 0; val = NULL; } }	//!< Destructor.

	//////////////////////////////////////////////////////////////////////////
	// Index operators
	inline T &get(unsigned int i);				//!< Getter by reference.
												//!< @param i - index for i_th component.
	inline T  get(unsigned int i) const;		//!< Getter by value.
												//!< @param i - index for i_th component.

	inline T &operator[](unsigned int i);		//!< Getter by reference.
												//!< @param i - index for i_th component.
	inline T  operator[](unsigned int i) const;	//!< Getter by value.
												//!< @param i - index for i_th component.

	inline T &operator()(unsigned int i);		//!< Getter by reference.
												//!< @param i - index for i_th component.
	inline T  operator()(unsigned int i) const;	//!< Getter by value.
												//!< @param i - index for i_th component.

	//////////////////////////////////////////////////////////////////////////
	// Size
	inline int dim() { return n; };	//!< Returns Dimension of the VectorN.

	//////////////////////////////////////////////////////////////////////////
	// Set dimension
	inline void setDim(const int i);

	////////////////////////////////////////////////////////////////////////
	// Simple Setter
	inline void zero();

	//////////////////////////////////////////////////////////////////////////
	// Assignment and set
	inline VectorN &set(T d);					//!< Global Setter.
												//!< @param d - Constant value to be assigned.
	inline VectorN &set(const VectorN &d);	//!< Component-wise Setter sourced from VectorN.
												//!< @param d - VectorN to be assinged.
	inline VectorN &set(const T *d);			//!< Component-wise Setter sourced form data array.
												//!< @param d - pointer of data Array to be assinged.
			// d should point to a T[N] that will be copied

	inline VectorN &operator=(T d);					//!< Assignment operation.
													//!< @param d - value to be assigned.
	inline VectorN &operator=(const VectorN &d);	//!< Assignment operation.
													//!< @param d - constant VectorN to be assigned.
	inline VectorN &operator=(const T    *d);		//!< Assignment operation.
													//!< @param d - constant pointer of data array to be assigned.

	//////////////////////////////////////////////////////////////////////////
	// Comparison operators
	inline int isEqual(const VectorN &d) const;		//!< Equality comparision function.
														//!< @param d - constant VectorN to be compared.
	inline int isNotEqual(const VectorN &d) const;	//!< Inequality comparision function.
														//!< @param d - constant VectorN to be compared.

	inline int operator==(const VectorN &d) const;	//!< Equality comparision operator.
														//!< @param d - constant VectorN to be compared.
	inline int operator!=(const VectorN &d) const;	//!< Inequality comparision operator.
														//!< @param d - constant VectorN to be compared.

	inline int operator==(T d) const;					//!< Equality comparision operator.
														//!< @param d - value to be compared.
	inline int operator!=(T d) const;					//!< Inequality comparision operator.
														//!< @param d - value to be compared.

	//////////////////////////////////////////////////////////////////////////
	// In place arithmetic
	inline VectorN &add(T d);			//!< Global addition function.
										//!< @param d - value to be added.
	inline VectorN &subtract(T d);		//!< Global subtraction function.
										//!< @param d - value to be subtractd.
	inline VectorN &multiply(T d);		//!< Global multiplication function.
										//!< @param d - value to be multiplicated.
	inline VectorN &divide(T d);		//!< Global division function.
										//!< @param d - divider value.

	inline VectorN &operator+=(T d);	//!< Global addition operator.
										//!< @param d - value to be added.
	inline VectorN &operator-=(T d);	//!< Global subtraction operator.
										//!< @param d - value to be subtractd.
	inline VectorN &operator*=(T d);	//!< Global multiplication operator.
										//!< @param d - value to be multiplicated.
	inline VectorN &operator/=(T d);	//!< Global division operator.
										//!< @param d - divider value.

	//////////////////////////////////////////////////////////////////////////
	// In place arithmetic; componentwise operations
	inline VectorN &operator+=(const VectorN &d);		//!< Component-wise addition operator.
														//!< @param d - constant VectorN to be added.
	inline VectorN &operator-=(const VectorN &d);		//!< Component-wise subtraction operator.
														//!< @param d - constant VectorN to be subtracted.
	inline VectorN &operator*=(const VectorN &d);		//!< Component-wise multiplication operator.
														//!< @param d - constant VectorN to be multiplied.
	inline VectorN &operator/=(const VectorN &d);		//!< Component-wise division operator.
														//!< @param d - constant divider VectorN.


	//------------------------------------------------------------------------
	// IO
	void printVector() const;
	void printVector(char *filename) const;

	T *val;

private:
	int n;
};

//-------------------------------------------------------------------
// 3.2. Inline implementation of VectorN
// Constructors
template <typename T> inline VectorN<T>::VectorN() {
	n = 0;
	val = NULL;
}

template <typename T> inline VectorN<T>::VectorN(int N) {
	n = N;
	val = new T[N];

	for (int i = 0; i < n; i++)
		val[i] = 0;
}

template <typename T> inline VectorN<T>::VectorN(int N, T d) {
	n = N;
	val = new T[N];

	for (int i = 0; i < n; i++)
		val[i] = d;
}

template <typename T> inline VectorN<T>::VectorN(VectorN &d) {
	n = d.dim();
	val = new T[n];

	for (int i = 0; i < n; i++)
		val[i] = d[i];
}

template <typename T> inline VectorN<T>::VectorN(int N, const T *d) {
	//DB_CHECK( N <= d.dim() );
	n = N;
	val = new T[n];

	for (int i = 0; i < n; i++)
		val[i] = d[i];
}

// Index operators
template <typename T> inline T &VectorN<T>::get(unsigned int i) {
	//DB_CHECK( i >= n );
	return val[i];
}

template <typename T> inline T  VectorN<T>::get(unsigned int i) const {
	//DB_CHECK( (int)i >= n );
	return val[i];
}

template <typename T> inline T &VectorN<T>::operator[](unsigned int i) {
	//DB_CHECK( (int)i >= (int)n );
	return val[i];
}

template <typename T> inline T  VectorN<T>::operator[](unsigned int i) const {
	//DB_CHECK( i >= n );

	return val[i];
}

template <typename T> inline T &VectorN<T>::operator()(unsigned int i) {
	//DB_CHECK( i >= n );
	return val[i];
}

template <typename T> inline T  VectorN<T>::operator()(unsigned int i) const {
	//DB_CHECK( i >= n );
	return val[i];
}

// Simple Setter
template <typename T> inline void VectorN<T>::zero() {
	for (int i = 0; i < n; i++)
		val[i] = (T)0.0;
}

// set Dimension
template <typename T> inline void VectorN<T>::setDim(const int N) {
	if (n > 0)
		delete[] val;

	n = N;
	val = new T[N];

	for (int i = 0; i < n; i++)
		val[i] = 0;
}


// Assignment and set
template <typename T> inline VectorN<T> &VectorN<T>::set(T d) {
	for (int i = 0; i < n; i++)
		val[i] = d;

	return (*this);
}

template <typename T> inline VectorN<T> &VectorN<T>::set(const VectorN &d) {
	//DB_CHECK( n != d.dim() );
	for (int i = 0; i < n; i++)
		val[i] = d[i];

	return (*this);
}

template <typename T> inline VectorN<T> &VectorN<T>::operator=(T d) {
	return set(d);
}

template <typename T> inline VectorN<T> &VectorN<T>::operator=(const VectorN &d) {
	return set(d);
}

	
// In place arithmetic
template <typename T> inline VectorN<T> &VectorN<T>::add(T d) {
	for (int i = 0; i < n; i++)
		val[i] += d;
	

	return (*this);
}

template <typename T> inline VectorN<T> &VectorN<T>::subtract(T d) {
	for (int i = 0; i < n; i++)
		val[i] -= d;

	return (*this);
}

template <typename T> inline VectorN<T> &VectorN<T>::multiply(T d) {
	for (int i = 0; i < n; i++)
		val[i] *= d;

	return (*this);
}

template <typename T> inline VectorN<T> &VectorN<T>::divide(T d) {
	for (int i = 0; i < n; i++)
		val[i] /= d;

	return (*this);
}

template <typename T> inline VectorN<T> &VectorN<T>::operator+=(T d) {
	for (int i = 0; i < n; i++)
		val[i] += d;

	return (*this);
}

template <typename T> inline VectorN<T> &VectorN<T>::operator-=(T d) {
	for (int i = 0; i < n; i++)
		val[i] -= d;

	return (*this);
}

template <typename T> inline VectorN<T> &VectorN<T>::operator*=(T d) {
	for (int i = 0; i < n; i++)
		val[i] *= d;

	return (*this);
}

template <typename T> inline VectorN<T> &VectorN<T>::operator/=(T d) {
	for (int i = 0; i < n; i++)
		val[i] /= d;

	return (*this);
}

// In place arithmetic; componentwise operations
template <typename T> inline VectorN<T> &VectorN<T>::operator+=(const VectorN &d) {
	//DB_CHECK( n != d.dim() );

	for (int i = 0; i < n; i++)
		val[i] += d[i];

	return (*this);
}

template <typename T> inline VectorN<T> &VectorN<T>::operator-=(const VectorN &d) {
	//DB_CHECK( n != d.dim() );
	for (int i = 0; i < n; i++)
		val[i] -= d[i];

	return (*this);
}

template <typename T> inline VectorN<T> &VectorN<T>::operator*=(const VectorN &d) {
	//DB_CHECK( n != d.dim() );
	for (int i = 0; i < n; i++)
		val[i] *= d[i];

	return (*this);
}

template <typename T> inline VectorN<T> &VectorN<T>::operator/=(const VectorN &d) {
	//DB_CHECK( n != d.dim() );
	for (int i = 0; i < n; i++)
		val[i] /= d[i];

	return (*this);
}



//-------------------------------------------------------------------
// 3.3. Inline implementation of functions related to VectorN
//! Addition of two vectors.
//! @param a - VectorN to be added.
//! @param b - VectorN to be added.
//! @param c - VectorN to store a+b.
template <typename T> inline void add(VectorN<T> &a, VectorN<T> &b, VectorN<T> &c) {
	//DB_CHECK( a.dim() != b.dim() || b.dim() != c.dim() );
	int N = a.dim();
	for (int i = 0; i < N; i++)
		c[i] = a[i] + b[i];

}

template <typename T> inline void dif(VectorN<T> &a, VectorN<T> &b, VectorN<T> &c) {
	//DB_CHECK( a.dim() != b.dim() || b.dim() != c.dim() );
	int N = a.dim();
	for (int i = 0; i < N; i++)
		c[i] = a[i] - b[i];

}

template <typename T> inline void mul(T &a, VectorN<T> &b, VectorN<T> &c) {
	//DB_CHECK( b.dim() != c.dim() );
	int N = b.dim();
	for (int i = 0; i < N; i++)
		c[i] = a * b[i];

}

template <typename T> inline T dot(VectorN<T> &a, VectorN<T> &b) {
	//DB_CHECK( a.dim() != b.dim() );
	int N = a.dim();
	T tmp = 0;
	for (int i = 0; i < N; i++)
		tmp += a[i] * b[i];

	return tmp;
}

template <typename T> std::ostream &operator<<(std::ostream &strm, VectorN<T> &v) {
	strm << "[";
	int N = v.dim();
	for (int i = 0; i < N - 1; i++) strm << v[i] << ", ";
	strm << v[N - 1] << "]";
	return strm;
}



//-------------------------------------------------------------------
// IO
template <typename T> inline void VectorN<T>::printVector() const {
	for (int i = 0; i < n; i++)
		cout << val[i] << "  ";

	cout << "\n";
}

template <typename T> inline void VectorN<T>::printVector(char *filename) const {
	FILE *fp = fopen(filename, "w");

	fprintf(fp, "[ ");
	for (int i = 0; i < n; i++) {
		fprintf(fp, "%lf ", val[i]);
		if (i != n - 1)
			fprintf(fp, " ;");
		else
			fprintf(fp, " ]\n");

	}
	fclose(fp);
}

typedef VectorN<float>				VectorNf;
typedef VectorN<double>				VectorNd;

#endif // __VECTORN_H__
