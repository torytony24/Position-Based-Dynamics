#pragma once

//============================================================================
//1. Constants

//! Constant pi 
const double PI = 3.14159265358979323846;		// ¥ğ

//! Constant 2*pi
const double TWO_PI = 6.28318530717959;

//! Constant 1/pi 
const double INV_PI = 0.31830988618379;

//! Constant 1/(2*pi) 
const double INV_TWO_PI = 0.1591549430919;

//! Constant pi/180 
const double D2R = 0.01745329251994;

//! Constant 180/pi
const double R2D = 57.2957795130931;

//! Constant 1/sqrt(2)
const double ONE_OVER_ROOT2 = 0.70710678118654;

//! Constant epsilon
const double EPS = 1.0e-6;

//! Constant epsilon
const double EPS4 = 1.0e-4;

//! Constant epsilon
const double EPS3 = 1.0e-3;

//! Constant epsilon
const double EPS2 = 1.0e-2;

//! Constant epsilon
const double EPS1 = 1.0e-1;

//! Constant ROOT2
const float SQRT2 = 1.4142135623730951f;
const float ROOT2f = (float)1.4142135623730951;
const double ROOT2d = 1.4142135623730951;

const float RadToDeg = 180.0f / PI;
const float DegToRad = PI / 180.0f;

#define DTOR		0.0174532925
#define RTOD		57.324840
#define DC_EPS			1.0e-4
#define ERROR_BOUND		0
#define     GLH_ZERO            float(0.0)
#define     GLH_EPSILON         float(10e-6)
#define		GLH_EPSILON_2		float(10e-12)

//============================================================================
// 2. Arithmatics

#ifndef ABS
#define ABS(x)		( ((x)>=0.0) ? (x) :-(x) )
#endif
#define ABSMIN(a,b)	( ((a*a)<(b*b)) ? (a) : (b) )
#define ABSMAX(a,b)	( ((a*a)>(b*b)) ? (a) : (b) )
#define SQ(x)		( (x)*(x) )
#define CB(x)		( (x)*(x)*(x) )
#define SQRT(x)		( sqrt(ABS(x)) )
#define MAG(x,y,z)	( sqrt( (x)*(x) + (y)*(y) + (z)*(z) ) )
#define MAG2(x)		( sqrt( (x[0])*(x[0]) + (x[1])*(x[1]) ) )
#define MAG3(x)		( sqrt( (x[0])*(x[0]) + (x[1])*(x[1]) + (x[2])*(x[2]) ) )
#define DIST2(x1,y1,x2,y2) ( sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) ) )
#define DOT2(x,y)	( x[0]*y[0] + x[1]*y[1] )
#define DOT3(x,y)	( x[0]*y[0] + x[1]*y[1] + x[2]*y[2] )
#define AVG2(x1,x2) ( 0.5*( (x1)+(x2) ) )
#define AVG4(x1,x2,x3,x4) ( 0.25*( (x1)+(x2)+(x3)+(x4) ) )
#define AVG8(x1,x2,x3,x4,x5,x6,x7,x8) ( 0.125*( (x1)+(x2)+(x3)+(x4)+(x5)+(x6)+(x7)+(x8) ) )
#define RANDOM(x)	( ((float)(rand()%(x)))/((float)x) )
#define RANDOM_SIGN	( (rand()%2 == 1) ? 1 : -1 )
#define ACOS(x)		( ((x)>1.0) ? (0) : ( ((x)<-1.0) ? (PI) : (acos(x)) ) )
#define ASIN(x)  ( ((x)>1.0) ? (PI/2.0) : ( ((x)<-1.0) ? (-PI/2.0) : (asin(x)) ) )
#define SWAP(a,b)	a^=b^a^=b

template <typename T>
inline void Set(T *&v, T val, int size) {
	for (int i = 0; i < size; ++i)
		v[i] = val;
}

template <typename T> inline void Swap(T &x, T &y) { T tmp = x; x = y; y = tmp; }
template <typename T> inline void Release(T *&x) { if (x != NULL) delete(x); x = NULL; }

template <typename T> inline T Rad(const T x) { return (float)DegToRad * x; }
template <typename T> inline T Deg(const T x) { return (float)RadToDeg * x; }

template <typename T> inline T MIN(const T a, const T b) { return (a < b) ? a : b; }
template <typename T> inline T MAX(const T a, const T b) { return (a > b) ? a : b; }
template <typename T> inline T MIN3(const T a, const T b, const T c) { return MIN(a, MIN(b, c)); }
template <typename T> inline T MAX3(const T a, const T b, const T c) { return MAX(a, MAX(b, c)); }

template <typename T> inline T SQUARE(const T a) { return a * a; }
template <typename T> inline T POW2(const T a) { return a * a; }
template <typename T> inline T CUBE(const T a) { return a * a*a; }
template <typename T> inline T POW3(const T a) { return a * a*a; }

template <typename T>
inline bool IsAlmostSame(const T x, const T y) {
	T a = x - y;
	if (a < 0)
		a = -a;
	if (a <= EPS) return true;
	else return false;
}

template <typename T>
inline bool IsAlmostZero(const T x) {
	T a = x;
	if (a < 0)
		a = -a;
	if (a <= EPS) return true;
	else return false;
}

template <typename T>
inline int SIGN(const T x) {
	if (x > 0)	  return  1;
	else		  return -1;
}

template<class T>
inline const T SIGN(const T &a, const T &b) {
	return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
}

inline const float SIGN(const float &a, const double &b) {
	return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
}

inline const float SIGN(const double &a, const float &b) {
	return (float)(b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a));
}

template <typename T>
inline int Sign(const T x) {
	if (IsAlmostZero(x)) return  0;
	else if (x > 0)	  return  1;
	else				  return -1;
}

template <typename T>
inline T clamp(T val, T min, T max) {
	if (val < min)		return min;
	else if (val > max)	return max;
	return val;
}

inline double pow2(int p) {
	switch (p) {
	case -20: return 9.53674e-07;
	case -19: return 1.90735e-06;
	case -18: return 3.8147e-06;
	case -17: return 7.62939e-06;
	case -16: return 1.52588e-05;
	case -15: return 3.05176e-05;
	case -14: return 6.10352e-05;
	case -13: return 0.0001220703125;
	case -12: return 0.000244140625;
	case -11: return 0.00048828125;
	case -10: return 0.0009765625;
	case -9: return 0.001953125;
	case -8: return 0.00390625;
	case -7: return 0.0078125;
	case -6: return 0.015625;
	case -5: return 0.03125;
	case -4: return 0.0625;
	case -3: return 0.125;
	case -2: return 0.25;
	case -1: return 0.5;
	case 0: return 1;
	case 1: return 2;
	case 2: return 4;
	case 3: return 8;
	case 4: return 16;
	case 5: return 32;
	case 6: return 64;
	case 7: return 128;
	case 8: return 256;
	case 9: return 512;
	case 10: return 1024;
	case 11: return 2048;
	case 12: return 4096;
	case 13: return 8192;
	case 14: return 16384;
	case 15: return 32768;
	case 16: return 65536;
	case 17: return 131072;
	case 18: return 262144;
	case 19: return 524288;
	case 20: return 1048576;
	default:
		double ret = 1;
		if (ABS(p) == p)
			for (int i = 0; i < ABS(p); i++) ret *= 2.0;
		else
			for (int i = 0; i < ABS(p); i++) ret *= 0.5;
		return ret;
	}
}

inline int Obin(int p) {
	int ret = 2;
	for (int i = 0; i < p * 2; i *= 2) {
		if (p < ret) return ret;
		else ret *= 2;
	}
	return ret;
}

inline int Odec(int p) {
	int ret = 10;
	for (int i = 0; i < p * 10; i *= 10) {
		if (p < ret) return ret;
		else ret *= 10;
	}
	return ret;
}

//============================================================================
// 3. Debug

#define PRINT_INT(src) printf(#src ": %d\n", src)
#define PRINT_DOUBLE(src) printf(#src ": %f\n", src)
#define PRINT_ERROR(src) printf("ERROR: %s\n", src)
#define FLAG printf("FLAG\n")
