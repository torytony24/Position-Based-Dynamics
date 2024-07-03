#ifndef __VECTOR_H__
#define __VECTOR_H__

#ifdef DEBUG
#define DB_CHECK( C ) { if ( ! (C) ) { abort(); } }
#else
#define DB_CHECK( C ) { }
#endif

#include <assert.h>
#include "vector2.h"
#include "vector3.h"
#include "vector4.h"
#include "vectorN.h"

typedef Vector2<int>				Vector2i;
typedef Vector2<float>				Vector2f;
typedef Vector2<double>				Vector2d;

typedef Vector3<int>				Vector3i;
typedef Vector3<unsigned int>		Vector3ui;
typedef Vector3<long int>			Vector3li;
typedef Vector3<unsigned long int>	Vector3uli;
typedef Vector3<float>				Vector3f;
typedef Vector3<double>				Vector3d;

typedef Vector4<float>				Vector4f;
typedef Vector4<double>				Vector4d;

typedef VectorN<int>				VectorNi;
typedef VectorN<unsigned int>		VectorNui;
typedef VectorN<float>				VectorNf;
typedef VectorN<double>				VectorNd;

#endif // __VECTOR_H__