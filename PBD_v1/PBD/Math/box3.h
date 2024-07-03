//============================================================================
//
//  Author: Inyong Jeon, jeon@graphics.snu.ac.kr
//  Copyright (c) Inyong Jeon. All rights reserved.
//	Last update: 2014. 02. 26
//
//============================================================================

#pragma once

#include "vector3.h"
//============================================================================
//	Box3 class definition

template <typename T> class Box3 {
public:
	// constructor
	Box3();
	Box3(Vector3<T> p1, Vector3<T> p2);
	Box3(Vector3<T> p1, Vector3<T> p2, Vector3<T> p3);

	// setter
	Box3&	set(const Box3& ab);
	Box3&	set(Vector3<T> p1);
	Box3&	set(Vector3<T> p1, Vector3<T> p2);
	Box3&	set(Vector3<T> p1, Vector3<T> p2, Vector3<T> p3);
	void	setMin(Vector3<T> p);
	void	setMax(Vector3<T> p);

	// getter
	inline T lengthX() const { return std::abs(max[0] - min[0]); }
	inline T lengthY() const { return std::abs(max[1] - min[1]); }
	inline T lengthZ() const { return std::abs(max[2] - min[2]); }
	inline T minX() const { return min[0]; }
	inline T minY() const { return min[1]; }
	inline T minZ() const { return min[2]; }
	inline T maxX() const { return max[0]; }
	inline T maxY() const { return max[1]; }
	inline T maxZ() const { return max[2]; }
	inline T centerX() const { return 0.5 * (min[0] + max[0]); }
	inline T centerY() const { return 0.5 * (min[1] + max[1]); }
	inline T centerZ() const { return 0.5 * (min[2] + max[2]); }

	inline Vector3<T> center() const { return Vector3<T>((T)0.5 * (min[0] + max[0]), (T)0.5 * (min[1] + max[1]), (T)0.5 * (min[2] + max[2])); }

	// expand
	inline void expandX(T x) { min[0] -= x; max[0] += x; }
	inline void expandY(T y) { min[1] -= y; max[1] += y; }
	inline void expandZ(T z) { min[2] -= z; max[2] += z; }
	inline void expand(T t) { for (int i = 0; i < 3; i++) { min[i] -= t; max[i] += t; } }

	inline Box3& operator=(const Box3& ab);

	inline Box3& add(const Box3& ab);
	inline Box3& add(const Vector3<T>& v);

	inline void	clear() { min.set((T)0); max.set((T)0); }

	inline bool contain(const Vector3<T> p) const;
	inline bool contain(const Vector3<T> p, const T tolr) const;

public:
	Vector3<T> min;
	Vector3<T> max;
};


template <typename T> inline Box3<T>::Box3() {
	clear();
}

template <typename T> inline Box3<T>::Box3(Vector3<T> p1, Vector3<T> p2) {
	set(p1, p2);
}

template <typename T> inline Box3<T>::Box3(Vector3<T> p1, Vector3<T> p2, Vector3<T> p3) {
	set(p1, p2, p3);
}

template <typename T> inline Box3<T>& Box3<T>::set(const Box3& ab) {
	min = ab.min;
	max = ab.max;
	return *this;
}

template <typename T> inline Box3<T>& Box3<T>::set(Vector3<T> p1) {
	min = p1;
	max = p1;
	return *this;
}

template <typename T> inline Box3<T>& Box3<T>::set(Vector3<T> p1, Vector3<T> p2) {
	for (int i = 0; i < 3; ++i) {
		min[i] = MIN(p1[i], p2[i]);
		max[i] = MAX(p1[i], p2[i]);
	}
	return *this;
}

template <typename T> inline Box3<T>& Box3<T>::set(Vector3<T> p1, Vector3<T> p2, Vector3<T> p3) {
	for (int i = 0; i < 3; ++i) {
		min[i] = MIN(p1[i], p2[i], p3[i]);
		max[i] = MAX(p1[i], p2[i], p3[i]);
	}
	return *this;
}

template <typename T> inline void Box3<T>::setMin(Vector3<T> p) {
	min = p;
}

template <typename T> inline void Box3<T>::setMax(Vector3<T> p) {
	max = p;
}

template <typename T> inline Box3<T> &Box3<T>::operator=( const Box3& ab ) {
	return set(ab);
}

template <typename T> inline Box3<T> &Box3<T>::add(const Box3& ab) {
	for (int i = 0; i < 3; ++i) {
		if (min[i] > ab.min[i])
			min[i] = ab.min[i];
		if (max[i] < ab.max[i])
			max[i] = ab.max[i];
	}
	return *this;
}

template <typename T> inline Box3<T> &Box3<T>::add(const Vector3<T>& v) {
	for (int i = 0; i < 3; ++i) {
		if (min[i] > v[i])
			min[i] = v[i];
		if (max[i] < v[i])
			max[i] = v[i];
	}
	return *this;
}

template <typename T> inline bool Box3<T>::contain(const Vector3<T> p) const {
	for (int i = 0; i < 3; ++i) {
		if (p[i] < min[i])
			return false;
		if (p[i] > max[i])
			return false;
	}

	return true;
}

template <typename T> inline bool Box3<T>::contain(const Vector3<T> p, const T tolr) const {
	for (int i = 0; i < 3; ++i) {
		if (p[i] < min[i] - tolr)
			return false;
		if (p[i] > max[i] + tolr)
			return false;
	}

	return true;
}

template <typename T> inline bool intersect(const Box3<T>& ab1, const Box3<T>& ab2) {
	if (ab1.min[0] > ab2.max[0]) return false;
	if (ab2.min[0] > ab1.max[0]) return false;

	if (ab1.min[1] > ab2.max[1]) return false;
	if (ab2.min[1] > ab1.max[1]) return false;

	if (ab1.min[2] > ab2.max[2]) return false;
	if (ab2.min[2] > ab1.max[2]) return false;

	return true;
}

template <typename T> inline bool sameSize(const Box3<T>& ab1, const Box3<T>& ab2, const T& tolr) {
	if (!almostEqual(ab1.lengthX(), ab2.lengthX(), tolr))
		return false;
	if (!almostEqual(ab1.lengthY(), ab2.lengthY(), tolr))
		return false;
	if (!almostEqual(ab1.lengthZ(), ab2.lengthZ(), tolr))
		return false;

	return true;
}

template <typename T> inline bool almostEqualSize(const Box3<T>& ab1, const Box3<T>& ab2) {
	if (!almostEqual(ab1.lengthX(), ab2.lengthX()))
		return false;
	if (!almostEqual(ab1.lengthY(), ab2.lengthY()))
		return false;
	if (!almostEqual(ab1.lengthZ(), ab2.lengthZ()))
		return false;

	return true;
}