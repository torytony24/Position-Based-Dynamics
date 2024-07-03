#pragma once

#include "vector2.h"
//============================================================================
//	Box2 class definition

template <typename T> class Box2 {
public:
	// constructor
	Box2();
	Box2(Vector2<T> p1, Vector2<T> p2);
	Box2(Vector2<T> p1, Vector2<T> p2, Vector2<T> p3);

	// setter
	Box2&	set(const Box2& ab);
	Box2&	set(Vector2<T> p1);
	Box2&	set(Vector2<T> p1, Vector2<T> p2);
	Box2&	set(Vector2<T> p1, Vector2<T> p2, Vector2<T> p3);
	void	setMin(Vector2<T> p);
	void	setMax(Vector2<T> p);

	// getter
	inline T lengthX() const { return std::abs(max[0] - min[0]); }
	inline T lengthY() const { return std::abs(max[1] - min[1]); }
	inline T minX() const { return min[0]; }
	inline T minY() const { return min[1]; }
	inline T maxX() const { return max[0]; }
	inline T maxY() const { return max[1]; }
	inline T centerX() const { return 0.5 * (min[0] + max[0]); }
	inline T centerY() const { return 0.5 * (min[1] + max[1]); }
	inline Vector2<T> center() const { return Vector2<T>((T)0.5 * (min[0] + max[0]), (T)0.5 * (min[1] + max[1])); }

	// expand
	inline void expandX(T x) { min[0] -= x; max[0] += x; }
	inline void expandY(T y) { min[1] -= y; max[1] += y; }
	inline void expand(T t) { for (int i = 0; i < 2; i++) { min[i] -= t; max[i] += t; } }

	// scale
	inline void scaleX(T x) { min[0] *= x; max[0] *= x; }
	inline void scaleY(T y) { min[1] *= y; max[1] *= y; }
	inline void scale(T t) { for (int i = 0; i < 2; i++) { min[i] *= t; max[i] *= t; } }

	inline Box2& operator=(const Box2& ab);

	inline Box2& add(const Box2& ab);
	inline Box2& add(const Vector2<T>& v);

	inline void clear() { min.set((T)0); max.set((T)0); }

	inline bool contain(const Vector2<T> p) const;
	inline bool contain(const Vector2<T> p, const T tolr) const;

public:
	Vector2<T> min;
	Vector2<T> max;
};


template <typename T> inline Box2<T>::Box2() {
	clear();
}

template <typename T> inline Box2<T>::Box2(Vector2<T> p1, Vector2<T> p2) {
	set(p1, p2);
}

template <typename T> inline Box2<T>::Box2(Vector2<T> p1, Vector2<T> p2, Vector2<T> p3) {
	set(p1, p2, p3);
}

template <typename T> inline Box2<T>& Box2<T>::set(const Box2& ab) {
	min = ab.min;
	max = ab.max;
	return *this;
}

template <typename T> inline Box2<T>& Box2<T>::set(Vector2<T> p1) {
	min = p1;
	max = p1;
	return *this;
}

template <typename T> inline Box2<T>& Box2<T>::set(Vector2<T> p1, Vector2<T> p2) {
	for (int i = 0; i < 2; ++i) {
		min[i] = MIN(p1[i], p2[i]);
		max[i] = MAX(p1[i], p2[i]);
	}
	return *this;
}

template <typename T> inline Box2<T>& Box2<T>::set(Vector2<T> p1, Vector2<T> p2, Vector2<T> p3) {
	for (int i = 0; i < 2; ++i) {
		min[i] = MIN(p1[i], p2[i], p3[i]);
		max[i] = MAX(p1[i], p2[i], p3[i]);
	}
	return *this;
}

template <typename T> inline void Box2<T>::setMin(Vector2<T> p) {
	min = p;
}

template <typename T> inline void Box2<T>::setMax(Vector2<T> p) {
	max = p;
}

template <typename T> inline Box2<T> &Box2<T>::operator=( const Box2& ab ) {
	return set(ab);
}

template <typename T> inline Box2<T> &Box2<T>::add(const Box2& ab) {
	for (int i = 0; i < 2; ++i) {
		if (min[i] > ab.min[i])
			min[i] = ab.min[i];
		if (max[i] < ab.max[i])
			max[i] = ab.max[i];
	}
	return *this;
}

template <typename T> inline Box2<T> &Box2<T>::add(const Vector2<T>& v) {
	for (int i = 0; i < 2; ++i) {
		if (min[i] > v[i])
			min[i] = v[i];
		if (max[i] < v[i])
			max[i] = v[i];
	}

	return *this;
}

template <typename T> inline bool Box2<T>::contain(const Vector2<T> p) const {
	for (int i = 0; i < 2; ++i) {
		if (p[i] < min[i])
			return false;
		if (p[i] > max[i])
			return false;
	}

	return true;
}

template <typename T> inline bool Box2<T>::contain(const Vector2<T> p, const T tolr) const {
	for (int i = 0; i < 2; ++i) {
		if (p[i] < min[i] - tolr)
			return false;
		if (p[i] > max[i] + tolr)
			return false;
	}

	return true;
}

template <typename T> inline bool intersect(const Box2<T>& ab1, const Box2<T>& ab2) {
	if (ab1.min[0] > ab2.max[0]) return false;
	if (ab2.min[0] > ab1.max[0]) return false;

	if (ab1.min[1] > ab2.max[1]) return false;
	if (ab2.min[1] > ab1.max[1]) return false;

	return true;
}

template <typename T> inline bool sameSize(const Box2<T>& ab1, const Box2<T>& ab2, const T& tolr) {
	if (!almostEqual(ab1.lengthX(), ab2.lengthX(), tolr))
		return false;
	if (!almostEqual(ab1.lengthY(), ab2.lengthY(), tolr))
		return false;

	return true;
}

template <typename T> inline bool almostEqualSize(const Box2<T>& ab1, const Box2<T>& ab2) {
	if (!almostEqual(ab1.lengthX(), ab2.lengthX()))
		return false;
	if (!almostEqual(ab1.lengthY(), ab2.lengthY()))
		return false;

	return true;
}