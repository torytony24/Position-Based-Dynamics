#ifndef __QUATERNION_H__
#define __QUATERNION_H__

#include "vector.h"
#include "matrix.h"
#include <math.h>
#include <iostream>


template <typename T> class Quaternion {

public:
	///////////////////////////////////////////////////////
	// Constructors
	Quaternion() { q[0] = q[1] = q[2] = 0.0;  q[3] = 1.0; }									//! Default constructor. Builds an identity rotation.
	Quaternion(const Vector3<T>& axis, T angle) { setAxisAngle(axis, angle); }	//! Constructor from rotation axis (non null) and angle (in radians).
	Quaternion(T q0, T q1, T q2, T q3) { q[0] = q0; q[1] = q1; q[2] = q2;  q[3] = q3; }
	Quaternion(const Vector3<T>& from, const Vector3<T>& to);
	Quaternion(const Quaternion& Q) { for (int i = 0; i < 4; ++i) q[i] = Q.q[i]; }

	////////////////////////////////////////////////////////////////////////////
	// Index operators
	inline T& operator[](int i) { return q[i]; }
	inline const T& operator[](int i) const { return q[i]; }


	///////////////////////////////////////////////////////
	// Setter 
	inline Quaternion& operator=(const Quaternion& Q) {
		memcpy(q, Q.q, 4 * sizeof(T));
		return (*this);
	}

	inline void setAxisAngle(const Vector3<T>& axis, T angle) {
		const T norm = l2Norm(axis);
		if (norm < 1e-8) {
			// Null rotation
			q[0] = 0.0; q[1] = 0.0; q[2] = 0.0;
			q[3] = 1.0;
		}
		else {
			const T sin_half_angle = (T)sin(0.5f*angle);
			q[0] = sin_half_angle * axis[0] / norm;
			q[1] = sin_half_angle * axis[1] / norm;
			q[2] = sin_half_angle * axis[2] / norm;
			q[3] = (T)cos(angle / 2.0);
		}
	}
	inline void setValue(T q0, T q1, T q2, T q3) { q[0] = q0; q[1] = q1; q[2] = q2;  q[3] = q3; }
	inline void setFromRotationMatrix(const T m[3][3]);
	inline void setFromRotatedBasis(const Vector3<T>& X, const Vector3<T>& Y, const Vector3<T>& Z);
	inline void setFromSO3(const Matrix3<T>& d); // Special Orthogonal Group
	inline void setFromSE3(const Matrix4<T>& d); // Special Euclidean  Group 
	inline void identity();

	///////////////////////////////////////////////////////
	// Getter 
	inline Vector3<T> axis() const;
	inline T angle() const;
	inline void getAxisAngle(Vector3<T>& axis, T& angle) const;

	//////////////////////////////////////////////////////////////////////////
	// In place arithmetic	
	inline friend Quaternion operator*(const Quaternion& a, const Quaternion& b) {
		return Quaternion(a.q[3] * b.q[0] + b.q[3] * a.q[0] + a.q[1] * b.q[2] - a.q[2] * b.q[1],
			a.q[3] * b.q[1] + b.q[3] * a.q[1] + a.q[2] * b.q[0] - a.q[0] * b.q[2],
			a.q[3] * b.q[2] + b.q[3] * a.q[2] + a.q[0] * b.q[1] - a.q[1] * b.q[0],
			a.q[3] * b.q[3] - b.q[0] * a.q[0] - a.q[1] * b.q[1] - a.q[2] * b.q[2]);
	}		// Returns the composition of the \p a and \p b rotations.
			// The order is important. When applied to a Vec \c v (see operator*(const Quaternion&, const Vector3<T>&)
			// and rotate()) the resulting Quaternion acts as if \p b was applied first and then \p a was
			// applied. This is obvious since the image \c v' of \p v by the composited rotation satisfies: \code
			// v'= (a*b) * v = a * (b*v) \endcode
			// Note that a*b usually differs from b*a.
			// For efficiency reasons, the resulting Quaternion is not normalized. Use normalize() in
			// case of numerical drift with small rotation composition.

	inline Quaternion& operator*=(const Quaternion &q) {
		*this = (*this)*q;
		return *this;
	}		// Quaternion rotation is composed with \p q.
			// See operator*(), since this is equivalent to \c this = \c this * \p q.
			// \note For efficiency reasons, the resulting Quaternion is not normalized.
			// You may normalize() it after each application in case of numerical drift. */


	inline friend Vector3<T> operator*(const Quaternion& q, const Vector3<T>& v) {
		return q.rotate(v);
	}		// Returns the image of \p v by the rotation \p q.
			// Same as q.rotate(v). See rotate() and inverseRotate().


	inline Vector3<T> rotate(const Vector3<T>& v) const;
	inline Vector3<T> inverseRotate(const Vector3<T>& v) const;

	inline Quaternion inverse() const { return Quaternion(-q[0], -q[1], -q[2], q[3]); }
	// Returns the inverse Quaternion (inverse rotation).
	// Result has a negated axis() direction and the same angle(). A composition (see operator*()) of a
	// Quaternion and its inverse() results in an identity function.
	// Use invert() to actually modify the Quaternion.


	inline void invert() { q[0] = -q[0]; q[1] = -q[1]; q[2] = -q[2]; }
	// Inverses the Quaternion (same rotation angle(), but negated axis()).
	// See also inverse().


	inline void negate() { invert(); q[3] = -q[3]; }
	// Negates all the coefficients of the Quaternion.
	// This results in an other representation of the \e same rotation (opposite rotation angle, but with
	// a negated axis direction: the two cancel out). However, note that the results of axis() and
	// angle() are unchanged after a call to this method since angle() always returns a value in [0,pi].
	// This method is mainly useful for Quaternion interpolation, so that the spherical
	// interpolation takes the shortest path on the unit sphere. See slerp() for details.


	inline T normalize() {
		const T norm = sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
		for (int i = 0; i < 4; ++i) q[i] /= norm;
		return norm;
	}		// Normalizes the Quaternion coefficients.
			// This method should not need to be called since we only deal with unit Quaternions. This is however
			// useful to prevent numerical drifts, especially with small rotational increments. See also
			// normalized().

	inline Quaternion normalized() const {
		T Q[4];
		const T norm = sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
		for (int i = 0; i < 4; ++i) Q[i] = q[i] / norm;
		return Quaternion(Q[0], Q[1], Q[2], Q[3]);
	}

	inline void getMatrix(T m[4][4]) const {
		const T q00 = (T)(2.0 * q[0] * q[0]);
		const T q11 = (T)(2.0 * q[1] * q[1]);
		const T q22 = (T)(2.0 * q[2] * q[2]);

		const T q01 = (T)(2.0 * q[0] * q[1]);
		const T q02 = (T)(2.0 * q[0] * q[2]);
		const T q03 = (T)(2.0 * q[0] * q[3]);

		const T q12 = (T)(2.0 * q[1] * q[2]);
		const T q13 = (T)(2.0 * q[1] * q[3]);

		const T q23 = (T)(2.0 * q[2] * q[3]);

		m[0][0] = (T)(1.0l - q11 - q22);
		m[1][0] = (T)(q01 + q23);
		m[2][0] = (T)(q02 - q13);

		m[0][1] = (T)(q01 - q23);
		m[1][1] = (T)(1.0l - q22 - q00);
		m[2][1] = (T)(q12 + q03);

		m[0][2] = (T)(q02 + q13);
		m[1][2] = (T)(q12 - q03);
		m[2][2] = (T)(1.0l - q11 - q00);

		m[0][3] = (T)(0.0l);
		m[1][3] = (T)(0.0l);
		m[2][3] = (T)(0.0l);

		m[3][0] = (T)(0.0l);
		m[3][1] = (T)(0.0l);
		m[3][2] = (T)(0.0l);
		m[3][3] = (T)(1.0l);
	}

	inline void getInverseMatrix(T m[4][4]) const {
		inverse().getMatrix(m);
	}

	inline void getRotationMatrix(T m[3][3]) const {
		static T mat[4][4];
		getMatrix(mat);
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j)
				m[i][j] = mat[j][i]; // Beware of transposition
	}

	inline Matrix3<T> getRotationMatrix() const {
		static T mat[4][4];
		getMatrix(mat);
		return Matrix3<T>(mat[0][0], mat[1][0], mat[2][0],
			mat[0][1], mat[1][1], mat[2][1],
			mat[0][2], mat[1][2], mat[2][2]);
	}

	inline void getInverseRotationMatrix(T m[3][3]) const {
		static T mat[4][4];
		getInverseMatrix(mat);
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j)
				m[i][j] = mat[j][i]; // Beware of transposition
	}

	inline Matrix3<T> getInverseRotationMatrix() const {
		T mat[4][4];
		getInverseMatrix(mat);
		return Matrix3<T>(mat[0][0], mat[1][0], mat[2][0],
			mat[0][1], mat[1][1], mat[2][1],
			mat[0][2], mat[1][2], mat[2][2]);
	}

	inline Matrix4<T> getSE3() {
		T mat[4][4];
		getMatrix(mat);
		return Matrix4<T>(mat[0][0], mat[1][0], mat[2][0], mat[3][0], // col 0
			mat[0][1], mat[1][1], mat[2][1], mat[3][1], // col 1
			mat[0][2], mat[1][2], mat[2][2], mat[3][2], // col 2
			mat[0][3], mat[1][3], mat[2][3], mat[3][3]);// col 3
	}

	inline T* getSE4() const {
		T* mat;
		//getMatrix(mat);
		return mat;
	}

	static Quaternion slerp(const Quaternion& a, const Quaternion& b, T t, bool allowFlip = true);
	static Quaternion squad(const Quaternion& a, const Quaternion& tgA, const Quaternion& tgB, const Quaternion& b, T t);
	static T dot(const Quaternion& a, const Quaternion& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3]; }
	// Returns the "dot" product of \p a and \p b: a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3].

	Quaternion log();
	Quaternion exp();
	static Quaternion lnDif(const Quaternion& a, const Quaternion& b);
	static Quaternion squadTangent(const Quaternion& before, const Quaternion& center, const Quaternion& after);

	static Quaternion randomQuaternion();

	template <typename S> friend std::ostream& operator<<(std::ostream& o, const Quaternion<S> & q);

private:

	T q[4];
};


// OBJ: Constructs a Quaternion that will rotate from the \p from direction to the \p to direction.
//		Note that this rotation is not uniquely defined. The selected axis is usually orthogonal to \p from
//		and \p to. However, this method is robust and can handle small or almost identical vectors.
template <typename T> Quaternion<T>::Quaternion(const Vector3<T>& from, const Vector3<T>& to) {
	const T epsilon = 1e-10f;

	const T fromSqNorm = sqrMag(from);
	const T toSqNorm = sqrMag(to);
	// Identity Quaternion when one vector is null
	if ((fromSqNorm < epsilon) || (toSqNorm < epsilon)) {
		q[0] = q[1] = q[2] = 0.0;
		q[3] = 1.0;
	}
	else {
		Vector3<T> axis = cross(from, to);
		const T axisSqNorm = sqrMag(axis);

		// Aligned vectors, pick any axis, not aligned with from or to
		if (axisSqNorm < epsilon)
			axis = tangential(from);
		//axis = Vector3<T>(1,0,0);

		T angle = asin(sqrt(axisSqNorm / (fromSqNorm * toSqNorm)));

		if (from[0] * to[0] + from[1] * to[1] + from[2] * to[2] < 0.0)
			angle = PI - angle;

		setAxisAngle(axis, angle);
	}
}

// OBJ: Returns the image of \p v by the Quaternion inverse() rotation.
//		rotate() performs an inverse transformation. Same as inverse().rotate(v).
template <typename T> Vector3<T> Quaternion<T>::inverseRotate(const Vector3<T>& v) const {
	return inverse().rotate(v);
}

// OBJ: Returns the image of \p v by the Quaternion rotation.
//		See also inverseRotate() and operator*(const Quaternion&, const Vector3<T>&).
template <typename T> Vector3<T> Quaternion<T>::rotate(const Vector3<T>& v) const {
	const T q00 = (T)(2.0l * q[0] * q[0]);
	const T q11 = (T)(2.0l * q[1] * q[1]);
	const T q22 = (T)(2.0l * q[2] * q[2]);

	const T q01 = (T)(2.0l * q[0] * q[1]);
	const T q02 = (T)(2.0l * q[0] * q[2]);
	const T q03 = (T)(2.0l * q[0] * q[3]);

	const T q12 = (T)(2.0l * q[1] * q[2]);
	const T q13 = (T)(2.0l * q[1] * q[3]);

	const T q23 = (T)(2.0l * q[2] * q[3]);

	return Vector3<T>(
		(T)((1.0 - q11 - q22)*v[0] + (q01 - q23)*v[1] + (q02 + q13)*v[2]),
		(T)((q01 + q23)*v[0] + (1.0 - q22 - q00)*v[1] + (q12 - q03)*v[2]),
		(T)((q02 - q13)*v[0] + (q12 + q03)*v[1] + (1.0 - q11 - q00)*v[2])
		);
}

// OBJ: Set the Quaternion from a (supposedly correct) 3x3 rotation matrix.
//		The matrix is expressed in European format: its three \e columns are the images by the rotation of
//		the three vectors of an orthogonal basis. Note that OpenGL uses a symmetric representation for its
//		matrices.
//		setFromRotatedBasis() sets a Quaternion from the three axis of a rotated frame. It actually fills
//		the three columns of a matrix with these rotated basis vectors and calls this method. */
template <typename T> void Quaternion<T>::setFromRotationMatrix(const T m[3][3]) {
	// Compute one plus the trace of the matrix
	const T onePlusTrace = 1.0 + m[0][0] + m[1][1] + m[2][2];

	if (onePlusTrace > 1E-5) {
		// Direct computation
		const T s = sqrt(onePlusTrace) * 2.0;
		q[0] = (m[2][1] - m[1][2]) / s;
		q[1] = (m[0][2] - m[2][0]) / s;
		q[2] = (m[1][0] - m[0][1]) / s;
		q[3] = 0.25 * s;
	}
	else {
		// Computation depends on major diagonal term
		if ((m[0][0] > m[1][1])&(m[0][0] > m[2][2])) {
			const T s = sqrt(1.0 + m[0][0] - m[1][1] - m[2][2]) * 2.0;
			q[0] = 0.25 * s;
			q[1] = (m[0][1] + m[1][0]) / s;
			q[2] = (m[0][2] + m[2][0]) / s;
			q[3] = (m[1][2] - m[2][1]) / s;
		}
		else if (m[1][1] > m[2][2]) {
			const T s = sqrt(1.0 + m[1][1] - m[0][0] - m[2][2]) * 2.0;
			q[0] = (m[0][1] + m[1][0]) / s;
			q[1] = 0.25 * s;
			q[2] = (m[1][2] + m[2][1]) / s;
			q[3] = (m[0][2] - m[2][0]) / s;
		}
		else {
			const T s = sqrt(1.0 + m[2][2] - m[0][0] - m[1][1]) * 2.0;
			q[0] = (m[0][2] + m[2][0]) / s;
			q[1] = (m[1][2] + m[2][1]) / s;
			q[2] = 0.25 * s;
			q[3] = (m[0][1] - m[1][0]) / s;
		}
	}
	normalize();
}

template <typename T> void Quaternion<T>::setFromSO3(const Matrix3<T>& d) {
	// Compute one plus the trace of the matrix
	const T onePlusTrace = 1.0 + d(0, 0) + d(1, 1) + d(2, 2);

	if (onePlusTrace > 1E-5) {
		// Direct computation
		const T s = sqrt(onePlusTrace) * 2.0;
		q[0] = (d(2, 1) - d(1, 2)) / s;
		q[1] = (d(0, 2) - d(2, 0)) / s;
		q[2] = (d(1, 0) - d(0, 1)) / s;
		q[3] = 0.25 * s;
	}
	else {
		// Computation depends on major diagonal term
		if ((d(0, 0) > d(1, 1))&(d(0, 0) > d(2, 2))) {
			const T s = sqrt(1.0 + d(0, 0) - d(1, 1) - d(2, 2)) * 2.0;
			q[0] = 0.25 * s;
			q[1] = (d(0, 1) + d(1, 0)) / s;
			q[2] = (d(0, 2) + d(2, 0)) / s;
			q[3] = (d(1, 2) - d(2, 1)) / s;
		}
		else if (d(1, 1) > d(2.2)) {
			const T s = sqrt(1.0 + d(1, 1) - d(0, 0) - d(2, 2)) * 2.0;
			q[0] = (d(0, 1) + d(1, 0)) / s;
			q[1] = 0.25 * s;
			q[2] = (d(1, 2) + d(2, 1)) / s;
			q[3] = (d(0, 2) - d(2, 0)) / s;
		}
		else {
			const T s = sqrt(1.0 + d(2, 2) - d(0, 0) - d(1, 1)) * 2.0;
			q[0] = (d(0, 2) + d(2, 0)) / s;
			q[1] = (d(1, 2) + d(2, 1)) / s;
			q[2] = 0.25 * s;
			q[3] = (d(0, 1) - d(1, 0)) / s;
		}
	}
	normalize();
}

template <typename T> void Quaternion<T>::setFromSE3(const Matrix4<T>& d) {
	// Compute one plus the trace of the matrix
	const T onePlusTrace = 1.0 + d(0, 0) + d(1, 1) + d(2, 2);

	if (onePlusTrace > 1E-5) {
		// Direct computation
		const T s = sqrt(onePlusTrace) * 2.0;
		q[0] = (d(2, 1) - d(1, 2)) / s;
		q[1] = (d(0, 2) - d(2, 0)) / s;
		q[2] = (d(1, 0) - d(0, 1)) / s;
		q[3] = 0.25 * s;
	}
	else {
		// Computation depends on major diagonal term
		if ((d(0, 0) > d(1, 1))&(d(0, 0) > d(2, 2))) {
			const T s = sqrt(1.0 + d(0, 0) - d(1, 1) - d(2, 2)) * 2.0;
			q[0] = 0.25 * s;
			q[1] = (d(0, 1) + d(1, 0)) / s;
			q[2] = (d(0, 2) + d(2, 0)) / s;
			q[3] = (d(1, 2) - d(2, 1)) / s;
		}
		else if (d(1, 1) > d(2.2)) {
			const T s = sqrt(1.0 + d(1, 1) - d(0, 0) - d(2, 2)) * 2.0;
			q[0] = (d(0, 1) + d(1, 0)) / s;
			q[1] = 0.25 * s;
			q[2] = (d(1, 2) + d(2, 1)) / s;
			q[3] = (d(0, 2) - d(2, 0)) / s;
		}
		else {
			const T s = sqrt(1.0 + d(2, 2) - d(0, 0) - d(1, 1)) * 2.0;
			q[0] = (d(0, 2) + d(2, 0)) / s;
			q[1] = (d(1, 2) + d(2, 1)) / s;
			q[2] = 0.25 * s;
			q[3] = (d(0, 1) - d(1, 0)) / s;
		}
	}
	normalize();
}

template <typename T> void Quaternion<T>::identity() {
	q[0] = q[1] = q[2] = 0.0;  q[3] = 1.0;
}

// OBJ: Sets the Quaternion from the three rotated vectors of an orthogonal basis.
//		The three vectors do not have to be normalized but must be orthogonal.
//		<code>
//			Quaternion q;
//			q.setFromRotatedBasis(X, Y, Z);
//			Now q.rotate(Vector3(1,0,0)) == X and q.inverseRotate(X) == Vector3(1,0,0)
//			Same goes for Y and Z with Vector3(0,1,0) and Vector3(0,0,1).
//		</code>
//		See also setFromRotationMatrix() and Quaternion(const Vector3<T>&, const Vector3<T>&). */
template <typename T> void Quaternion<T>::setFromRotatedBasis(const Vector3<T>& X, const Vector3<T>& Y, const Vector3<T>& Z) {
	T m[3][3];
	T normX = X.norm();
	T normY = Y.norm();
	T normZ = Z.norm();

	for (int i = 0; i < 3; ++i) {
		m[i][0] = X[i] / normX;
		m[i][1] = Y[i] / normY;
		m[i][2] = Z[i] / normZ;
	}

	setFromRotationMatrix(m);
}

// OBJ: Returns the axis vector and the angle (in radians) of the rotation represented by the Quaternion.
//		See the axis() and angle() documentations.
template <typename T> void Quaternion<T>::getAxisAngle(Vector3<T> &axis, T &angle) const {
	angle = 2.0*acos(q[3]);
	axis = Vector3<T>(q[0], q[1], q[2]);
	const T sinus = mag(axis);
	if (sinus > 1E-8)
		axis /= sinus;

	if (angle > PI) {
		angle = 2.0*PI - angle;
		axis = -axis;
	}
}

// OBJ: Returns the normalized axis direction of the rotation represented by the Quaternion.
//		It is null for an identity Quaternion. See also angle() and getAxisAngle().
template <typename T> Vector3<T> Quaternion<T>::axis() const {
	Vector3<T> res(q[0], q[1], q[2]);
	const T sinus = mag(res);
	if (sinus > 1e-8) res /= sinus;
	return (acos(q[3]) <= PI / 2.0) ? res : -res;
}

// OBJ: Returns the angle (in radians) of the rotation represented by the Quaternion.
//		This value is always in the range [0-pi]. Larger rotational angles are obtained by inverting the
//		axis() direction.
//		See also axis() and getAxisAngle().
template <typename T> T Quaternion<T>::angle() const {
	const T angle = 2.0 * acos(q[3]);
	return (angle <= PI) ? angle : 2.0*PI - angle;
}

// OBJ: Returns the slerp interpolation of Quaternions \p a and \p b, at time \p t.
//		\p t should range in [0,1]. Result is \p a when \p t=0 and \p b when \p t=1.
//		When \p allowFlip is \c true (default) the slerp interpolation will always use the "shortest path"
//		between the Quaternions' orientations, by "flipping" the source Quaternion if needed (see negate()).
template <typename T> Quaternion<T> Quaternion<T>::slerp(const Quaternion& a, const Quaternion& b, T t, bool allowFlip) {
	T cosAngle = Quaternion<T>::dot(a, b);
	T c1, c2;

	// Linear interpolation for close orientations
	if ((1.0 - fabs(cosAngle)) < 0.01) {
		c1 = (T)1.0 - t;
		c2 = t;
	}
	else {
		// Spherical interpolation
		T angle = acos(fabs(cosAngle));
		T sinAngle = sin(angle);
		c1 = (T)sin(angle * (1.0 - t)) / sinAngle;
		c2 = (T)sin(angle * t) / sinAngle;
	}

	// Use the shortest path
	if (allowFlip && (cosAngle < 0.0))
		c1 = -c1;

	return Quaternion<T>(c1*a[0] + c2 * b[0], c1*a[1] + c2 * b[1], c1*a[2] + c2 * b[2], c1*a[3] + c2 * b[3]);
}

// OBJ: Returns the slerp interpolation of the two Quaternions \p a and \p b, at time \p t, using
//		tangents \p tgA and \p tgB.
//		The resulting Quaternion is "between" \p a and \p b (result is \p a when \p t=0 and \p b for \p t=1).
//		Use squadTangent() to define the Quaternion tangents \p tgA and \p tgB. */
template <typename T> Quaternion<T> Quaternion<T>::squad(const Quaternion& a, const Quaternion& tgA, const Quaternion& tgB, const Quaternion& b, T t) {
	Quaternion ab = Quaternion<T>::slerp(a, b, t);
	Quaternion tg = Quaternion<T>::slerp(tgA, tgB, t, false);
	return Quaternion<T>::slerp(ab, tg, 2.0*t*(1.0 - t), false);
}

// OBJ: Returns the logarithm of the Quaternion. See also exp().
template <typename T> Quaternion<T> Quaternion<T>::log() {
	T len = sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2]);

	if (len < 1e-6)
		return Quaternion(q[0], q[1], q[2], 0.0);
	else {
		T coef = acos(q[3]) / len;
		return Quaternion(q[0] * coef, q[1] * coef, q[2] * coef, 0.0);
	}
}

// OBJ: Returns the exponential of the Quaternion. See also log().
template <typename T> Quaternion<T> Quaternion<T>::exp() {
	T theta = sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2]);

	if (theta < 1E-6)
		return Quaternion(q[0], q[1], q[2], cos(theta));
	else {
		T coef = sin(theta) / theta;
		return Quaternion(q[0] * coef, q[1] * coef, q[2] * coef, cos(theta));
	}
}

// OBJ: Returns log(a. inverse() * b). Useful for squadTangent().
template <typename T> Quaternion<T> Quaternion<T>::lnDif(const Quaternion& a, const Quaternion& b) {
	Quaternion dif = a.inverse()*b;
	dif.normalize();
	return dif.log();
}

// OBJ: Returns a tangent Quaternion for \p center, defined by \p before and \p after Quaternions.
//		Useful for smooth spline interpolation of Quaternion with squad() and slerp().
template <typename T> Quaternion<T> Quaternion<T>::squadTangent(const Quaternion& before, const Quaternion& center, const Quaternion& after) {
	Quaternion l1 = Quaternion<T>::lnDif(center, before);
	Quaternion l2 = Quaternion<T>::lnDif(center, after);
	Quaternion e;
	for (int i = 0; i < 4; ++i)
		e.q[i] = -0.25 * (l1.q[i] + l2.q[i]);
	e = center * (e.exp());

	// if (Quaternion<T>::dot(e,b) < 0.0)
	// e.negate();

	return e;
}

// OBJ: Returns a random unit Quaternion.
//		You can create a randomly directed unit vector using:
//		<code>
//			Vec randomDir = Quaternion<T>::randomQuaternion() * Vec(1.0, 0.0, 0.0); // or any other Vec
//		</code>
//		\note This function uses rand() to create pseudo-random numbers and the random number generator can
//		be initialized using srand().
template <typename T> Quaternion<T> Quaternion<T>::randomQuaternion() {
	// The rand() function is not very portable and may not be available on your system.
	// Add the appropriate include or replace by an other random function in case of problem.
	T seed = rand() / (T)RAND_MAX;
	T r1 = sqrt(1.0 - seed);
	T r2 = sqrt(seed);
	T t1 = 2.0 * PI * (rand() / (T)RAND_MAX);
	T t2 = 2.0 * PI * (rand() / (T)RAND_MAX);
	return Quaternion(sin(t1)*r1, cos(t1)*r1, sin(t2)*r2, cos(t2)*r2);
}

template <typename T> std::ostream& operator<<(std::ostream& o, const Quaternion<T> &Q) {
	return o << "[" << Q[0] << ", " << Q[1] << ", " << Q[2] << ", " << Q[3] << "]";
}

typedef Quaternion<float>	Quaternionf;
typedef Quaternion<double>	Quaterniond;

#endif // __QUATERNION_H__

