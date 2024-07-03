#pragma once

#include <gl/glut.h>
#include "../Math/Util.h"
#include "../Math/vector.h"
#include "../Math/quaternion.h"

#define TRANSLATE_SENSITIVITY   0.25f
#define DOLLY_SENSITIVITY       0.5f
#define DOLLY_MINIMUM           1.f
#define DOLLY_MAXIMUM           100.f
#define ROTATE_SENSITIVITY      1.f
#define ZOOMING_SENSITIVITY     3.f
#define ZOOMING_MINIMUM         120.f
#define ZOOMING_MAXIMUM         30.f

class CCamera {
public:
	// Constructors & Destructors
	inline  CCamera();
	inline  ~CCamera();

	inline void            initialize();

	// Getters & setters
	inline Vector3f        getPosition();
	inline void            setPosition(Vector3f _src);
	inline void            setPosition(float x, float y, float z);
	inline Vector3f        getLookAt();
	inline void            setLookAt(Vector3f _src);
	inline void            setLookAt(float x, float y, float z);
	inline Vector3f        getUpVec();
	inline void            setUpVec(Vector3f _src);
	inline void            setUpVec(float x, float y, float z);
	inline Quaternionf     getOrientation();
	inline void            setOrientation(Quaternionf _src);
	inline void            setOrientation(float x, float y, float z, float theta);
	inline void            setOrientation(Vector3f axis, float theta);
	inline float           getNear();
	inline void            setNear(float _near);
	inline float           getFar();
	inline void            setFar(float _far);
	inline float           getViewingAngle();
	inline void            setViewingAngle(float _vangle);
	inline Vector3f        getLocalXAxis();
	inline Vector3f        getLocalYAxis();
	inline Vector3f        getLocalZAxis();

	// Member Functions
	inline void            glRender();
	inline void            translate(float x, float y, float z);
	inline void			   translate(Vector3f vec);
	inline void			   translate(int prevX, int prevY, int curX, int curY);
	inline void            rotate(Vector3f axis, float theta);
	inline void            rotate(Quaternionf rotation);
	inline void			   rotate(int prevX, int prevY, int curX, int curY);
	inline void            zoomIn();
	inline void            zoomOut();
	inline void            dollyIn();
	inline void            dollyOut();
	inline Vector3f		   viewPortPositionToTrackBall(int x, int y);

	// Private Member Variables
private:
	// Projection
	float zNear;
	float zFar;
	float viewingAngle;

	// ModelView
	Vector3f position;
	Vector3f lookAt;
	Vector3f upVec;
	Quaternionf orientation;
};

// Constructors & Destructors
inline CCamera::CCamera()
	:position(Vector3f(0.f, 0.f, 30.f)), lookAt(Vector3f(0.f, 0.f, 0.f)), upVec(Vector3f(0.f, 1.f, 0.f)), zNear(0.1f), zFar(1000.f), viewingAngle(45.f), orientation(Quaternionf(0.f, 0.f, 0.f, 1.f)) {

}

inline CCamera::~CCamera() {

}

inline void CCamera::initialize() {
	position = Vector3f(0.f, 50.f, 80.f);
	lookAt = Vector3f(0.f, 0.f, 0.f);
	upVec = Vector3f(0.f, 1.f, 0.f);
	zNear = 0.1f;
	zFar = 1000.f;
	viewingAngle = 45.f;
	orientation = Quaternionf(0.f, 0.f, 0.f, 1.f);
}

// Getters & setters
inline Vector3f CCamera::getPosition() {
	return position;
}

inline void CCamera::setPosition(Vector3f _src) {
	position = _src;
}

inline void CCamera::setPosition(float x, float y, float z) {
	position[0] = x;
	position[1] = y;
	position[2] = z;
}

inline Vector3f CCamera::getLookAt() {
	return lookAt;
}

inline void CCamera::setLookAt(Vector3f _src) {
	lookAt = _src;
}

inline void CCamera::setLookAt(float x, float y, float z) {
	lookAt[0] = x;
	lookAt[1] = y;
	lookAt[2] = z;
}

inline Vector3f CCamera::getUpVec() {
	return upVec;
}

inline void CCamera::setUpVec(Vector3f _src) {
	upVec = _src;
}

inline void CCamera::setUpVec(float x, float y, float z) {
	upVec[0] = x;
	upVec[1] = y;
	upVec[2] = z;
}

inline Quaternionf CCamera::getOrientation() {
	return orientation;
}

inline void CCamera::setOrientation(Quaternionf _src) {
	orientation = _src;
}

inline void CCamera::setOrientation(float x, float y, float z, float theta) {
	Vector3f axis = Vector3f(x, y, z);
	normalize(axis);
	orientation = Quaternionf(axis, theta);
}

inline void CCamera::setOrientation(Vector3f axis, float theta) {
	orientation = Quaternionf(axis, theta);
}

inline float CCamera::getNear() {
	return zNear;
}

inline void CCamera::setNear(float _near) {
	zNear = _near;
}

inline float CCamera::getFar() {
	return zFar;
}

inline void CCamera::setFar(float _far) {
	zFar = _far;
}

inline float CCamera::getViewingAngle() {
	return viewingAngle;
}

inline void CCamera::setViewingAngle(float _vangle) {
	viewingAngle = _vangle;
}

inline Vector3f CCamera::getLocalXAxis() {
	return normalized(cross(getLocalYAxis(), getLocalZAxis()));
}

inline Vector3f CCamera::getLocalYAxis() {
	return normalized(getUpVec());
}

inline Vector3f CCamera::getLocalZAxis() {
	return normalized(getPosition() - getLookAt());
}

// Member Functions
inline void CCamera::glRender() {
	GLint m_viewport[4];
	glGetIntegerv(GL_VIEWPORT, m_viewport);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	gluPerspective(viewingAngle, (float)m_viewport[2] / (float)m_viewport[3], zNear, zFar);
	gluLookAt(position[0], position[1], position[2], lookAt[0], lookAt[1], lookAt[2], upVec[0], upVec[1], upVec[2]);
}

inline void CCamera::translate(float x, float y, float z) {
	Vector3f trans(x, y, z);
	normalize(trans);
	position += TRANSLATE_SENSITIVITY * trans;
	lookAt += TRANSLATE_SENSITIVITY * trans;
}

inline void CCamera::translate(Vector3f vec) {
	normalize(vec);
	position += TRANSLATE_SENSITIVITY * vec;
	lookAt += TRANSLATE_SENSITIVITY * vec;
}

inline void CCamera::translate(int prevX, int prevY, int curX, int curY) {
	float dx = curX - prevX;
	float dy = curY - prevY;
	translate(-dx * getLocalXAxis() + dy * getLocalYAxis());
}

inline void CCamera::rotate(Vector3f axis, float theta) {
	Quaternionf rotateQuat(axis, -theta);

	Vector3f lookAtVec = lookAt - position;
	lookAtVec = rotateQuat.rotate(lookAtVec);
	position = lookAt - lookAtVec;
	upVec = rotateQuat.rotate(upVec);
	orientation = rotateQuat * orientation;
}

inline void CCamera::rotate(Quaternionf rotation) {
	Vector3f lookAtVec = lookAt - position;
	lookAtVec = rotation.rotate(lookAtVec);
	position = lookAt - lookAtVec;
	upVec = rotation.rotate(upVec);
	orientation = rotation * orientation;
}

inline void CCamera::rotate(int prevX, int prevY, int curX, int curY) {
	Vector3f prevPos = viewPortPositionToTrackBall(prevX, prevY);
	Vector3f currentPos = viewPortPositionToTrackBall(curX, curY);
	Quaternionf rotation(currentPos, prevPos);
	rotate(rotation);
}

inline void CCamera::zoomIn() {
	if (viewingAngle - ZOOMING_SENSITIVITY > ZOOMING_MAXIMUM)
		viewingAngle -= ZOOMING_SENSITIVITY;
}

inline void CCamera::zoomOut() {
	if (viewingAngle + ZOOMING_SENSITIVITY < ZOOMING_MINIMUM)
		viewingAngle += ZOOMING_SENSITIVITY;
}

inline void CCamera::dollyIn() {
	Vector3f vec = lookAt - position;
	normalize(vec);
	if ((lookAt - position).norm() > DOLLY_MINIMUM)
		position += DOLLY_SENSITIVITY * vec;
	//lookAt += DOLLY_SENSITIVITY * vec;
}

inline void CCamera::dollyOut() {
	Vector3f vec = position - lookAt;
	normalize(vec);
	if ((lookAt - position).norm() < DOLLY_MAXIMUM)
		position += DOLLY_SENSITIVITY * vec;
	//lookAt += DOLLY_SENSITIVITY * vec;
}

Vector3f CCamera::viewPortPositionToTrackBall(int x, int y) {
	// 2D - Viewport -> 3D - sphere
	GLint m_viewport[4];

	glGetIntegerv(GL_VIEWPORT, m_viewport);

	Vector3f tempVec;

	// Translate of Origin Point
	tempVec[0] = (float)x - m_viewport[2] / 2.f;
	tempVec[1] = (float)y - m_viewport[3] / 2.f;

	// Scaling
	tempVec[0] /= (float)m_viewport[2] / SQRT2;
	tempVec[1] /= (float)(-m_viewport[3]) / SQRT2;
	tempVec[2] = 0.f;

	// Normalize
	if (tempVec.norm() >= 1.f)
		normalize(tempVec);
	else
		tempVec[2] = sqrt(1.f - tempVec[0] * tempVec[0] - tempVec[1] * tempVec[1]);

	// Consider Orientation of Trackball
	tempVec = getOrientation().rotate(tempVec);

	return tempVec;
}