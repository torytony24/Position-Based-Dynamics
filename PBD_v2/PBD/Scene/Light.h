#pragma once

#include <GL/freeglut.h>

#include "../Math\vector.h"
#include "../Math\quaternion.h"

#define DEFALT_AMBIENT			0.1f
#define DEFALT_DIFFUSE			0.7f
#define DEFALT_SPECULAR			0.8f


class CLight {
public:
	// Constructors & Destructors
	inline CLight();
	inline CLight(GLenum _ID);
	inline ~CLight();

	inline void            initialize();

	// Getters & setters
	inline GLenum		   getID();
	inline void			   setID(GLenum _ID);
	inline float*		   getAmbient();
	inline void			   setAmbient(float r, float g, float b, float a);
	inline void			   setAmbient(float* _src);
	inline float*		   getDiffuse();
	inline void			   setDiffuse(float r, float g, float b, float a);
	inline void			   setDiffuse(float* _src);
	inline float*		   getSpecular();
	inline void			   setSpecular(float r, float g, float b, float a);
	inline void			   setSpecular(float* _src);
	inline Vector3f        getPosition();
	inline void            setPosition(Vector3f _src);
	inline void            setPosition(float x, float y, float z);
	inline Quaternionf     getOrientation();
	inline void            setOrientation(Quaternionf _src);
	inline void            setOrientation(float x, float y, float z, float theta);
	inline void            setOrientation(Vector3f axis, float theta);
	inline void			   setLightAttrib(float amb_r, float amb_g, float amb_b, float amb_a,
							float dif_r, float dif_g, float dif_b, float dif_a,
							float spe_r, float spe_g, float spe_b, float spe_a,
							float matSpe_r, float matSpe_g, float matSpe_b, float matSpe_a,
							float shi);
	inline void			   setLightAttrib(float* amb, float* dif, float* spe, float* matSpe, float shi);

	// Member Functions
	inline void            glRender();
	inline void		       glDrawLight();
	inline void			   on();
	inline void			   off();
	inline void			   toggle();
	//inline void            translate(float x, float y, float z);
	//inline void            translate(Vector3f vec);
	//inline void            rotate(Vector3f axis, float theta);
	//inline void            rotate(Quaternionf rotation);

	// Private Member Variables
private:
	// Attribute
	GLfloat ambient[4];
	GLfloat diffuse[4];
	GLfloat specular[4];
	Vector3f position;
	Quaternionf orientation;

	GLenum ID;
	bool state;
};

// Constructors & Destructors
inline CLight::CLight() {
	ambient[0] = ambient[1] = ambient[2] = DEFALT_AMBIENT;
	diffuse[0] = diffuse[1] = diffuse[2] = DEFALT_DIFFUSE;
	specular[0] = specular[1] = specular[2] = DEFALT_SPECULAR;
	ambient[3] = diffuse[3] = specular[3] = 1.f;

	position = Vector3f(0.f, 0.f, 10.f);
	orientation = Quaternionf(0.f, 0.f, 0.f, 1.f);
	ID = GL_LIGHT0;
	state = false;
}

inline CLight::CLight(GLenum _ID) {
	ambient[0] = ambient[1] = ambient[2] = DEFALT_AMBIENT;
	diffuse[0] = diffuse[1] = diffuse[2] = DEFALT_DIFFUSE;
	specular[0] = specular[1] = specular[2] = DEFALT_SPECULAR;
	ambient[3] = diffuse[3] = specular[3] = 1.f;

	position = Vector3f(0.f, 0.f, 10.f);
	orientation = Quaternionf(0.f, 0.f, 0.f, 1.f);
	ID = _ID;
	state = false;
}

inline CLight::~CLight() {

}

inline void CLight::initialize() {
	glRender();
	glEnable(ID);
}

// Getters & setters
inline GLenum CLight::getID() {
	return ID;
}

inline void CLight::setID(GLenum _ID) {
	ID = _ID;
	glRender();
}

inline float* CLight::getAmbient() {
	return ambient;
}

inline void	CLight::setAmbient(float r, float g, float b, float a) {
	ambient[0] = r;
	ambient[1] = g;
	ambient[2] = b;
	ambient[3] = a;
	glLightfv(ID, GL_AMBIENT, ambient);
}

inline void	CLight::setAmbient(float* _src) {
	ambient[0] = _src[0];
	ambient[1] = _src[1];
	ambient[2] = _src[2];
	ambient[3] = _src[3];
	glLightfv(ID, GL_AMBIENT, ambient);
}

inline float* CLight::getDiffuse() {
	return diffuse;
}

inline void	CLight::setDiffuse(float r, float g, float b, float a) {
	diffuse[0] = r;
	diffuse[1] = g;
	diffuse[2] = b;
	diffuse[3] = a;
	glLightfv(ID, GL_DIFFUSE, diffuse);
}

inline void	CLight::setDiffuse(float* _src) {
	diffuse[0] = _src[0];
	diffuse[1] = _src[1];
	diffuse[2] = _src[2];
	diffuse[3] = _src[3];
	glLightfv(ID, GL_DIFFUSE, diffuse);
}

inline float* CLight::getSpecular() {
	return specular;
}

inline void	CLight::setSpecular(float r, float g, float b, float a) {
	specular[0] = r;
	specular[1] = g;
	specular[2] = b;
	specular[3] = a;
	glLightfv(ID, GL_SPECULAR, specular);
}

inline void	CLight::setSpecular(float* _src) {
	specular[0] = _src[0];
	specular[1] = _src[1];
	specular[2] = _src[2];
	specular[3] = _src[3];
	glLightfv(ID, GL_SPECULAR, specular);
}

inline Vector3f CLight::getPosition() {
	return position;
}

inline void CLight::setPosition(Vector3f _src) {
	position = _src;
	GLfloat pos[4] = { position.x, position.y, position.z, 1 };
	glLightfv(ID, GL_POSITION, pos);
}

inline void CLight::setPosition(float x, float y, float z) {
	position[0] = x;
	position[1] = y;
	position[2] = z;
	GLfloat pos[4] = { position.x, position.y, position.z, 1 };
	glLightfv(ID, GL_POSITION, pos);
}

inline Quaternionf CLight::getOrientation() {
	return orientation;
}

inline void CLight::setOrientation(Quaternionf _src) {
	orientation = _src;
}

inline void CLight::setOrientation(float x, float y, float z, float theta) {
	Vector3f axis = Vector3f(x, y, z);
	normalize(axis);
	orientation = Quaternionf(axis, theta);
}

inline void CLight::setOrientation(Vector3f axis, float theta) {
	orientation = Quaternionf(axis, theta);
}

inline void	CLight::setLightAttrib(float amb_r, float amb_g, float amb_b, float amb_a,
	float dif_r, float dif_g, float dif_b, float dif_a,
	float spe_r, float spe_g, float spe_b, float spe_a,
	float matSpe_r, float matSpe_g, float matSpe_b, float matSpe_a,
	float shi) {
	setAmbient(amb_r, amb_g, amb_b, amb_a);
	setDiffuse(dif_r, dif_g, dif_b, dif_a);
	setSpecular(spe_r, spe_g, spe_b, spe_a);
}

inline void	CLight::setLightAttrib(float* amb, float* dif, float* spe, float* matSpe, float shi) {
	setAmbient(amb);
	setDiffuse(dif);
	setSpecular(spe);
}

// Member Functions
inline void CLight::glRender() {
	glLightfv(ID, GL_AMBIENT, ambient);
	glLightfv(ID, GL_DIFFUSE, diffuse);
	glLightfv(ID, GL_SPECULAR, specular);
	GLfloat pos[4] = { position.x, position.y, position.z, 1 };
	glLightfv(ID, GL_POSITION, pos);
	glEnable(ID);
}

inline void CLight::glDrawLight() {
	glDisable(GL_LIGHTING);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
	glTranslatef(position.x, position.y, position.z);
	glColor3f(0.f, 1.f, 0.f);
	glutSolidSphere(5.f, 10, 10);
	glPopMatrix();
	glEnable(GL_LIGHTING);
}

inline void	CLight::on() {
	glEnable(ID);
	state = true;
}

inline void	CLight::off() {
	glDisable(ID);
	state = false;
}

inline void	CLight::toggle() {
	state ? off() : on();
}

//inline void CLight::translate(float x, float y, float z) {
//	Vector3f trans(x, y, z);
//	normalize(trans);
//	position += TRANSLATE_SENSITIVITY * trans;
//	lookAt += TRANSLATE_SENSITIVITY * trans;
//}
//
//inline void CLight::translate(Vector3f vec) {
//	normalize(vec);
//	position += TRANSLATE_SENSITIVITY * vec;
//	lookAt += TRANSLATE_SENSITIVITY * vec;
//}
//
//inline void CLight::rotate(Vector3f axis, float theta) {
//	Quaternionf rotateQuat(axis, -theta);
//
//	Vector3f lookAtVec = lookAt - position;
//	lookAtVec = rotateQuat.rotate(lookAtVec);
//	position = lookAt - lookAtVec;
//	upVec = rotateQuat.rotate(upVec);
//	orientation = rotateQuat * orientation;
//}
//
//inline void CLight::rotate(Quaternionf rotation) {
//	Vector3f lookAtVec = lookAt - position;
//	lookAtVec = rotation.rotate(lookAtVec);
//	position = lookAt - lookAtVec;
//	upVec = rotation.rotate(upVec);
//	orientation = rotation * orientation;
//}
