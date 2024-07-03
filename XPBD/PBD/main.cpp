#include <GL/freeglut.h>

#include "../Scene/Object.h"
#include "../Scene/Light.h"
#include "../Scene/Camera.h"
#include "../Simulator/ClothSimulator.h"
#include "../Simulator/TriangularMesh.h"

const int width = 1024, height = 1024;

CLight light;
CCamera camera;
ClothSimulator clothSim;

void init() {
	light.setPosition(100, 100, 100);
	clothSim.initPBD();
}

void idle() {
	clothSim.updatePBD();
	glutPostRedisplay();
}

void renderScene() {
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(.9f, .9f, .9f, 1.f);

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_LIGHTING);

	light.initialize();
	camera.initialize();

	camera.glRender();
	light.glRender();

	clothSim.drawMesh();
	clothSim.drawVertices();

	glutSwapBuffers();
}

void main(int argc, char** argv) {
	string path = "./OBJ/cl.obj";
	TriangularMesh* tm = new TriangularMesh(path.c_str());
	clothSim.initWithClothMesh(tm);
	init();

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
	glutInitWindowPosition(500, 50);
	glutInitWindowSize(width, height);
	glutCreateWindow("Position based Dynamics");

	glutDisplayFunc(renderScene);
	glutIdleFunc(idle);

	glutMainLoop();
}