#include <GL/freeglut.h>

#include "ClothSimulator.h"

ClothSimulator::ClothSimulator() {
	clothMesh = NULL;
	nV = nE = nT = nB = 0;
}

ClothSimulator::ClothSimulator(TriangularMesh* clmesh) :clothMesh(clmesh) {
	initWithClothMesh(clmesh);
}

ClothSimulator::~ClothSimulator() {

}

void ClothSimulator::toggleFixedConstraint(int index) {
	list<int>::iterator itr = fixedConstraintVertices.begin();
	for (; itr != fixedConstraintVertices.end(); itr++) {
		if (*itr < index)
			continue;
		else if (*itr > index) {
			fixedConstraintVertices.insert(itr, index);
			return;
		}
		else if (index == *itr) {
			fixedConstraintVertices.erase(itr);
			return;
		}
	}

	fixedConstraintVertices.push_back(index);
}

void ClothSimulator::addCollider(TriangularMesh* tm) {
	colliders.push_back(tm);
}

void ClothSimulator::drawMesh() {
	glPushMatrix();

	// Material setting (pearl)
	float mat_emission[] = { 0.0, 0.0, 0.0, 1.0 };
	float mat_ambient[] = { 0.25, 0.20725, 0.20725, 0.922 };
	float mat_diffuse[] = { 1.0, 0.829, 0.829, 0.922 };
	float mat_specular[] = { 0.296648, 0.296648, 0.296648, 0.922 };
	float mat_shininess[] = { 11.264 };

	glShadeModel(GL_SMOOTH);
	glMaterialfv(GL_FRONT, GL_EMISSION, mat_emission);
	glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
	glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);

	glBegin(GL_TRIANGLES);
	for (int i = 0; i < nT; i++) {
		//glColor3f(1, 0, 0);
		glVertex3f(pos[triangles[i][0]][0], pos[triangles[i][0]][1], pos[triangles[i][0]][2]);
		glVertex3f(pos[triangles[i][1]][0], pos[triangles[i][1]][1], pos[triangles[i][1]][2]);
		glVertex3f(pos[triangles[i][2]][0], pos[triangles[i][2]][1], pos[triangles[i][2]][2]);
	}
	glEnd();

	glPopMatrix();
}

void ClothSimulator::drawVertices() {
	glPushMatrix();

	glDisable(GL_LIGHTING);

	// vertices
	glPointSize(8.0);
	glBegin(GL_POINTS);
	glColor3f(1, 0.3, 0.3);
	for (int i = 0; i < nV; i++)
		glVertex3f(pos[i][0], pos[i][1], pos[i][2]);
	glEnd();

	// edges
	glLineWidth(4.0);
	glBegin(GL_LINES);
	glColor3f(0.5, 0, 0);
	for (int i = 0; i < nE; i++) {
		glVertex3f(pos[edges[i][0]][0], pos[edges[i][0]][1], pos[edges[i][0]][2]);
		glVertex3f(pos[edges[i][1]][0], pos[edges[i][1]][1], pos[edges[i][1]][2]);
	}
	glEnd();
	//glEnable(GL_LIGHTING);

	glPopMatrix();
}

void ClothSimulator::drawBVH() {
	// argument meaning (1: root, 2: children of root, ...)
	for (auto col : CCD_colliders)
		if(col.getBVH() != NULL)
			col.visualizeBVH(3);
}

void ClothSimulator::initPBD() {
	calForce();
}

double itr = Parameter::getInstance().getIteration();

void ClothSimulator::updatePBD() {
	explicitEuler();

	for (int n = 0; n < itr; n++) {
		calDistanceConstraint();
		calBendingConstraint();
		calVolumeConstraint();
	}
	
	integrationForPBD();
}

void ClothSimulator::initRestAngles() {
	restAngles.resize(nE);
	fill(restAngles.begin(), restAngles.end(), 3.141592);
}

void ClothSimulator::initRestLenghtsOfBending() {
	for (int i = 0; i < nE; i++) {
		Vector3d p1 = pos[edges[i][0]];
		Vector3d p2 = pos[edges[i][1]];
		double l0 = (p1 - p2).norm();
		restLengthsOfBending.push_back(l0);
	}
}

void ClothSimulator::initRestVolume() {
	for (int i = 0; i < nT; i++) {
		Vector3d p1 = pos[triangles[i][0]];
		Vector3d p2 = pos[triangles[i][1]];
		Vector3d p3 = pos[triangles[i][2]];
		restVolume += dot(cross(p1, p2), p3);
	}
}

void ClothSimulator::initTotalMass() {
	for (int i = 0; i < nV; i++) {
		totalMass += inv_mass[i];
	}
}

void ClothSimulator::initVertexTriangles()
{
	vertexTriangles.resize(nV);

	for (int i = 0; i < nT; i++) {
		int idx1 = triangles[i][0];
		int idx2 = triangles[i][1];
		int idx3 = triangles[i][2];

		vertexTriangles[idx1].push_back(i);
		vertexTriangles[idx2].push_back(i);
		vertexTriangles[idx3].push_back(i);
	}
}

void ClothSimulator::calForce() {
	double g = Parameter::getInstance().getGravity();
	for (int i = 0; i < nV; i++) {
		force[i] = { 0.0, mass[i] * g, 0.0 };
	}
}

void ClothSimulator::calDistanceConstraint() {
	for (int i = 0; i < nE; i++) {
		Vector3d p1 = pred_pos[edges[i][0]];
		Vector3d p2 = pred_pos[edges[i][1]];
		Vector3d n = (p1 - p2) / (p1 - p2).norm();

		double w1 = inv_mass[edges[i][0]];
		double w2 = inv_mass[edges[i][1]];
		double C = (p1 - p2).norm() - restLengthsOfBending[i];

		Vector3d dp1 = -(w1 / (w1 + w2)) * C * n;
		Vector3d dp2 = (w2 / (w1 + w2)) * C * n;

		pred_pos[edges[i][0]] += dp1;
		pred_pos[edges[i][1]] += dp2;
	}
}

void ClothSimulator::calBendingConstraint() {
	for (int i = 0; i < bendingElementIdx.size() / 6; i++) {
		int idx1 = bendingElementIdx[6 * i + 2];
		int idx2 = bendingElementIdx[6 * i + 3];
		int idx3 = bendingElementIdx[6 * i + 0];
		int idx4 = bendingElementIdx[6 * i + 1];

		Vector3d p2 = pred_pos[idx2] - pred_pos[idx1];
		Vector3d p3 = pred_pos[idx3] - pred_pos[idx1];
		Vector3d p4 = pred_pos[idx4] - pred_pos[idx1];
		Vector3d n1 = cross(p2, p3) / cross(p2, p3).norm();
		Vector3d n2 = cross(p2, p4) / cross(p2, p4).norm();
		double d = dot(n1, n2);
		if (d < -1) d = -1;
		else if (d > 1)d = 1;

		Vector3d q3 = (cross(p2, n2) + cross(n1, p2) * d) / cross(p2, p3).norm();
		Vector3d q4 = (cross(p2, n1) + cross(n2, p2) * d) / cross(p2, p4).norm();
		Vector3d q2 = -(cross(p3, n2) + cross(n1, p3) * d) / cross(p2, p3).norm() - (cross(p4, n1) + cross(n2, p4) * d) / cross(p2, p4).norm();
		Vector3d q1 = -q2 - q3 - q4;

		double w1 = inv_mass[idx1];
		double w2 = inv_mass[idx2];
		double w3 = inv_mass[idx3];
		double w4 = inv_mass[idx4];

		double s = 4 / (w1 + w2 + w3 + w4);
		double dp = -s * sqrt(1 - d * d) * (acos(d) - restAngles[i]) / (q1.norm() * q1.norm() + q2.norm() * q2.norm() + q3.norm() * q3.norm() + q4.norm() * q4.norm());

		pred_pos[idx1] += dp * w1 * q1;
		pred_pos[idx2] += dp * w2 * q2;
		pred_pos[idx3] += dp * w3 * q3;
		pred_pos[idx4] += dp * w4 * q4;

	}
}

void ClothSimulator::calVolumeConstraint() {

	// calculate all gradC_pi
	double gradSum = 0;
	vector<Vector3d> gradC_pi;

	for (int i = 0; i < nV; i++) {
		Vector3d gradC;
		vector<int> adjTriangles = vertexTriangles[i];
		for (int j : adjTriangles) {
			int idx1 = triangles[j][0];
			int idx2 = triangles[j][1];
			int idx3 = triangles[j][2];
			Vector3d p1 = pred_pos[idx1];
			Vector3d p2 = pred_pos[idx2];
			Vector3d p3 = pred_pos[idx3];
			if (idx1 == i) {
				gradC += cross(p2, p3);
			}
			else if (idx2 == i) {
				gradC += cross(p3, p1);
			}
			else if (idx3 == i) {
				gradC += cross(p1, p2);
			}
		}
		gradC_pi.push_back(gradC);
		gradSum += gradC.norm() * gradC.norm();
	}
	

	// calculating C
	double C = 0;
	for (int j = 0; j < nT; j++) {
		Vector3d p1 = pred_pos[triangles[j][0]];
		Vector3d p2 = pred_pos[triangles[j][1]];
		Vector3d p3 = pred_pos[triangles[j][2]];
		C += dot(cross(p1, p2), p3);
	}
	C -= restVolume;

	// calculating dp
	for (int i = 0; i < nV; i++) {
		Vector3d dp = -(C / gradSum) * (nV * inv_mass[i] / totalMass) * gradC_pi[i];
		pred_pos[i] += dp;
	}
}