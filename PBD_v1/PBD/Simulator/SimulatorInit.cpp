#include "ClothSimulator.h"

void ClothSimulator::initWithClothMesh(TriangularMesh* clmesh) {
	clothMesh = clmesh;

	nV = clothMesh->getNumOfVertices();
	nE = clothMesh->getNumOfEdges();
	nT = clothMesh->getNumOfTriangles();
	nB = clothMesh->getNumOfTrianglePairs();

	pos = clothMesh->vertices;

	triangles.resize(nT);
	for (int i = 0; i < nT; i++)
		triangles[i] = clothMesh->triangles[i].vindices;

	edges = clothMesh->edges;
	bendingElementIdx = clothMesh->adjacentTrianglesIdx;

	trianglesEdges = clothMesh->trianglesEdges;

	if (clothMesh->texcoords.empty())
		return;

	restPos.resize(nV);
	for (int i = 0; i < nT; i++) {
		restPos[clothMesh->triangles[i].vindices[0]] = clothMesh->texcoords[clothMesh->triangles[i].tindices[0]];
		restPos[clothMesh->triangles[i].vindices[1]] = clothMesh->texcoords[clothMesh->triangles[i].tindices[1]];
		restPos[clothMesh->triangles[i].vindices[2]] = clothMesh->texcoords[clothMesh->triangles[i].tindices[2]];
	}

	restLengths.resize(nE);
	for (int i = 0; i < nE; i++)
		restLengths[i] = (restPos[edges[i][0]] - restPos[edges[i][1]]).norm();

	restTriangleAreas.resize(nT);
	for (int i = 0; i < nT; i++) {
		double area;
		double a = (restPos[triangles[i][0]] - restPos[triangles[i][1]]).norm();
		double b = (restPos[triangles[i][0]] - restPos[triangles[i][2]]).norm();
		double c = (restPos[triangles[i][1]] - restPos[triangles[i][2]]).norm();

		double l = (a + b + c) / 2.0f;
		area = sqrt(l*(l - a)*(l - b)*(l - c));
		restTriangleAreas[i] = area;
	}

	force.resize(nV);
	prev_pos.resize(nV);
	pred_pos.resize(nV);
	vel.resize(nV);
	prev_vel.resize(nV);

	// Mass
	mass.resize(nV);
	inv_mass.resize(nV);
	double density = Parameter::getInstance().getMass();
	if (nV > 0) memset(&(mass[0]), 0, sizeof(double)*nV);
	for (int i = 0; i < nT; ++i) {
		double areaDensity = restTriangleAreas[i] * (density) / 3.0;
		for (int j = 0; j < 3; j++) {
			mass[triangles[i][j]] += areaDensity;
			inv_mass[triangles[i][j]] = 1 / mass[triangles[i][j]];
		}
	}

	updateFacetNormals();
	updateVertexNormals();
	initRestAngles();
	initRestLenghtsOfBending();

	CCD_colliders.push_back(CCD_Object());
}

void ClothSimulator::updateFacetNormals() {
	int  i;
	Vector3d u;
	Vector3d v;

	facetNormals.clear();
	facetNormals.resize(nT);

#pragma omp parallel for
	for (i = 0; i < nT; i++) {
		Vector3d a = pos[triangles[i][0]];
		Vector3d b = pos[triangles[i][1]];
		Vector3d c = pos[triangles[i][2]];

		u = b - a;
		v = c - a;

		facetNormals[i] = cross(u, v);
		normalize(facetNormals[i]);
	}
}

void ClothSimulator::updateVertexNormals() {
	vertexNormals.clear();
	vertexNormals.resize(nV);

#pragma omp parallel for
	for (int i = 0; i < nT; ++i) {
		Vector3d n = facetNormals[i];
		vertexNormals[triangles[i][0]] += n;
		vertexNormals[triangles[i][1]] += n;
		vertexNormals[triangles[i][2]] += n;
	}

#pragma omp parallel for
	for (int i = 0; i < nV; ++i) {
		normalize(vertexNormals[i]);
	}
}

void ClothSimulator::reset() {
	if (!clothMesh)
		return;

	pos = clothMesh->vertices;

	force.resize(nV);
	prev_pos.resize(nV);
	pred_pos.resize(nV);
	vel.resize(nV);
	prev_vel.resize(nV);

	updateFacetNormals();
	updateVertexNormals();
}