#pragma once

#include <list>

#include "TriangularMesh.h"
#include "../layout/Parameter.h"
#include "../Math/matrix.h"
#include "../Math/performance.h"

#include "../OpenCCD/BVH.h"
#include "../OpenCCD/CCD.h"

class ClothSimulator {
public:
	ClothSimulator();
	ClothSimulator(TriangularMesh* clmesh);
	~ClothSimulator();

	void initWithClothMesh(TriangularMesh* clmesh);
	void reset();
	void toggleFixedConstraint(int index);

	void addCollider(TriangularMesh* tm);

	// draw
	void drawMesh();
	void drawVertices();
	void drawBVH();

	// for PBD
	void initPBD();
	void updatePBD();

private:
	// Original Mesh
	TriangularMesh*				clothMesh;

	vector<TriangularMesh*>		colliders;
	vector<CCD_Object>			CCD_colliders;

	// Geometry
	unsigned int				nV;
	unsigned int				nE;
	unsigned int				nT;
	unsigned int				nB;

	vector<Vector3i>			triangles;
	vector<Vector2i>			edges;
	vector<Vector3i>			trianglesEdges;
	vector<int>					bendingElementIdx;

	vector<Vector2d>			restPos;
	vector<Vector3i>			restTri;
	vector<double>				restLengths;
	vector<double>				restTriangleAreas;

	vector<Vector3d>			vertexNormals;
	vector<Vector3d>			facetNormals;

	// Physical Datas
	vector<double>				mass;
	vector<Vector3d>			pos;
	vector<Vector3d>			prev_pos;
	vector<Vector3d>			vel;
	vector<Vector3d>			prev_vel;
	vector<Vector3d>			force;

	list<int>					fixedConstraintVertices;

	void updateFacetNormals();
	void updateVertexNormals();

	// for PBD
	vector<double>				inv_mass;
	vector<Vector3d>			pred_pos;	// prediction
	vector<double>				restAngles;
	vector<double>				restLengthsOfBending;
	void initRestAngles();
	void initRestLenghtsOfBending();
	void calForce();
	void explicitEuler();
	void calDistanceConstraint();
	void calBendingConstraint();
	void integrationForPBD();

};
