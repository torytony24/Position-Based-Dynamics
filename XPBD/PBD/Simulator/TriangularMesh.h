#pragma once

#include "../Scene/Object.h"

class TriangularMesh : public Object {
public:
	TriangularMesh();
	TriangularMesh(const char* filePath);
	~TriangularMesh();

	void clear() {};

	void initGeometryInformation();
	void initEdges();
	void initAdjacentTrianglesElements();
	void printMeshInfo();

	int getNumOfEdges() { return edges.size(); };
	int getNumOfTrianglePairs() { return adjacentTrianglesIdx.size() / 6; };

	// Additional Mesh Information
	vector<Vector2i>			edges;
	vector<Vector3i>			trianglesEdges;
	vector<int>					adjacentTrianglesIdx;
	int							numbendings;
};