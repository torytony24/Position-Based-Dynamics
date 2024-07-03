#include <unordered_map>

#include "TriangularMesh.h"
#include "../Math/Util.h"

typedef unsigned long long EdgeKey;
typedef std::unordered_map<EdgeKey, int> EdgeMapInt;
typedef std::unordered_map<EdgeKey, Vector2i> EdgeMapVec2i;

static int findEdge(EdgeMapInt& edgeMap, const Vector2i& e) {
	Vector2i tmpE = e;
	tmpE[0] = MIN(e[0], e[1]);
	tmpE[1] = MAX(e[0], e[1]);
	EdgeKey eKey = tmpE[0];
	eKey = eKey << 32;
	eKey |= tmpE[1];

	EdgeMapInt::iterator itr = edgeMap.find(eKey);

	if (itr != edgeMap.end()) // find 
		return itr->second;
	else // do not find 
		return -1;
}

static void addEdge(EdgeMapVec2i& edgeMap, std::vector<Vector2i>& edgeList, const Vector2i& e) {
	Vector2i tmpE = e;
	tmpE[0] = MIN(e[0], e[1]);
	tmpE[1] = MAX(e[0], e[1]);
	EdgeKey eKey = tmpE[0];
	eKey = eKey << 32;
	eKey |= tmpE[1];
	if (edgeMap.count(eKey) == 0) { // edgeMap doesn't have this "eKey", then add Edge
		edgeMap.insert(std::pair<EdgeKey, Vector2i>(eKey, e));
		edgeList.push_back(e);
	}
}


TriangularMesh::TriangularMesh() : Object() {

}

TriangularMesh::TriangularMesh(const char* filePath) : Object(filePath) {
	initGeometryInformation();
}

TriangularMesh::~TriangularMesh() {

}

void TriangularMesh::initGeometryInformation() {
	initEdges();
	initAdjacentTrianglesElements();
	//printMeshInfo();
}

void TriangularMesh::printMeshInfo() {
	cout << "===================== Mesh Information ====================" << endl;
	cout << "Number of Vertices       : " << vertices.size() << endl;
	cout << "Number of Edges          : " << edges.size() << endl;
	cout << "Number of Triangles      : " << triangles.size() << endl;
	cout << "Number of Triangle Pairs : " << adjacentTrianglesIdx.size() / 6 << endl << endl;
	cout << "============================================================" << endl;
}

void TriangularMesh::initEdges() {
	edges.clear();

	EdgeMapVec2i edgeMap;
	for (int i = 0; i < (int)triangles.size(); ++i) {
		Vector2i edge1, edge2, edge3;
		edge1[0] = triangles[i].vindices[0];
		edge1[1] = triangles[i].vindices[1];

		edge2[0] = triangles[i].vindices[1];
		edge2[1] = triangles[i].vindices[2];

		edge3[0] = triangles[i].vindices[2];
		edge3[1] = triangles[i].vindices[0];

		addEdge(edgeMap, edges, edge1);
		addEdge(edgeMap, edges, edge2);
		addEdge(edgeMap, edges, edge3);
	}
}


void TriangularMesh::initAdjacentTrianglesElements() {
	EdgeMapInt edgeMapInt;
	int numE = (int)edges.size();
	for (int i = 0; i < numE; ++i) {
		Vector2i tmpE = edges[i];
		tmpE[0] = MIN(edges[i][0], edges[i][1]);
		tmpE[1] = MAX(edges[i][0], edges[i][1]);
		EdgeKey eKey = tmpE[0];
		eKey = eKey << 32;
		eKey |= tmpE[1];
		edgeMapInt.insert(std::make_pair(eKey, i));
	}

	vector<vector<int>> edgeSharedTriangle;
	edgeSharedTriangle.resize(edges.size());

	trianglesEdges.clear();
	trianglesEdges.resize(numtriangles);

	for (int i = 0; i < numtriangles; i++) {
		Vector2i e1(triangles[i].vindices[0], triangles[i].vindices[1]);
		Vector2i e2(triangles[i].vindices[1], triangles[i].vindices[2]);
		Vector2i e3(triangles[i].vindices[2], triangles[i].vindices[0]);
		int e1Idx = findEdge(edgeMapInt, e1);
		int e2Idx = findEdge(edgeMapInt, e2);
		int e3Idx = findEdge(edgeMapInt, e3);
		edgeSharedTriangle[e1Idx].push_back(i);
		edgeSharedTriangle[e2Idx].push_back(i);
		edgeSharedTriangle[e3Idx].push_back(i);
		trianglesEdges[i] = Vector3i(e1Idx, e2Idx, e3Idx);
	}

	for (int i = 0; i < numE; i++) {
		if (edgeSharedTriangle[i].size() > 1) {
			int idx1 = edgeSharedTriangle[i][0];
			int idx2 = edgeSharedTriangle[i][1];

			int v[6];
			for (int i = 0; i < 3; i++) {
				bool bSharedVertex = false;
				for (int j = 0; j < 3; j++)
					if (triangles[idx1].vindices[i] == triangles[idx2].vindices[j])
						bSharedVertex = true;

				if (!bSharedVertex) {
					// v[0], v[1] : unshared vertex
					// v[2], v[3] : shared vertex
					// v[4], v[5] : triagle idx
					v[0] = triangles[idx1].vindices[i];
					v[2] = triangles[idx1].vindices[(i + 2) % 3];
					v[3] = triangles[idx1].vindices[(i + 1) % 3];
					v[4] = idx1;
					v[5] = idx2;
					for (int p = 0; p < 3; ++p) {
						if (triangles[idx2].vindices[p] != v[2] && 
							triangles[idx2].vindices[p] != v[3]) {
							v[1] = triangles[idx2].vindices[p];
							break;
						}
					}
					break;
				}
			}
			for (int j = 0; j < 6; j++)
				adjacentTrianglesIdx.push_back(v[j]);
		}
	}
}
