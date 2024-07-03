#include "Object.h"
#include <gl/glut.h>

Object::Object() {

}

// Construct with *.obj file
Object::Object(const char* filePath) {
	readOBJFile(filePath);
	updateFacetNormals();
	updateVertexNormals(0.f);
}

Object::~Object() {
	vertices.clear();
	normals.clear();
	texcoords.clear();
	facetnorms.clear();
	triangles.clear();
	materials.clear();
	groups.clear();
}

OBJGroup* Object::findGroup(char* name) {
	for (int i = 0; i < groups.size(); i++)
		if (!strcmp(name, groups[i].name))
			return &groups[i];

	return NULL;
}

int Object::findMaterial(char* name) {
	int i;

	/* XXX doing a linear search on a string key'd list is pretty lame,
	but it works and is fast enough for now. */
	for (i = 0; i < nummaterials; i++) {
		if (!strcmp(materials[i].name, name))
			return i;
	}

	/* didn't find the name, so print a warning and return the default material (0). */
	printf("glmFindMaterial():  can't find material \"%s\".\n", name);
	return 0;
}

void Object::readOBJFile(const char* filePath) {
	FILE* file;
	errno_t err;

	/* open the file */
	err = fopen_s(&file, filePath, "r");
	if (err) {
		fprintf(stderr, "ReadOBJ failed: can't open data file \"%s\".\n", filePath);
		exit(1);
	}

	/* allocate a new model */
	pathname = _strdup(filePath);
	numvertices = 0;
	numnormals = 0;
	numtexcoords = 0;
	numfacetnorms = 0;
	numtriangles = 0;
	nummaterials = 0;
	numgroups = 0;

	/* make a first pass through the file to get a count of the number
	of vertices, normals, texcoords & triangles */
	firstPass(file);

	/* allocate memory */
	vertices.resize(numvertices);
	triangles.resize(numtriangles);
	if (numnormals) normals.resize(numnormals);
	if (numtexcoords) texcoords.resize(numtexcoords);

	/* rewind to beginning of file and read in the data this pass */
	rewind(file);

	secondPass(file);

	/* close the file */
	fclose(file);
}

void Object::firstPass(FILE* file) {
	int numvertices;        /* number of vertices in model */
	int numnormals;         /* number of normals in model */
	int numtexcoords;       /* number of texcoords in model */
	int numtriangles;       /* number of triangles in model */
	OBJGroup* group;           /* current group */
	int v, n, t;
	char buf[128];
	char buf2[128];
	char buf3[128];

	/* make a default group */
	groups.push_back(OBJGroup("default"));
	numgroups++;
	group = &(groups.front());

	numvertices = numnormals = numtexcoords = numtriangles = 0;
	while (fscanf_s(file, "%s", buf, sizeof(buf)) != EOF) {
		switch (buf[0]) {

		case '#':               /* comment */
			/* eat up rest of line */
			fgets(buf, sizeof(buf), file);
			break;
		case 'v':               /* v, vn, vt */
			switch (buf[1]) {
			case '\0':          /* vertex */
				/* eat up rest of line */
				fgets(buf, sizeof(buf), file);
				numvertices++;
				break;
			case 'n':           /* normal */
				/* eat up rest of line */
				fgets(buf, sizeof(buf), file);
				numnormals++;
				break;
			case 't':           /* texcoord */
				/* eat up rest of line */
				fgets(buf, sizeof(buf), file);
				numtexcoords++;
				break;
			default:
				printf("glmFirstPass(): Unknown token \"%s\".\n", buf);
				exit(1);
				break;
			}
			break;
		case 'm':
			fgets(buf, sizeof(buf), file);
			sscanf_s(buf, "%s", buf, sizeof(buf));
			mtllibname = _strdup(buf);
			readMTLFile(buf);
			break;
		case 'u':
			/* eat up rest of line */
			fgets(buf, sizeof(buf), file);
			break;
		case 'g':               /* group */
			/* eat up rest of line */
			fgets(buf, sizeof(buf), file);
			buf[strlen(buf) - 1] = '\0';  /* nuke '\n' */
			if (!(group = findGroup(buf))) {
				groups.push_back(OBJGroup(buf));
				group = &(groups.back());
			}
			break;
		case 'f':               /* face */
			v = n = t = 0;
			fscanf_s(file, "%s", buf, sizeof(buf));
			/* can be one of %d, %d//%d, %d/%d, %d/%d/%d %d//%d */
			if (strstr(buf, "//")) {
				/* v//n */
				sscanf_s(buf, "%d//%d", &v, &n);
				fscanf_s(file, "%d//%d", &v, &n);
				fscanf_s(file, "%d//%d", &v, &n);
				numtriangles++;
				group->numtriangles++;
				while (fscanf_s(file, "%d//%d", &v, &n) > 0) {
					numtriangles++;
					group->numtriangles++;
				}
			}
			else if (sscanf_s(buf, "%d/%d/%d", &v, &t, &n) == 3) {
				/* v/t/n */
				fscanf_s(file, "%d/%d/%d", &v, &t, &n);
				fscanf_s(file, "%d/%d/%d", &v, &t, &n);
				numtriangles++;
				group->numtriangles++;
				while (fscanf_s(file, "%d/%d/%d", &v, &t, &n) > 0) {
					numtriangles++;
					group->numtriangles++;
				}
			}
			else if (sscanf_s(buf, "%d/%d", &v, &t) == 2) {
				/* v/t */
				fscanf_s(file, "%d/%d", &v, &t);
				fscanf_s(file, "%d/%d", &v, &t);
				numtriangles++;
				group->numtriangles++;
				while (fscanf_s(file, "%d/%d", &v, &t) > 0) {
					numtriangles++;
					group->numtriangles++;
				}
			}
			else {
				/* v */
				fscanf_s(file, "%d", &v);
				fscanf_s(file, "%d", &v);
				numtriangles++;
				group->numtriangles++;
				while (fscanf_s(file, "%d", &v) > 0) {
					numtriangles++;
					group->numtriangles++;
				}
			}
			break;

		default:
			/* eat up rest of line */
			fgets(buf, sizeof(buf), file);
			break;
		}
	}

	/* set the stats in the model structure */
	this->numvertices = numvertices;
	this->numnormals = numnormals;
	this->numtexcoords = numtexcoords;
	this->numtriangles = numtriangles;

	/* allocate memory for the triangles in each group */
	for (int i = 0; i < groups.size(); i++) {
		groups[i].triangles.resize(groups[i].numtriangles);
		groups[i].numtriangles = 0;
	}
}

void Object::secondPass(FILE* file) {
	int numvertices;        /* number of vertices in model */
	int numnormals;         /* number of normals in model */
	int numtexcoords;       /* number of texcoords in model */
	int numtriangles;       /* number of triangles in model */
	OBJGroup* group;        /* current group pointer */
	int material;           /* current material */
	int v, n, t;
	char buf[128];

	/* set the pointer shortcuts */
	group = &(groups.back());

	/* on the second pass through the file, read all the data into the
	allocated arrays */
	numvertices = numnormals = numtexcoords = 0;
	numtriangles = 0;
	material = 0;
	while (fscanf_s(file, "%s", buf, sizeof(buf)) != EOF) {
		switch (buf[0]) {
		case '#':               /* comment */
			/* eat up rest of line */
			fgets(buf, sizeof(buf), file);
			break;
		case 'v':               /* v, vn, vt */
			switch (buf[1]) {
			case '\0':          /* vertex */
				fscanf_s(file, "%lf %lf %lf",
					&vertices[numvertices][0],
					&vertices[numvertices][1],
					&vertices[numvertices][2]);
				numvertices++;
				break;
			case 'n':           /* normal */
				fscanf_s(file, "%lf %lf %lf",
					&normals[numnormals][0],
					&normals[numnormals][1],
					&normals[numnormals][2]);
				numnormals++;
				break;
			case 't':           /* texcoord */
				fscanf_s(file, "%lf %lf",
					&texcoords[numtexcoords][0],
					&texcoords[numtexcoords][1]);
				numtexcoords++;
				break;
			}
			break;
		case 'u':
			fgets(buf, sizeof(buf), file);
			sscanf_s(buf, "%s", buf, sizeof(buf));
			group->material = material = findMaterial(buf);
			break;
		case 'g':               /* group */
			/* eat up rest of line */
			fgets(buf, sizeof(buf), file);
			buf[strlen(buf) - 1] = '\0';  /* nuke '\n' */
			group = findGroup(buf);
			group->material = material;
			break;
		case 'f':               /* face */
			v = n = t = 0;
			fscanf_s(file, "%s", buf, sizeof(buf));
			/* can be one of %d, %d//%d, %d/%d, %d/%d/%d %d//%d */
			if (strstr(buf, "//")) {
				/* v//n */
				sscanf_s(buf, "%d//%d", &v, &n);
				triangles[numtriangles].vindices[0] = v < 0 ? v + numvertices - 1 : v - 1;
				triangles[numtriangles].nindices[0] = n < 0 ? n + numnormals - 1 : n - 1;
				fscanf_s(file, "%d//%d", &v, &n);
				triangles[numtriangles].vindices[1] = v < 0 ? v + numvertices : v - 1;
				triangles[numtriangles].nindices[1] = n < 0 ? n + numnormals : n - 1;
				fscanf_s(file, "%d//%d", &v, &n);
				triangles[numtriangles].vindices[2] = v < 0 ? v + numvertices - 1 : v - 1;
				triangles[numtriangles].nindices[2] = n < 0 ? n + numnormals - 1 : n - 1;
				group->triangles[group->numtriangles++] = numtriangles;
				numtriangles++;
				while (fscanf_s(file, "%d//%d", &v, &n) > 0) {
					triangles[numtriangles].vindices[0] = triangles[numtriangles - 1].vindices[0];
					triangles[numtriangles].nindices[0] = triangles[numtriangles - 1].nindices[0];
					triangles[numtriangles].vindices[1] = triangles[numtriangles - 1].vindices[2];
					triangles[numtriangles].nindices[1] = triangles[numtriangles - 1].nindices[2];
					triangles[numtriangles].vindices[2] = v < 0 ? v + numvertices - 1 : v - 1;
					triangles[numtriangles].nindices[2] = n < 0 ? n + numnormals - 1 : n - 1;
					group->triangles[group->numtriangles++] = numtriangles;
					numtriangles++;
				}
			}
			else if (sscanf_s(buf, "%d/%d/%d", &v, &t, &n) == 3) {
				/* v/t/n */
				triangles[numtriangles].vindices[0] = v < 0 ? v + numvertices - 1 : v - 1;
				triangles[numtriangles].tindices[0] = t < 0 ? t + numtexcoords - 1 : t - 1;
				triangles[numtriangles].nindices[0] = n < 0 ? n + numnormals - 1 : n - 1;
				fscanf_s(file, "%d/%d/%d", &v, &t, &n);
				triangles[numtriangles].vindices[1] = v < 0 ? v + numvertices - 1 : v - 1;
				triangles[numtriangles].tindices[1] = t < 0 ? t + numtexcoords - 1 : t - 1;
				triangles[numtriangles].nindices[1] = n < 0 ? n + numnormals - 1 : n - 1;
				fscanf_s(file, "%d/%d/%d", &v, &t, &n);
				triangles[numtriangles].vindices[2] = v < 0 ? v + numvertices - 1 : v - 1;
				triangles[numtriangles].tindices[2] = t < 0 ? t + numtexcoords - 1 : t - 1;
				triangles[numtriangles].nindices[2] = n < 0 ? n + numnormals - 1 : n - 1;
				group->triangles[group->numtriangles++] = numtriangles;
				numtriangles++;
				while (fscanf_s(file, "%d/%d/%d", &v, &t, &n) > 0) {
					triangles[numtriangles].vindices[0] = triangles[numtriangles - 1].vindices[0];
					triangles[numtriangles].tindices[0] = triangles[numtriangles - 1].tindices[0];
					triangles[numtriangles].nindices[0] = triangles[numtriangles - 1].nindices[0];
					triangles[numtriangles].vindices[1] = triangles[numtriangles - 1].vindices[2];
					triangles[numtriangles].tindices[1] = triangles[numtriangles - 1].tindices[2];
					triangles[numtriangles].nindices[1] = triangles[numtriangles - 1].nindices[2];
					triangles[numtriangles].vindices[2] = v < 0 ? v + numvertices - 1 : v - 1;
					triangles[numtriangles].tindices[2] = t < 0 ? t + numtexcoords - 1 : t - 1;
					triangles[numtriangles].nindices[2] = n < 0 ? n + numnormals - 1 : n - 1;
					group->triangles[group->numtriangles++] = numtriangles;
					numtriangles++;
				}
			}
			else if (sscanf_s(buf, "%d/%d", &v, &t) == 2) {
				/* v/t */
				triangles[numtriangles].vindices[0] = v < 0 ? v - 1 + numvertices : v - 1;
				triangles[numtriangles].tindices[0] = t < 0 ? t - 1 + numtexcoords : t - 1;
				fscanf_s(file, "%d/%d", &v, &t);
				triangles[numtriangles].vindices[1] = v < 0 ? v - 1 + numvertices : v - 1;
				triangles[numtriangles].tindices[1] = t < 0 ? t - 1 + numtexcoords : t - 1;
				fscanf_s(file, "%d/%d", &v, &t);
				triangles[numtriangles].vindices[2] = v < 0 ? v - 1 + numvertices : v - 1;
				triangles[numtriangles].tindices[2] = t < 0 ? t - 1 + numtexcoords : t - 1;
				group->triangles[group->numtriangles++] = numtriangles;
				numtriangles++;
				while (fscanf_s(file, "%d/%d", &v, &t) > 0) {
					triangles[numtriangles].vindices[0] = triangles[numtriangles - 1].vindices[0];
					triangles[numtriangles].tindices[0] = triangles[numtriangles - 1].tindices[0];
					triangles[numtriangles].vindices[1] = triangles[numtriangles - 1].vindices[2];
					triangles[numtriangles].tindices[1] = triangles[numtriangles - 1].tindices[2];
					triangles[numtriangles].vindices[2] = v < 0 ? v + numvertices - 1 : v - 1;
					triangles[numtriangles].tindices[2] = t < 0 ? t + numtexcoords - 1 : t - 1;
					group->triangles[group->numtriangles++] = numtriangles;
					numtriangles++;
				}
			}
			else {
				/* v */
				sscanf_s(buf, "%d", &v);
				triangles[numtriangles].vindices[0] = v < 0 ? v + numvertices - 1 : v - 1;
				fscanf_s(file, "%d", &v);
				triangles[numtriangles].vindices[1] = v < 0 ? v + numvertices - 1 : v - 1;
				fscanf_s(file, "%d", &v);
				triangles[numtriangles].vindices[2] = v < 0 ? v + numvertices - 1 : v - 1;
				group->triangles[group->numtriangles++] = numtriangles;
				numtriangles++;
				while (fscanf_s(file, "%d", &v) > 0) {
					triangles[numtriangles].vindices[0] = triangles[numtriangles - 1].vindices[0];
					triangles[numtriangles].vindices[1] = triangles[numtriangles - 1].vindices[2];
					triangles[numtriangles].vindices[2] = v < 0 ? v + numvertices - 1 : v - 1;
					group->triangles[group->numtriangles++] = numtriangles;
					numtriangles++;
				}
			}
			break;

		default:
			/* eat up rest of line */
			fgets(buf, sizeof(buf), file);
			break;
		}
	}
}

void Object::writeOBJFile(char* filePath) {
	//GLuint i;
	//FILE* file;
	//OBJGroup* group;

	/* do a bit of warning */
	//if (mode & GLM_FLAT && !model->facetnorms) {
	//    printf("glmWriteOBJ() warning: flat normal output requested "
	//        "with no facet normals defined.\n");
	//    mode &= ~GLM_FLAT;
	//}
	//if (mode & GLM_SMOOTH && !model->normals) {
	//    printf("glmWriteOBJ() warning: smooth normal output requested "
	//        "with no normals defined.\n");
	//    mode &= ~GLM_SMOOTH;
	//}
	//if (mode & GLM_TEXTURE && !model->texcoords) {
	//    printf("glmWriteOBJ() warning: texture coordinate output requested "
	//        "with no texture coordinates defined.\n");
	//    mode &= ~GLM_TEXTURE;
	//}
	//if (mode & GLM_FLAT && mode & GLM_SMOOTH) {
	//    printf("glmWriteOBJ() warning: flat normal output requested "
	//        "and smooth normal output requested (using smooth).\n");
	//    mode &= ~GLM_FLAT;
	//}
	//if (mode & GLM_COLOR && !model->materials) {
	//    printf("glmWriteOBJ() warning: color output requested "
	//        "with no colors (materials) defined.\n");
	//    mode &= ~GLM_COLOR;
	//}
	//if (mode & GLM_MATERIAL && !model->materials) {
	//    printf("glmWriteOBJ() warning: material output requested "
	//        "with no materials defined.\n");
	//    mode &= ~GLM_MATERIAL;
	//}
	//if (mode & GLM_COLOR && mode & GLM_MATERIAL) {
	//    printf("glmWriteOBJ() warning: color and material output requested "
	//        "outputting only materials.\n");
	//    mode &= ~GLM_COLOR;
	//}
	//
	//
	///* open the file */
	//file = fopen(filename, "w");
	//if (!file) {
	//    fprintf(stderr, "glmWriteOBJ() failed: can't open file \"%s\" to write.\n",
	//        filename);
	//    exit(1);
	//}
	//
	///* spit out a header */
	//fprintf(file, "#  \n");
	//fprintf(file, "#  Wavefront OBJ generated by GLM library\n");
	//fprintf(file, "#  \n");
	//fprintf(file, "#  GLM library\n");
	//fprintf(file, "#  Nate Robins\n");
	//fprintf(file, "#  ndr@pobox.com\n");
	//fprintf(file, "#  http://www.pobox.com/~ndr\n");
	//fprintf(file, "#  \n");
	//
	//if (mode & GLM_MATERIAL && model->mtllibname) {
	//    fprintf(file, "\nmtllib %s\n\n", model->mtllibname);
	//    glmWriteMTL(model, filename, model->mtllibname);
	//}
	//
	///* spit out the vertices */
	//fprintf(file, "\n");
	//fprintf(file, "# %d vertices\n", model->numvertices);
	//for (i = 1; i <= model->numvertices; i++) {
	//    fprintf(file, "v %f %f %f\n", 
	//        model->vertices[3 * i + 0],
	//        model->vertices[3 * i + 1],
	//        model->vertices[3 * i + 2]);
	//}
	//
	///* spit out the smooth/flat normals */
	//if (mode & GLM_SMOOTH) {
	//    fprintf(file, "\n");
	//    fprintf(file, "# %d normals\n", model->numnormals);
	//    for (i = 1; i <= model->numnormals; i++) {
	//        fprintf(file, "vn %f %f %f\n", 
	//            model->normals[3 * i + 0],
	//            model->normals[3 * i + 1],
	//            model->normals[3 * i + 2]);
	//    }
	//} else if (mode & GLM_FLAT) {
	//    fprintf(file, "\n");
	//    fprintf(file, "# %d normals\n", model->numfacetnorms);
	//    for (i = 1; i <= model->numnormals; i++) {
	//        fprintf(file, "vn %f %f %f\n", 
	//            model->facetnorms[3 * i + 0],
	//            model->facetnorms[3 * i + 1],
	//            model->facetnorms[3 * i + 2]);
	//    }
	//}
	//
	///* spit out the texture coordinates */
	//if (mode & GLM_TEXTURE) {
	//    fprintf(file, "\n");
	//    fprintf(file, "# %d texcoords\n", model->numtexcoords);
	//    for (i = 1; i <= model->numtexcoords; i++) {
	//        fprintf(file, "vt %f %f\n", 
	//            model->texcoords[2 * i + 0],
	//            model->texcoords[2 * i + 1]);
	//    }
	//}
	//
	//fprintf(file, "\n");
	//fprintf(file, "# %d groups\n", model->numgroups);
	//fprintf(file, "# %d faces (triangles)\n", model->numtriangles);
	//fprintf(file, "\n");
	//
	//group = model->groups;
	//while(group) {
	//    fprintf(file, "g %s\n", group->name);
	//    if (mode & GLM_MATERIAL)
	//        fprintf(file, "usemtl %s\n", model->materials[group->material].name);
	//    for (i = 0; i < group->numtriangles; i++) {
	//        if (mode & GLM_SMOOTH && mode & GLM_TEXTURE) {
	//            fprintf(file, "f %d/%d/%d %d/%d/%d %d/%d/%d\n",
	//                T(group->triangles[i]).vindices[0], 
	//                T(group->triangles[i]).tindices[0],
	//                T(group->triangles[i]).nindices[0], 
	//                T(group->triangles[i]).vindices[1],
	//                T(group->triangles[i]).tindices[1],
	//                T(group->triangles[i]).nindices[1],
	//                T(group->triangles[i]).vindices[2],
	//                T(group->triangles[i]).tindices[2],
	//                T(group->triangles[i]).nindices[2]);
	//        } else if (mode & GLM_FLAT && mode & GLM_TEXTURE) {
	//            fprintf(file, "f %d/%d %d/%d %d/%d\n",
	//                T(group->triangles[i]).vindices[0],
	//                T(group->triangles[i]).findex,
	//                T(group->triangles[i]).vindices[1],
	//                T(group->triangles[i]).findex,
	//                T(group->triangles[i]).vindices[2],
	//                T(group->triangles[i]).findex);
	//        } else if (mode & GLM_TEXTURE) {
	//            fprintf(file, "f %d/%d %d/%d %d/%d\n",
	//                T(group->triangles[i]).vindices[0],
	//                T(group->triangles[i]).tindices[0],
	//                T(group->triangles[i]).vindices[1],
	//                T(group->triangles[i]).tindices[1],
	//                T(group->triangles[i]).vindices[2],
	//                T(group->triangles[i]).tindices[2]);
	//        } else if (mode & GLM_SMOOTH) {
	//            fprintf(file, "f %d//%d %d//%d %d//%d\n",
	//                T(group->triangles[i]).vindices[0],
	//                T(group->triangles[i]).nindices[0],
	//                T(group->triangles[i]).vindices[1],
	//                T(group->triangles[i]).nindices[1],
	//                T(group->triangles[i]).vindices[2], 
	//                T(group->triangles[i]).nindices[2]);
	//        } else if (mode & GLM_FLAT) {
	//            fprintf(file, "f %d//%d %d//%d %d//%d\n",
	//                T(group->triangles[i]).vindices[0], 
	//                T(group->triangles[i]).findex,
	//                T(group->triangles[i]).vindices[1],
	//                T(group->triangles[i]).findex,
	//                T(group->triangles[i]).vindices[2],
	//                T(group->triangles[i]).findex);
	//        } else {
	//            fprintf(file, "f %d %d %d\n",
	//                T(group->triangles[i]).vindices[0],
	//                T(group->triangles[i]).vindices[1],
	//                T(group->triangles[i]).vindices[2]);
	//        }
	//    }
	//    fprintf(file, "\n");
	//    group = group->next;
	//}
	//
	//fclose(file);
}

void Object::readMTLFile(char* name) {
	FILE* file;
	errno_t err;
	char* dir;
	char* filename;
	char buf[128];
	GLuint nummaterials, i;

	char* temp = _strdup(pathname);
	char* s;
	s = strrchr(temp, '/');
	if (s)
		s[1] = '\0';
	else
		temp[0] = '\0';

	dir = temp;
	//dir = glmDirName(pathname);
	size_t file_len = strlen(dir) + strlen(name) + 1;
	filename = new char[file_len];
	strcpy_s(filename, file_len, dir);
	strcat_s(filename, file_len, name);
	free(dir);

	err = fopen_s(&file, filename, "r");
	if (err) {
		fprintf(stderr, "glmReadMTL() failed: can't open material file \"%s\".\n",
			filename);
		exit(1);
	}
	free(filename);

	/* count the number of materials in the file */
	nummaterials = 1;
	while (fscanf_s(file, "%s", buf, sizeof(buf)) != EOF) {
		switch (buf[0]) {
		case '#':               /* comment */
			/* eat up rest of line */
			fgets(buf, sizeof(buf), file);
			break;
		case 'n':               /* newmtl */
			fgets(buf, sizeof(buf), file);
			nummaterials++;
			sscanf_s(buf, "%s", buf, sizeof(buf));
			break;
		default:
			/* eat up rest of line */
			fgets(buf, sizeof(buf), file);
			break;
		}
	}

	rewind(file);

	materials.resize(nummaterials);
	this->nummaterials = nummaterials;

	/* set the default material */
	for (i = 0; i < nummaterials; i++) {
		materials[i].name = NULL;
		materials[i].shininess = 65.0f;
		materials[i].diffuse[0] = 0.8f;
		materials[i].diffuse[1] = 0.8f;
		materials[i].diffuse[2] = 0.8f;
		materials[i].diffuse[3] = 1.0f;
		materials[i].ambient[0] = 0.2f;
		materials[i].ambient[1] = 0.2f;
		materials[i].ambient[2] = 0.2f;
		materials[i].ambient[3] = 1.0f;
		materials[i].specular[0] = 0.0f;
		materials[i].specular[1] = 0.0f;
		materials[i].specular[2] = 0.0f;
		materials[i].specular[3] = 1.0f;
	}
	materials[0].name = _strdup("default");

	/* now, read in the data */
	nummaterials = 0;
	while (fscanf_s(file, "%s", buf, sizeof(buf)) != EOF) {
		switch (buf[0]) {
		case '#':               /* comment */
			/* eat up rest of line */
			fgets(buf, sizeof(buf), file);
			break;
		case 'n':               /* newmtl */
			fgets(buf, sizeof(buf), file);
			sscanf_s(buf, "%s", buf, sizeof(buf));
			nummaterials++;
			materials[nummaterials].name = _strdup(buf);
			break;
		case 'N':
			fscanf_s(file, "%f", &materials[nummaterials].shininess);
			/* wavefront shininess is from [0, 1000], so scale for OpenGL */
			materials[nummaterials].shininess /= 1000.0;
			materials[nummaterials].shininess *= 128.0;
			break;
		case 'K':
			switch (buf[1]) {
			case 'd':
				fscanf_s(file, "%f %f %f",
					&materials[nummaterials].diffuse[0],
					&materials[nummaterials].diffuse[1],
					&materials[nummaterials].diffuse[2]);
				break;
			case 's':
				fscanf_s(file, "%f %f %f",
					&materials[nummaterials].specular[0],
					&materials[nummaterials].specular[1],
					&materials[nummaterials].specular[2]);
				break;
			case 'a':
				fscanf_s(file, "%f %f %f",
					&materials[nummaterials].ambient[0],
					&materials[nummaterials].ambient[1],
					&materials[nummaterials].ambient[2]);
				break;
			default:
				/* eat up rest of line */
				fgets(buf, sizeof(buf), file);
				break;
			}
			break;
		default:
			/* eat up rest of line */
			fgets(buf, sizeof(buf), file);
			break;
		}
	}
}

void Object::writeMTLFile(char* filePath) {
	//FILE* file;
 //   char* dir;
 //   char* filename;
 //   GLMmaterial* material;
 //   GLuint i;
 //   
 //   dir = glmDirName(modelpath);
 //   filename = (char*)malloc(sizeof(char) * (strlen(dir)+strlen(mtllibname)));
 //   strcpy(filename, dir);
 //   strcat(filename, mtllibname);
 //   free(dir);
 //   
 //   /* open the file */
 //   file = fopen(filename, "w");
 //   if (!file) {
 //       fprintf(stderr, "glmWriteMTL() failed: can't open file \"%s\".\n",
 //           filename);
 //       exit(1);
 //   }
 //   free(filename);
 //   
 //   /* spit out a header */
 //   fprintf(file, "#  \n");
 //   fprintf(file, "#  Wavefront MTL generated by GLM library\n");
 //   fprintf(file, "#  \n");
 //   fprintf(file, "#  GLM library\n");
 //   fprintf(file, "#  Nate Robins\n");
 //   fprintf(file, "#  ndr@pobox.com\n");
 //   fprintf(file, "#  http://www.pobox.com/~ndr\n");
 //   fprintf(file, "#  \n\n");
 //   
 //   for (i = 0; i < model->nummaterials; i++) {
 //       material = &model->materials[i];
 //       fprintf(file, "newmtl %s\n", material->name);
 //       fprintf(file, "Ka %f %f %f\n", 
 //           material->ambient[0], material->ambient[1], material->ambient[2]);
 //       fprintf(file, "Kd %f %f %f\n", 
 //           material->diffuse[0], material->diffuse[1], material->diffuse[2]);
 //       fprintf(file, "Ks %f %f %f\n", 
 //           material->specular[0],material->specular[1],material->specular[2]);
 //       fprintf(file, "Ns %f\n", material->shininess / 128.0 * 1000.0);
 //       fprintf(file, "\n");
 //   }
}

void Object::updateFacetNormals() {
	GLuint  i;
	Vector3d u;
	Vector3d v;

	if (!(facetnorms.empty()))
		facetnorms.clear();

	numfacetnorms = numtriangles;
	facetnorms.resize(numfacetnorms);

	for (i = 0; i < numtriangles; i++) {
		triangles[i].findex = i;

		Vector3d a = vertices[triangles[i].vindices[0]];
		Vector3d b = vertices[triangles[i].vindices[1]];
		Vector3d c = vertices[triangles[i].vindices[2]];

		u = b - a;
		v = c - a;

		facetnorms[i] = cross(u, v);
		normalize(facetnorms[i]);
	}
}

void Object::updateVertexNormals(float angle) {
	if (normals.size() > 0)
		normals.clear();

	vector<vector<int>> facetIndices;

	facetIndices.resize(numvertices);

	for (int i = 0; i < numtriangles; i++) {
		int facetIndex = triangles[i].findex;
		facetIndices[triangles[i].vindices[0]].push_back(facetIndex);
		facetIndices[triangles[i].vindices[1]].push_back(facetIndex);
		facetIndices[triangles[i].vindices[2]].push_back(facetIndex);
	}

	normals.resize(numvertices);

	for (int i = 0; i < numvertices; i++) {
		Vector3d normal;
		int numOfNorm = facetIndices[i].size();
		for (int j = 0; j < numOfNorm; j++)
			normal += facetnorms[facetIndices[i][j]];
		normals[i] = normal / (double)numOfNorm;
	}
}