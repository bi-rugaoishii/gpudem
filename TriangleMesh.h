#pragma once
#ifndef _TRIANGLEMESH_H_
#define _TRIANGLEMESH_H_
#include "Vec3.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* =====================================
   double precision triangle mesh
===================================== */
typedef struct TriangleMesh {
    /* vertices */
    double* mx; 
    double* my;
    double* mz;

    /* normals */
    double* nx;
    double* ny;
    double* nz;

    /* used for aabb */
    double *minx, *maxx;
    double *miny, *maxy;
    double *minz, *maxz;

    /* triangle indices */
    int* tri_i0;
    int* tri_i1;
    int* tri_i2;

    int nVert; /* number of vertices */
    int nTri; /* number of triangles */

} TriangleMesh;

static int count_ascii_stl_triangles(FILE* fp);
int load_ascii_stl_double(const char* filename, TriangleMesh* mesh);

#endif
