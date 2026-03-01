#pragma once
#ifndef _TRIANGLEMESH_H_
#define _TRIANGLEMESH_H_
#include "Vec3.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ======== used for merging ========== */
typedef struct {
    long long ix, iy, iz;  // 量子化キー
    int index;             // unique vertex index
    int used;              // 0=empty, 1=occupied
} HashEntry;

typedef struct {
    HashEntry* table;
    size_t size;
} VertexHash;

void vertex_hash_init(VertexHash* h, size_t size);
int vertex_hash_get_or_insert(VertexHash* h, double x, double y, double z, double eps, int* nextIndex);

int get_edge_index(int v0, int v1, int* edge,int eIndex);
void sort3(int *x, int *y, int *z);
/* =====================================
   Device Triangle Mesh
===================================== */

typedef struct DeviceTriangleMesh {
    /* vertices */
    double* mx; 
    double* my;
    double* mz;

    /* normals */
    double* nx;
    double* ny;
    double* nz;
    double *d;

    /* edge  */
    /* defined by pair of vertex indices */
    int *edge;


    /* edge vectors */
    double* e01x;
    double* e01y;
    double* e01z;

    double* e02x;
    double* e02y;
    double* e02z;

    double* e12x;
    double* e12y;
    double* e12z;

    /* barycentric */
    double* d00;
    double* d00inv;
    double* d01;
    double* d11;
    double* d11inv;
    double* d22;
    double* d22inv;
    double* denom;

    /* used for aabb */
    double *minx, *maxx;
    double *miny, *maxy;
    double *minz, *maxz;

    /* triangle vertex indices */
    int* tri_i0;
    int* tri_i1;
    int* tri_i2;

    /* triangle edge indices */
    int* tri_e0;
    int* tri_e1;
    int* tri_e2;

    int nVert; /* number of vertices */
    int nTri; /* number of triangles */
    int nShift; /*  shift of id  from vertex to edge  */

    /* bounding box as a whole triangle mesh */

    double gminx, gmaxx;
    double gminy, gmaxy;
    double gminz, gmaxz;

} DeviceTriangleMesh;
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
    double *d;

    /* edge  */
    /* defined by pair of vertex indices */
    int *edge;


    /* edge vectors */
    double* e01x;
    double* e01y;
    double* e01z;

    double* e02x;
    double* e02y;
    double* e02z;

    double* e12x;
    double* e12y;
    double* e12z;

    /* barycentric */
    double* d00;
    double* d00inv;
    double* d01;
    double* d11;
    double* d11inv;
    double* d22;
    double* d22inv;
    double* denom;

    /* used for aabb */
    double *minx, *maxx;
    double *miny, *maxy;
    double *minz, *maxz;

    /* triangle vertex indices */
    int* tri_i0;
    int* tri_i1;
    int* tri_i2;

    /* triangle edge indices */
    int* tri_e0;
    int* tri_e1;
    int* tri_e2;

    int nVert; /* number of vertices */
    int nTri; /* number of triangles */
    int nShift; /*  shift of id  from vertex to edge  */

    /* bounding box as a whole triangle mesh */

    double gminx, gmaxx;
    double gminy, gmaxy;
    double gminz, gmaxz;

    DeviceTriangleMesh d_mesh;
    DeviceTriangleMesh *d_meshPtr;

} TriangleMesh;

/* ==== device related functions ====*/
#if USE_GPU
void deviceMallocCopyTriangleMesh(TriangleMesh *mesh);
#endif

void free_TriangleMesh(TriangleMesh* mesh);
static int count_ascii_stl_triangles(FILE* fp);
int load_ascii_stl_double(const char* filename, TriangleMesh* mesh);

#endif
