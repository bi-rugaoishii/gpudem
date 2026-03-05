#pragma once
#ifndef _BVH_H_
#define _BVH_H_

#include "TriangleMesh.h"
#include "ParticleSystem.h"
#include "TriangleContactCache.h"

typedef struct BVH{
    double* minx;
    double* miny;
    double* minz;
    double* maxx;
    double* maxy;
    double* maxz;

    int* left;
    int* right;
    int* tri;      // leafなら三角形index, 内部ノードなら-1

    int nodeCount;
} BVH;

void initializeBVH(BVH *bvh, int numTriangles);

void update_neighborlist_wall_nobvh(ParticleSystem *p,TriangleMesh *mesh,BoundingBox *box,double skinR);
TriangleContactCache dist_triangle_neighbor(ParticleSystem* ps, int i, TriangleMesh* mesh, int j,double skinR);

void update_neighborlist_wall(ParticleSystem *p,TriangleMesh *mesh, BVH *bvh,double skinR);

void free_BVH(BVH *bvh);
void computeNodeAABB(int start, int end,TriangleMesh *mesh, BVH *bvh,int k);



int buildBVH(BVH* bvh,TriangleMesh *mesh);
#endif
