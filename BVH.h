#pragma once
#ifndef _BVH_H_
#define _BVH_H_

#include "TriangleMesh.h"

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

void free_BVH(BVH *bvh);
void computeNodeAABB(int start, int end,TriangleMesh *mesh, BVH *bvh,int k);



int buildBVH(BVH* bvh,TriangleMesh *mesh);
#endif
