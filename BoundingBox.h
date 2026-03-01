#ifndef _BOUNDINGBOX_H_
#define _BOUNDINGBOX_H_

#include <cuda_runtime.h>
#include <cub/cub.cuh>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include "ParticleSystem.h"
#include "TriangleMesh.h"
struct ParticleSystem;
struct DeviceParticleGroup;

/*
 * ==============================
 * for GPU
 * ==============================
 */


typedef struct DeviceBoundingBox{
    double minx,maxx,miny,maxy,minz,maxz;
    double dx,dy,dz;
    double invdx,invdy,invdz;

    double rangex,rangey,rangez;
    int sizex,sizey,sizez,N;

    int *pList;  //list of particle in cell CSR formal
    int *pNum; //number of particle in  cell
    int *pStart; //starting index of the cell in pList
    int *cellOffset; 

    /* for cub */

    void *tmpExSum; /* used for exclusive sum*/
    size_t scanTmpBytes;

    /* for triangles */
    int *tList;  //list of triangles in a cell
    int *tNum; //number of particle in  cell
    int MAX_TRI;

} DeviceBoundingBox;


/*
 * ==============================
 * for CPU
 * ==============================
 */


typedef struct BoundingBox{

    DeviceBoundingBox d_box;
    DeviceBoundingBox* d_boxPtr;

    double minx,maxx,miny,maxy,minz,maxz;
    double dx,dy,dz;
    double invdx,invdy,invdz;

    double rangex,rangey,rangez;
    int sizex,sizey,sizez,N;

    int *pList;  //list of particle in cell CSR format
    int *pNum; //number of particle in  cell
    int *usedCells; //id of cells with particles
    int numUsedCells; //number of cells with particles
    int *pStart; //starting index of the cell in pList
    int *cellOffset; 

    /* for triangles */
    int *tList;  //list of triangles in a cell
    int *tNum; //number of particle in  cell
    int MAX_TRI;
    
} BoundingBox;


/* =====================
 * functions
 * =====================
 */

int compare_int(const void *a, const void *b);

void radixSortUint32(
    uint32_t** keyPtr,
    int**      indexPtr,
    int N,
    uint32_t*  workKey,
    int*       workIndex);

void swap_ps(ParticleSystem *ps, ParticleSystem *tmpps);

void update_tList(BoundingBox *box, TriangleMesh *mesh);

void update_pList(ParticleSystem *p, BoundingBox *box);
void update_pList_fast(ParticleSystem *p, BoundingBox *box);
void update_pList_withSort(ParticleSystem *p, ParticleSystem *tmpPs,BoundingBox *box);
void update_pList_withSort_fast(ParticleSystem *p, ParticleSystem *tmpPs,BoundingBox *box);

void initialize_BoundingBox(ParticleSystem *p, BoundingBox *box,TriangleMesh* mesh, double minx, double maxx, double miny, double maxy, double minz, double maxz);
void free_BoundingBox(BoundingBox *box);

/* ============================
 * device related functions
 * ============================*/

__global__ void dk_build_cellCount(DeviceParticleGroup* p, DeviceBoundingBox* box);
__global__ void dk_build_pList(DeviceParticleGroup* p, DeviceBoundingBox* box);

__device__ int d_calcCellId(DeviceParticleGroup* p,int i, DeviceBoundingBox* box);
void d_update_pList(ParticleSystem *p, BoundingBox *box,int gridSize, int blockSize);
void copyToDeviceBox(BoundingBox *box, ParticleSystem *ps);

#endif 
