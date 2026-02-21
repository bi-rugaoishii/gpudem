#ifndef _BOUNDINGBOX_H_
#define _BOUNDINGBOX_H_

#include <cuda_runtime.h>
#include <cub/cub.cuh>
#include <math.h>
#include <stdlib.h>
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

} DeviceBoundingBox;


/*
 * ==============================
 * for CPU
 * ==============================
 */


typedef struct {

    DeviceBoundingBox d_box;

    double minx,maxx,miny,maxy,minz,maxz;
    double dx,dy,dz;
    double invdx,invdy,invdz;

    double rangex,rangey,rangez;
    int sizex,sizey,sizez,N;

    int *pList;  //list of particle in cell CSR formal
    int *pNum; //number of particle in  cell
    int *pStart; //starting index of the cell in pList
    int *cellOffset; 
    
} BoundingBox;


/* =====================
 * functions
 * =====================
 */

void update_pList(ParticleSystem *p, BoundingBox *box);
void initialize_BoundingBox(ParticleSystem *p, BoundingBox *box, double minx, double maxx, double miny, double maxy, double minz, double maxz);
void free_BoundingBox(BoundingBox *box);

/* ============================
 * device related functions
 * ============================*/

__global__ void dk_build_cellCount(DeviceParticleGroup p, DeviceBoundingBox box);
__global__ void dk_build_pList(DeviceParticleGroup p, DeviceBoundingBox box);

__device__ int d_calcCellId(DeviceParticleGroup p,int i, DeviceBoundingBox box);
void d_update_pList(ParticleSystem *p, BoundingBox *box,int gridSize, int blockSize);
void copyToDeviceBox(BoundingBox *box, ParticleSystem *ps);

#endif 
