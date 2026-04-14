#pragma once

#include <cuda_runtime.h>
#include <cub/cub.cuh>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include "ParticleSystem.h"
#include "hardCodedParameters.h"
#include "cJSON.h"
struct ParticleSystem;
struct HostMemory;
struct DeviceMemory;
struct Common;

template <typename Memory>// typename used later to distinguish gpu and cpu
struct ParticleSys;

template <>
struct ParticleSys<HostMemory>;

template <>
struct ParticleSys<DeviceMemory>;

struct DeviceParticleGroup;
struct TriangleMesh;

/*
 * ==============================
 * for GPU
 * ==============================
 */


typedef struct DeviceBoundingBox{
    double minx,maxx,miny,maxy,minz,maxz;
    double dx,dy,dz;
    double invdx,invdy,invdz;

    double skinR;
    double refreshThresh;
    double refreshThreshSq;

    double rangex,rangey,rangez;
    int sizex,sizey,sizez,N;

    int *pList;  //list of particle in cell CSR formal
    int *pNum; //number of particle in  cell
    int *usedCells; //id of cells with particles
    int numUsedCells; //number of cells with particles
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

    double skinR;
    double refreshThresh;
    double refreshThreshSq;

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

/* ============ Morton key related ================ */

__device__ __host__ __forceinline__ uint32_t expandBits(uint32_t v){
    v &= 0x000003ff;                 // 10bit
    v = (v | (v << 16)) & 0x030000FF;
    v = (v | (v <<  8)) & 0x0300F00F;
    v = (v | (v <<  4)) & 0x030C30C3;
    v = (v | (v <<  2)) & 0x09249249;
    return v;
}
__device__ __host__ __forceinline__  uint32_t morton3D(uint32_t ix, uint32_t iy,uint32_t iz){
    return (expandBits(ix) << 2)| (expandBits(iy) << 1)| (expandBits(iz));
}

void radixSortUint32(
    uint32_t** keyPtr,
    int**      indexPtr,
    int N,
    uint32_t*  workKey,
    int*       workIndex);

void swap_ps(ParticleSys<HostMemory> *ps, ParticleSys<HostMemory> *tmpps);
void swap_device_ps_pointer(DeviceParticleGroup **p, DeviceParticleGroup **tmp);

void update_tList(BoundingBox *box, TriangleMesh *mesh);

void update_neighborlist(ParticleSys<HostMemory> *p,ParticleSys<HostMemory> *tmpP, BoundingBox *box);

void update_pList(ParticleSystem *p, BoundingBox *box);
void update_pList_fast(ParticleSystem *p, BoundingBox *box);
void update_pList_withSort(ParticleSys<HostMemory> *p, ParticleSys<HostMemory> *tmpPs,BoundingBox *box);
void update_pList_withSort_fast(ParticleSys<HostMemory> *p, ParticleSys<HostMemory> *tmpPs,BoundingBox *box);

void calc_BoundingBoxLimits(BoundingBox *box, TriangleMesh *mesh, cJSON *json_inlet_type);

void initialize_BoundingBox(ParticleSys<HostMemory> *p, BoundingBox *box,TriangleMesh* mesh, cJSON *json_inlet, int isGPUon);
void free_BoundingBox(BoundingBox *box, int isGPUon);

/* ============================
 * device related functions
 * ============================*/

__device__ __forceinline__ void d_sort_neighborlist(int *neiList, int startInd, int numNei);
__global__ void d_swap_ps(DeviceParticleGroup *p, DeviceParticleGroup *tmp);
__global__ void dk_build_cellCount(DeviceParticleGroup* p, DeviceBoundingBox* box);
__global__ void dk_build_pList(DeviceParticleGroup* p, DeviceBoundingBox* box);

__global__ void k_update_neighborlist_endsort(Common *p,DeviceBoundingBox *box);
__global__ void k_update_neighborlist_endsort(DeviceParticleGroup *p,DeviceBoundingBox *box);
 __global__ void k_update_neighborlist(DeviceParticleGroup *p,DeviceBoundingBox *box);

void d_update_pList_withSort(ParticleSystem *p,ParticleSystem *tmpPs, BoundingBox *box,int gridSize, int blockSize);

__device__ int d_calcCellId(DeviceParticleGroup* p,int i, DeviceBoundingBox* box);
__device__ int d_calcCellId(Common* p,int i, DeviceBoundingBox* box);
void d_update_pList(ParticleSystem *p, BoundingBox *box,int gridSize, int blockSize);
void d_update_pList(ParticleSys<DeviceMemory> *p,int N, BoundingBox *box,int gridSize, int blockSize);
void copyToDeviceBox(BoundingBox *box, int N);

