#pragma once
#ifndef _DEVICE_DEM_H_
#define _DEVICE_DEM_H_

#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#include <math.h>
#include <chrono>
#include "ParticleSystem.h"
#include "BoundingBox.h"
#include "BVH.h"
#include "Vec3.h"
#include "ContactCache.h"
#include "TriangleContactCache.h"


/* =========== triangle related ====== */
__device__ __forceinline__
TriangleContactCache d_dist_triangle(DeviceParticleGroup* ps, int i, DeviceTriangleMesh* mesh, int j);

__device__ __forceinline__ void d_particle_collision_verlet(DeviceParticleGroup* p, int i ,DeviceBoundingBox *box);

__device__ __forceinline__
void d_wall_collision_triangles(DeviceParticleGroup* p,int i,DeviceBoundingBox *box, DeviceTriangleMesh* mesh);

__device__ __forceinline__
void d_wall_collision_verlet(DeviceParticleGroup* p,int i,DeviceTriangleMesh* mesh);


__device__ __forceinline__
void updateAcceleration(DeviceParticleGroup* p,int i);

__device__ __forceinline__
void updateVelocity(DeviceParticleGroup* p,int i);


/* 位置更新（オイラー法） */
__device__ __forceinline__
void updatePosition(DeviceParticleGroup* p,
                    int i);


/* 床衝突処理 */
__device__ __forceinline__
void resolveFloorCollision(DeviceParticleGroup* p,
                           int i,
                           double restitution);

__device__ __forceinline__
ContactCache d_calc_normal_force_wall(DeviceParticleGroup *p,int i,int j,Vec3 n,double delMag);

__device__ __forceinline__
ContactCache d_calc_normal_force(DeviceParticleGroup* p,int i,int j,Vec3 n,double delMag);

__device__ __forceinline__
void d_calc_tangential_force_wall(DeviceParticleGroup *p,int i,int j,ContactCache c);

__device__ __forceinline__
void d_calc_tangential_force(DeviceParticleGroup *p,int i,int j,ContactCache c);

__device__ __forceinline__
void d_update_history_wall(DeviceParticleGroup *p,int i);

__device__ __forceinline__
void d_wall_collision_naive(DeviceParticleGroup* ps,int i);

__device__ __forceinline__
void d_update_history(DeviceParticleGroup *p,int i);


__device__ __forceinline__
void d_particle_collision_cell_linked(DeviceParticleGroup* p, int i, DeviceBoundingBox* box);


__device__ __forceinline__
void d_particle_collision_cell_linked_noVec(DeviceParticleGroup* p, int i, DeviceBoundingBox* box);

__device__ __forceinline__
void d_particle_collision_naive(DeviceParticleGroup* ps, int i);
/* ============== check Out of Bounds ============  */
__global__ void dk_checkOoB(DeviceParticleGroup *p, DeviceBoundingBox* box);


/*
============================================================
カーネル
============================================================
*/
__global__ void integrateKernel(DeviceParticleGroup* p);

__global__ void check_g_kernel(DeviceParticleGroup* p,DeviceTriangleMesh *mesh);

__global__ void k_integrate(DeviceParticleGroup* p);

__global__ void k_shouldRefreshNeighborList(DeviceParticleGroup *p, DeviceBoundingBox* box);

__global__ void k_collision_verlet_verlet(DeviceParticleGroup* p, DeviceBoundingBox* box,DeviceTriangleMesh* mesh);

__global__ void k_collision_triangle(DeviceParticleGroup* p, DeviceBoundingBox* box,DeviceTriangleMesh* mesh);

__global__ void k_collision(DeviceParticleGroup* p, DeviceBoundingBox* box);

/*
   ======================================================
   main routine 
   ======================================================
*/

void device_dem_verlet_verlet(ParticleSystem *p, BoundingBox *box,TriangleMesh *mesh, BVH *bvh, int gridSize, int blockSize);

void device_dem_verlet_triangles(ParticleSystem *p, BoundingBox *box,TriangleMesh *mesh, int gridSize, int blockSize);


void device_dem_triangles(ParticleSystem *p, BoundingBox *box,TriangleMesh *mesh, int gridSize, int blockSize);

void device_dem(ParticleSystem *p, BoundingBox *box, int gridSize, int blockSize);
#endif
