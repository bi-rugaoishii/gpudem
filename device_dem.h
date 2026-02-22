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
#include "Vec3.h"
#include "ContactCache.h"


__device__ __forceinline__
void updateAcceleration(DeviceParticleGroup p,int i);

__device__ __forceinline__
void updateVelocity(DeviceParticleGroup p,int i);


/* 位置更新（オイラー法） */
__device__ __forceinline__
void updatePosition(DeviceParticleGroup p,
                    int i);


/* 床衝突処理 */
__device__ __forceinline__
void resolveFloorCollision(DeviceParticleGroup p,
                           int i,
                           double restitution);

/* normal force */

inline ContactCache d_calc_normal_force(DeviceParticleGroup p,int i,int j,Vec3 n,double delMag,double dist);

__device__ __forceinline__
void d_particle_collision_cell_linked(DeviceParticleGroup p, int i, DeviceBoundingBox box);

__device__ __forceinline__
ContactCache d_calc_normal_force(ParticleSystem *p,int i,int j,Vec3 n,double delMag,double dist);

__device__ __forceinline__
void d_particle_collision_cell_linked_noVec(DeviceParticleGroup p, int i, DeviceBoundingBox box);

__device__ __forceinline__
void particle_collision_naive(DeviceParticleGroup* ps, int i);

/*
============================================================
カーネル
============================================================
*/
__global__ void integrateKernel(DeviceParticleGroup p);

__global__ void check_g_kernel(DeviceParticleGroup p);

__global__ void k_integrate(DeviceParticleGroup p);

__global__ void k_collision(DeviceParticleGroup p, DeviceBoundingBox box);
/*
   ======================================================
   main routine 
   ======================================================
*/


void device_dem(ParticleSystem *p, BoundingBox *box, int gridSize, int blockSize);
#endif
