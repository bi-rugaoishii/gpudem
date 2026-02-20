#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#include <math.h>
#include <chrono>
#include "ParticleSystem.h"

#ifndef _DEVICE_DEM_H_
#define _DEVICE_DEM_H_

__device__ __forceinline__
void updateAcceleration(DeviceParticleGroup *p,int i);
__device__ __forceinline__
void updateVelocity(DeviceParticleGroup *p,int i);


/* 位置更新（オイラー法） */
__device__ __forceinline__
void updatePosition(DeviceParticleGroup *p,
                    int i);


/* 床衝突処理 */
__device__ __forceinline__
void resolveFloorCollision(DeviceParticleGroup p,
                           int i,
                           double restitution);


/*
============================================================
カーネル
============================================================
*/
__global__ void integrateKernel(DeviceParticleGroup p);

__global__ void check_g_kernel(DeviceParticleGroup p);


#endif
