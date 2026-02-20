#ifndef _CPU_DEM_H_
#define _CPU_DEM_H_
#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#include <math.h>
#include <chrono>
#include "ParticleSystem.h"
#define DIM 3

void integrateCPU(ParticleSystem* ps);
void particle_collision_naive(ParticleSystem* ps);
inline void velocity(ParticleSystem* ps);
inline void position(ParticleSystem* ps);
#endif
