#ifndef _CPU_DEM_H_
#define _CPU_DEM_H_
#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#include <math.h>
#include <chrono>
#include "ParticleSystem.h"
#include "BoundingBox.h"
#include "Vec3.h"
#include "ContactCache.h"
#define DIM 3

void integrateCPU(ParticleSystem *ps, BoundingBox *box);
void particle_collision_naive(ParticleSystem* ps);
void particle_collision_cell_linked(ParticleSystem* ps, BoundingBox *box);

void particle_collision_cell_linked_noVec3(ParticleSystem* ps, BoundingBox *box);

inline ContactCache calc_normal_force(ParticleSystem *p,int i,int j,Vec3 n,double delMag,double dist);

void wall_collision_naive(ParticleSystem* ps);

inline void velocity(ParticleSystem* ps);
inline void position(ParticleSystem* ps);
#endif
