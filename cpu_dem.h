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
#include "TriangleContactCache.h"
#define DIM 3

/* === triangle related === */
void wall_collision_triangles(ParticleSystem* ps, BoundingBox *box,TriangleMesh* mesh);
void cpu_dem_sort_triangles(ParticleSystem* ps, ParticleSystem *tmpPs, BoundingBox* box,TriangleMesh *mesh, int step);

/* calculated distance between a triangle */
/* i is index of a particle, j is index of a triangle*/
/* returns dist 1e10 if no collision */
TriangleContactCache dist_triangle(ParticleSystem* ps,int i, TriangleMesh* mesh, int j);

void integrateCPU(ParticleSystem *ps, BoundingBox *box);
void particle_collision_naive(ParticleSystem* ps);
void particle_collision_cell_linked(ParticleSystem* ps, BoundingBox *box);
void particle_collision_cell_linked_fastUpdate(ParticleSystem* p, BoundingBox *box);
void particle_collision_cell_linked_withSort_fastupdate(ParticleSystem* p,ParticleSystem* tmpPs, BoundingBox *box);



void cpu_dem_nosort(ParticleSystem* ps, ParticleSystem *tmpPs, BoundingBox* box);
void particle_collision_cell_linked_withSort(ParticleSystem* p,ParticleSystem* tmpPs, BoundingBox *box);

void particle_collision_cell_linked_noVec3(ParticleSystem* ps, BoundingBox *box);

inline ContactCache calc_normal_force(ParticleSystem *p,int i,int j,Vec3 n,double delMag,double dist);


inline ContactCache calc_normal_force_wall(ParticleSystem *p,int i,int j,Vec3 n,double delMag,double dist);

inline void calc_tangential_force_wall(ParticleSystem *p,int i,int j,ContactCache c);

inline void calc_tangential_force(ParticleSystem *p,int i,int j,ContactCache c);

inline void update_history(ParticleSystem *p,int i);
inline void update_history_wall(ParticleSystem *p,int i);

void cpu_dem_sort(ParticleSystem *ps, ParticleSystem *tmpPs, BoundingBox* box, int step);
void wall_collision_naive(ParticleSystem* ps);

inline void velocity(ParticleSystem* ps);
inline void position(ParticleSystem* ps);
#endif
