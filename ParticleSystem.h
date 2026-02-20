#ifndef _PARTICLESYSTEM_H_
#define _PARTICLESYSTEM_H_
#define DIM 3
#include "BoundingBox.h"

/*
============================================================
GPU専用構造体（SoA）
============================================================
*/
typedef struct {
    int N;

    /* assumed to be normalized */
    double *n;

    double *d;
} DeviceWallGroup;

typedef struct {
    int N;
    double* x;
    double* v;
    double* a;
    double* f;
    double* r;
    double* rsq;
    double* invr;
    double* m;
    double* sqrtm;
    double* invm;
    double* k;
    double* etaconst;
    double* g;
    double dt;
    DeviceWallGroup walls;
} DeviceParticleGroup;



/*
============================================================
CPU管理構造体
============================================================
*/

typedef struct {
    int N;

    /* assumed to be normalized */
    double *n;

    double *d;

} WallGroup;

typedef struct ParticleSystem{
    int N;

    double* x;
    double* v;
    double* r;
    double* rsq;
    double* invr;
    double* a;
    double* f;
    double* k;
    double* m;
    double* sqrtm;
    double* invm;
    double* etaconst;
    double* g;

    int* cellId;
    int* cellx; //coordinate index in structured cell
    
    double dt;
    /* nondimensionalize factor */
    double time_factor;
    double length_factor;
    double mass_factor;


    WallGroup walls;

    DeviceParticleGroup d_group;
    DeviceWallGroup d_walls;

} ParticleSystem;


/*
============================================================
free
============================================================
*/
void freeMemory(ParticleSystem* ps);

/*
============================================================
Host→Device転送
============================================================
*/
void copyToDevice(ParticleSystem* ps);


/*
============================================================
Device→Host転送
============================================================
*/
void copyFromDevice(ParticleSystem* ps);

/*
============================================================
メモリ確保
============================================================
*/
void allocateMemory(ParticleSystem* ps);

/*
============================================================
初期化
============================================================
*/
void initializeParticles(ParticleSystem* ps,double r,double m, double k, double res);

void nondimensionalize(ParticleSystem* ps, BoundingBox* box);



#endif
