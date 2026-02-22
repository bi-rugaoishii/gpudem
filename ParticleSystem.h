#pragma once
#ifndef _PARTICLESYSTEM_H_
#define _PARTICLESYSTEM_H_
#define DIM 3
#include "Vec3.h"
#include "BoundingBox.h"

/*
============================================================
GPU専用構造体（SoA）
============================================================
*/
typedef struct DeviceWallGroup{
    int N;

    /* assumed to be normalized */
    double *n;

    double *d;
} DeviceWallGroup;

typedef struct DeviceParticleGroup{
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

    int MAX_NEIGHBOR;
    int* cellId;
    int* cellx; //coordinate index in structured cell

    DeviceWallGroup walls;
} DeviceParticleGroup;



/*
============================================================
CPU管理構造体
============================================================
*/

typedef struct WallGroup{
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

    double* angv; /* angular velocity */
    double* anga; /* angular accelartion */
    double* moi; /* moment of intertia */
    double* invmoi;
    double* mom; /* torque */

    double* etaconst;
    double* g;

    /* ======== tangential force related =========*/
     /* history of delta tangent */
    double* deltHisx; 
    double* deltHisy;
    double* deltHisz;

    double* deltHisxWall;
    double* deltHisyWall;
    double* deltHiszWall;

    int* indHis; /* Index of collided particles */ 
    int* indHisWall; /* Index of collided particles */ 
    int* isContact; /* flag if a particle has contacted with the neighbor in history*/ 
    int* isContactWall;
    int MAX_NEI; //sets maximum number of neighbors
                 
    int* numCont; /* number of contact */
    int* numContWall; /* number of contact with walls */



    int* cellId;
    int* cellx; //coordinate index in structured cell
    
                 //
    double dt;
    double mu; //friction
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
