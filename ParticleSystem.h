#pragma once
#include "hardCodedParameters.h"
#include <stdint.h>
#include <cstdlib>
#include <cuda_runtime.h>
#include "Vec3.h"
#include "BoundingBox.h"
#include "TriangleMesh.h"
#include "cJSON.h"
struct TriangleMesh;
struct BoundingBox;
/* ==================== Memory Structs =================== */
struct HostMemory{
    template <typename T>
    static T* allocate(int N){
            printf("N = %d\n", N);
            printf("MAX_NEI = %d\n", MAX_NEI);
            printf("alloc size = %zu\n", (size_t)N * MAX_NEI * sizeof(double)); 
        return static_cast<T*>(malloc(sizeof(T)*N));
    }

    template <typename T>
    static void deallocate(T* ptr){
        free(ptr);
    }
};

struct DeviceMemory{
    template <typename T>
        static T* allocate(int N){
            T* ptr;
            printf("N = %d\n", N);
            printf("MAX_NEI = %d\n", MAX_NEI);
            printf("alloc size = %zu\n", (size_t)N * MAX_NEI * sizeof(double)); 
            cudaMalloc((void**)&ptr, sizeof(T)*N);
            return ptr;
        }

    template <typename T>
        static void deallocate(T* ptr){
            cudaFree(ptr);
        }

};

/* 
   =======================================
   interface 
   =======================================
   */

typedef struct Parameters{
    double dt;
    double mu; //friction
    int N;

    /* nondimensionalize factor */
    double time_factor;
    double length_factor;
    double mass_factor;
}Parameters;

typedef struct Common{
    #define MEMBER(type,name,Np,SAVE_FLAG) type* name;
    #include "ParticleSystemMember_common.def"
    #undef MEMBER
}Common;

typedef struct HostOnly{
}HostOnly;

typedef struct DeviceOnly{
    /* for sorting */
    void* tmp_storage;
    size_t tmp_bytes;
    int* refreshVerletFlag;
}DeviceOnly;


template <typename Memory>// typename used later to distinguish gpu and cpu
struct ParticleSys;

template <>
struct ParticleSys<DeviceMemory>{
    Parameters parameters;
    Common p;
    DeviceOnly deviceOnly;
};

template <>
struct ParticleSys<HostMemory>{
    Parameters parameters;
    Common p;
    HostOnly hostOnly;
};

/* ==================== Memory Allocation =================== */
void allocate(ParticleSys<HostMemory> *ps);
void allocate(ParticleSys<DeviceMemory> *ps);


void deallocate(ParticleSys<HostMemory> *ps);
void deallocate(ParticleSys<DeviceMemory> *ps);

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

    int* isActive; /* flag if particle is OoB */ 

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
    int* indHisWallNow; /* Index of collided walls in the current step*/ 

    int* indHisVorENow; /* number of contact with edge or vertex*/

    int* isContact; /* flag if a particle has contacted with the neighbor in history*/ 
    int* isContactWall;

    int* numCont; /* number of contact */
    int* numContWall; /* number of contact with walls */



    int* cellId;
    int* cellx; //coordinate index in structured cell

    //
    double dt;
    double mu; //friction

    int* pId;
    uint32_t *mortonKey;
    int* mortonOrder;

    uint32_t *tmpMortonKey;
    int* tmpMortonOrder;

    /* for sorting */
    void* tmp_storage;
    size_t tmp_bytes;


    /* ========== verlet list related ========= */
    double* refx;
    double* refy;
    double* refz; 
    int *neiList; /* neighbor list */
    int *numNei; /* neighbor list */
    int *neiListWall; /* neighbor list */
    int *numNeiWall; /* neighbor list */
    int *refreshVerletFlag;

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
    int* isActive; /* flag if particle is OoB */ 

    /* ======== tangential force related =========*/
    /* history of delta tangent */
    double* deltHisx; 
    double* deltHisy;
    double* deltHisz;

    double* deltHisxWall;
    double* deltHisyWall;
    double* deltHiszWall;

    int* indHis; /* Index of collided particles */ 
    int* indHisWall; /* Index of collided walls */ 
    int* indHisWallNow; /* Index of collided walls in the current step*/ 
    int* isContact; /* flag if a particle has contacted with the neighbor in history*/ 
    int* isContactWall;

    int* numCont; /* number of contact */
    int* numContWall; /* number of contact with walls */

    int* indHisVorENow; /* number of contact with edge or vertex*/


    int* cellId;
    int* cellx; //coordinate index in structured cell

    //
    double dt;
    double mu; //friction

    int* pId;
    uint32_t *mortonKey;
    int* mortonOrder;

    uint32_t *tmpMortonKey;
    int* tmpMortonOrder;
    /* nondimensionalize factor */
    double time_factor;
    double length_factor;
    double mass_factor;

    /* ========== verlet list related ========= */
    double* refx;
    double* refy;
    double* refz; 
    int *neiList; /* neighbor list */
    int *numNei; /* neighbor list */
    int *neiListWall; /* neighbor list */
    int *numNeiWall; /* neighbor list */


    WallGroup walls;

    DeviceParticleGroup d_group;
    DeviceParticleGroup *d_groupPtr;
    DeviceWallGroup *d_wallsPtr;

} ParticleSystem;


/*
   ============================================================
   free
   ============================================================
   */
void freeMemory(ParticleSystem* ps, int isGPUon);

/*
   ============================================================
   Host→Device転送
   ============================================================
   */
void copyToDevice(ParticleSystem* ps);
void copyToDevice(ParticleSys<HostMemory> *ps, ParticleSys<DeviceMemory> *d_ps);


/*
   ============================================================
   Device→Host転送
   ============================================================
   */
void copyFromDevice(ParticleSystem* ps);
void copyFromDevice(ParticleSys<DeviceMemory>* d_ps, ParticleSys<HostMemory>* ps);

/*
   ============================================================
   メモリ確保
   ============================================================
   */
void allocateMemory(ParticleSystem* ps, int isGPUon);

/*
   ============================================================
   初期化
   ============================================================
   */
void initializeTmpParticles(ParticleSys<HostMemory>* ps,cJSON *json_inlet, double r,double m,double k,double res);
void initializeParticles(ParticleSys<HostMemory>* ps,cJSON *json_inlet, double r,double m, double k, double res);

void nondimensionalize(ParticleSys<HostMemory>* ps, BoundingBox* box, TriangleMesh *mesh);



