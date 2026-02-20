// stage4_5_device_function_version.cu
#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#include <math.h>
#include <chrono>
#define DIM 3
#define USE_GPU 1
#define OUTPUT 1

/*
============================================================
GPU専用構造体（SoA）
============================================================
*/
typedef struct {
    int N;
    double* x;
    double* v;
    double* r;
    double g[3]={0.,-9.81,0.};
} DeviceParticleGroup;


/*
============================================================
CPU管理構造体
============================================================
*/
typedef struct {
    int N;

    double* h_x;
    double* h_v;
    double* h_r;

    DeviceParticleGroup d_group;

} ParticleSystem;

/*
============================================================
cpu functions
============================================================
*/
void integrateCPU(ParticleSystem* ps,
                  double dt,
                  double restitution)
{
    for (int i = 0; i < ps->N; i++)
    {
        // velocity
        ps->h_v[i*DIM+0] += ps->d_group.g[0] * dt;
        ps->h_v[i*DIM+1] += ps->d_group.g[1] * dt;
        ps->h_v[i*DIM+2] += ps->d_group.g[2] * dt;

        // position
        ps->h_x[i*DIM+0] += ps->h_v[i*DIM+0] * dt;
        ps->h_x[i*DIM+1] += ps->h_v[i*DIM+1] * dt;
        ps->h_x[i*DIM+2] += ps->h_v[i*DIM+2] * dt;

        // floor collision
        if (ps->h_x[i*DIM+1] < ps->h_r[i])
        {
            ps->h_x[i*DIM+1] = ps->h_r[i];
            ps->h_v[i*DIM+1] = -ps->h_v[i*DIM+1] * restitution;
        }
    }
}

/*
============================================================
__device__ 関数群
============================================================
*/

/* 速度更新（重力適用） */
__device__ __forceinline__
void updateVelocity(DeviceParticleGroup p,
                    int i,
                    double dt)
{
    p.v[i*DIM+0] += p.g[0] * dt;
    p.v[i*DIM+1] += p.g[1] * dt;
    p.v[i*DIM+2] += p.g[2] * dt;
}


/* 位置更新（オイラー法） */
__device__ __forceinline__
void updatePosition(DeviceParticleGroup p,
                    int i,
                    double dt)
{
    p.x[i*DIM+0] += p.v[i*DIM+0] * dt;
    p.x[i*DIM+1] += p.v[i*DIM+1] * dt;
    p.x[i*DIM+2] += p.v[i*DIM+2] * dt;
}


/* 床衝突処理 */
__device__ __forceinline__
void resolveFloorCollision(DeviceParticleGroup p,
                           int i,
                           double restitution)
{
    if (p.x[i*DIM+1] < p.r[i])
    {
        p.x[i*DIM+1] = p.r[i];
        p.v[i*DIM+1] = -p.v[i*DIM+1] * restitution;
    }
}


/*
============================================================
カーネル
============================================================
*/
__global__ void integrateKernel(DeviceParticleGroup p,
                                double dt,
                                double restitution)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= p.N) return;

    // 処理A
    updateVelocity(p, i,  dt);

    // 処理B
    updatePosition(p, i, dt);

    // 処理C
    resolveFloorCollision(p, i, restitution);
}


/*
============================================================
初期化
============================================================
*/
void initializeParticles(ParticleSystem* ps)
{
    for (int i = 0; i < ps->N; i++)
    {
        ps->h_x[i*DIM+0] = (double)rand() / RAND_MAX * 1.0;
        ps->h_x[i*DIM+1] = (double)rand() / RAND_MAX * 1.0 + 0.5;
        ps->h_x[i*DIM+2] = (double)rand() / RAND_MAX * 1.0;
        ps->h_v[i*DIM+0] = 0.0;
        ps->h_v[i*DIM+1] = 0.0;
        ps->h_v[i*DIM+2] = 0.0;
        ps->h_r[i] = 0.005;
    }
}


/*
============================================================
メモリ確保
============================================================
*/
void allocateMemory(ParticleSystem* ps)
{
    size_t size = ps->N * sizeof(double);

    // Host
    ps->h_x  = (double*)malloc(size*DIM);
    ps->h_v = (double*)malloc(size*DIM);
    ps->h_r = (double*)malloc(size);

    // Device
    cudaMalloc(&ps->d_group.x, size*DIM);
    cudaMalloc(&ps->d_group.v, size*DIM);
    cudaMalloc(&ps->d_group.r, size);

    ps->d_group.N = ps->N;
}


/*
============================================================
Host→Device転送
============================================================
*/
void copyToDevice(ParticleSystem* ps)
{
    size_t size = ps->N * sizeof(double);

    cudaMemcpy(ps->d_group.x,  ps->h_x,  size*DIM, cudaMemcpyHostToDevice);
    cudaMemcpy(ps->d_group.v, ps->h_v, size*DIM, cudaMemcpyHostToDevice);
    cudaMemcpy(ps->d_group.r, ps->h_r, size, cudaMemcpyHostToDevice);
}


/*
============================================================
Device→Host転送
============================================================
*/
void copyFromDevice(ParticleSystem* ps)
{
    size_t size = ps->N * sizeof(double);

    cudaMemcpy(ps->h_x, ps->d_group.x, size*DIM, cudaMemcpyDeviceToHost);
}

/*
============================================================
free
============================================================
*/
void freeMemory(ParticleSystem* ps)
{
    /* host*/
    free(ps->h_x);
    free(ps->h_v);
    free(ps->h_r);

    #if USE_GPU
        /* device */
        cudaFree(ps->d_group.x);
        cudaFree(ps->d_group.v);
        cudaFree(ps->d_group.r);
    #endif
}


/*
============================================================
出力
============================================================
*/
void writeParticles(ParticleSystem* ps, int step)
{
    char filename[256];
    sprintf(filename, "output_%05d.csv", step);

    FILE* fp = fopen(filename, "w");

    for (int i = 0; i < ps->N; i++)
    {
        fprintf(fp, "%f,%f %f\n", ps->h_x[i*DIM+0], ps->h_x[i*DIM+1],ps->h_x[i*DIM+2]);
    }

    fclose(fp);
}
/*
============================================================
VTK Binary出力（ParaView用・高速）
============================================================
*/

/* ============================================================
   エンディアン変換（little → big）
============================================================ */
void swapBytes(void* data, size_t size)
{
    char* p = (char*)data;
    size_t i;
    for (i = 0; i < size/2; i++)
    {
        char tmp = p[i];
        p[i] = p[size-1-i];
        p[size-1-i] = tmp;
    }
}
void writeParticlesVTKBinary(ParticleSystem* ps, int step)
{
    char filename[256];
    sprintf(filename, "vtkresults/particles_%05d.vtk", step);

    FILE* fp = fopen(filename, "wb");
    if (fp == NULL)
    {
        printf("Failed to open VTK file.\n");
        return;
    }

    /* ---------- ヘッダ（ASCII） ---------- */
    fprintf(fp, "# vtk DataFile Version 3.0\n");
    fprintf(fp, "DEM particle output\n");
    fprintf(fp, "BINARY\n");
    fprintf(fp, "DATASET POLYDATA\n");

    /* ---------- POINTS ---------- */
    fprintf(fp, "POINTS %d double\n", ps->N);

    int i, k;
    double value;

    for (i = 0; i < ps->N; i++)
    {
        for (k = 0; k < 3; k++)
        {
            value = ps->h_x[i*DIM + k];
            swapBytes(&value, sizeof(double));
            fwrite(&value, sizeof(double), 1, fp);
        }
    }

    fprintf(fp, "\n");

    /* ---------- VERTICES ---------- */
    fprintf(fp, "VERTICES %d %d\n", ps->N, ps->N*2);

    int one = 1;
    int index;

    for (i = 0; i < ps->N; i++)
    {
        int tmp = one;
        swapBytes(&tmp, sizeof(int));
        fwrite(&tmp, sizeof(int), 1, fp);

        index = i;
        swapBytes(&index, sizeof(int));
        fwrite(&index, sizeof(int), 1, fp);
    }

    fprintf(fp, "\n");

    /* ---------- 半径 ---------- */
    fprintf(fp, "POINT_DATA %d\n", ps->N);
    fprintf(fp, "SCALARS radius double 1\n");
    fprintf(fp, "LOOKUP_TABLE default\n");

    for (i = 0; i < ps->N; i++)
    {
        value = ps->h_r[i];
        swapBytes(&value, sizeof(double));
        fwrite(&value, sizeof(double), 1, fp);
    }

    fclose(fp);
}
/*
============================================================
VTK出力（ParaView用）
============================================================
*/
void writeParticlesVTK(ParticleSystem* ps, int step)
{
    char filename[256];

    /* フォルダ付き出力 */
    sprintf(filename, "vtkresults/particles_%05d.vtk", step);

    FILE* fp = fopen(filename, "w");
    if (fp == NULL)
    {
        printf("Failed to open VTK file.\n");
        return;
    }

    /* --- VTKヘッダ --- */
    fprintf(fp, "# vtk DataFile Version 3.0\n");
    fprintf(fp, "DEM particle output\n");
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET POLYDATA\n");

    /* --- 粒子位置 --- */
    fprintf(fp, "POINTS %d double\n", ps->N);
    for (int i = 0; i < ps->N; i++)
    {
        fprintf(fp, "%lf %lf %lf\n",
                ps->h_x[i*DIM+0],
                ps->h_x[i*DIM+1],
                ps->h_x[i*DIM+2]);
    }

    /* --- 各点を頂点セルとして定義 --- */
    fprintf(fp, "\nVERTICES %d %d\n", ps->N, ps->N * 2);
    for (int i = 0; i < ps->N; i++)
    {
        fprintf(fp, "1 %d\n", i);
    }

    /* --- 半径データ --- */
    fprintf(fp, "\nPOINT_DATA %d\n", ps->N);
    fprintf(fp, "SCALARS radius double 1\n");
    fprintf(fp, "LOOKUP_TABLE default\n");

    for (int i = 0; i < ps->N; i++)
    {
        fprintf(fp, "%lf\n", ps->h_r[i]);
    }

    fclose(fp);
}

/*
============================================================
メイン
============================================================
*/
int main()
{
    ParticleSystem ps;
    ps.N = 1000000;

    allocateMemory(&ps);
    initializeParticles(&ps);
    #if USE_GPU
        copyToDevice(&ps);
    #endif

    double dt = 1e-5;
    double restitution = 0.3;
    double out_time = 0.05;
    double end_time = 2.0;
    int outStep = (int)(out_time/dt);


    int steps = (int)(end_time/dt);

    #if USE_GPU
    int blockSize = 256;
    int gridSize = (ps.N + blockSize - 1) / blockSize;

        cudaEvent_t start, stop;
        float ms;
        cudaEventCreate(&start);
        cudaEventCreate(&stop);
        cudaEventRecord(start);
    #else
        clock_t start, end;
        double ms;
        start = clock();
    #endif

    for (int step = 0; step < steps; step++)
    {
        #if USE_GPU
            integrateKernel<<<gridSize, blockSize>>>(
                ps.d_group,
                dt,
                restitution);
            #if OUTPUT
                if (step % outStep == 0)
                {
                    copyFromDevice(&ps);
                    writeParticlesVTKBinary(&ps, step);
                    printf("Output step %d\n", step);
                }
            #endif
        #else
            integrateCPU(&ps,dt,restitution);
            #if OUTPUT
                if (step % outStep == 0)
                {
                    writeParticlesVTKBinary(&ps, step);
                    printf("Output step %d\n", step);
                }
            #endif
        #endif
    }

    #if USE_GPU
        cudaEventRecord(stop);
        cudaEventSynchronize(stop);
        cudaEventElapsedTime(&ms,start, stop);
        printf("GPU time: %f s\n",ms/1000.);
    #else
        end = clock();
        ms = (double)(end-start)*1000./CLOCKS_PER_SEC;;
        printf("CPU time: %f s\n",ms/1000.);
    #endif

        freeMemory(&ps);
        cudaDeviceReset();
    return 0;
}
