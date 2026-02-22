#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#include <math.h>
#include <chrono>
#include "ParticleSystem.h"
#include "BoundingBox.h"
#include "device_dem.h"
#include "cpu_dem.h"
#include "output.h"
#include "solver_output.h"
#include "Vec3.h"
#define DIM 3
#define OUTPUT 1
#define NONDIM 1





/*
============================================================
メイン
============================================================
*/
int main()
{
    setvbuf(stdout,NULL,_IOLBF,0);
    setvbuf(stderr,NULL,_IONBF,0);

    ParticleSystem ps;
    BoundingBox box;
    ps.N = 1000;
    ps.walls.N = 5;
    double r = 0.01;
    double res = 0.3; //CoR
    double density = 1000;
    double m = density*3.14*r*r*r*4./3.;
    double k = 1e6;

    double minx = 0.;
    double miny = 0.;
    double minz = 0.;

    double maxx = 0.5;
    double maxy = 3.0;
    double maxz = 0.5;
       
    printf("allocating memory\n");
    allocateMemory(&ps);
    printf("allocating memory done\n");
    printf("initalizing particles\n");
    initializeParticles(&ps,r,m,k,res);
    printf("initalizing particles done\n");

    initialize_BoundingBox(&ps, &box, minx, maxx, miny, maxy, minz, maxz);

    /* give gravity */
    ps.g[0] = 0.;
    ps.g[1] = -9.81;
    ps.g[2] = 0.;
    printf("%f %f %f\n", ps.g[0],ps.g[1],ps.g[2]);

    /* set time step */
    double dt = 2e-5;
    double out_time = 0.05;
    double end_time = 5.;
    int outStep = (int)(out_time/dt);

    ps.dt=dt;
    
    /* initialize Walls */
    double walllist[ps.walls.N*DIM]={0.,1.,0.,
    -1.,0.,0.,
    1.,0.,0.,
    0.,0.,1.0,
    0.,0.,-1.};

    double wallPoint[ps.walls.N*DIM]={0.,0.,0.,
    0.5,0.,0.,
    0.,0.,0.,
    0.,0.,0.0,
    0.,0.,0.5};

    
    printf("making walls\n");
    for (int i=0; i<ps.walls.N; i++){
    printf("%d\n",i);
        ps.walls.n[i*DIM+0]=walllist[i*DIM+0];
        ps.walls.n[i*DIM+1]=walllist[i*DIM+1];
        ps.walls.n[i*DIM+2]=walllist[i*DIM+2];

        ps.walls.d[i] = 0.;
        ps.walls.d[i] -= (ps.walls.n[i*DIM+0]*wallPoint[i*DIM+0]);
        ps.walls.d[i] -= (ps.walls.n[i*DIM+1]*wallPoint[i*DIM+1]);
        ps.walls.d[i] -= (ps.walls.n[i*DIM+2]*wallPoint[i*DIM+2]);
    }
    printf("making walls done\n");

    /* non dimensionalize */
    #if NONDIM
        printf("nondimensionalizing ...\n");
        nondimensionalize(&ps,&box);
        printf("nondimensionalizing done \n");
    #endif

    #if USE_GPU
        printf("copying memory to device\n");
        copyToDevice(&ps);
        copyToDeviceBox(&box,&ps);
        printf("copying memory to device done\n");
    #endif



    int steps = (int)(end_time/dt);

    #if USE_GPU
    int blockSize = 256;
    int gridSize = (ps.N + blockSize - 1) / blockSize;
    printf("grid=%d, block=%d\n", gridSize, blockSize);

    check_g_kernel<<<1, 1>>>(ps.d_group);
    cudaDeviceSynchronize();
    cudaEvent_t start, stop, now;
    float ms;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventCreate(&now);
    cudaEventRecord(start);
    #else
    clock_t start, end;
    double ms;
    start = clock();
    #endif

    /* initialization for file output */
    
    #if OUTPUT
    const char* outdir = "results";
    solver_output_init(outdir);

    
    #endif

    printf("starting \n");
    for (int step = 0; step < steps; step++)
    {
        #if USE_GPU
        /* if want naive collision*/
        //integrateKernel<<<gridSize, blockSize>>>(ps.d_group);
         
        device_dem(&ps, &box, gridSize, blockSize);
            #if OUTPUT
                if (step % outStep == 0)
                {
                    cudaDeviceSynchronize();
                    copyFromDevice(&ps);
                    //writeParticlesVTKBinary(&ps, step);
                    #if NONDIM
                        //writeParticlesDimensionalizeVTK(&ps, step);
                        write_frame_bin(outdir,step,ps.N,ps.x,ps.r,ps.length_factor);
                    #else
                        writeParticlesVTK(&ps, step);
                    #endif

                    cudaEventRecord(now);
                    cudaEventSynchronize(now);

                    cudaEventElapsedTime(&ms, start, now);
                    printf("Output step %d, current time: %f, GPU time: %f s\n", step, step*dt,ms/1000.0f);
                }
            #endif
        #else
            integrateCPU(&ps,&box);
            #if OUTPUT
            if (step % outStep == 0)
            {
                writeParticlesVTKBinary(&ps, step);
                    #if NONDIM
                        //writeParticlesDimensionalizeVTK(&ps, step);
                        write_frame_bin(outdir,step,ps.N,ps.x,ps.r,ps.length_factor);
                    #else
                        writeParticlesVTK(&ps, step);
                    #endif
                        end = clock();
                        ms = (double)(end-start)*1000./CLOCKS_PER_SEC;;
                        printf("Output step %d,current time: %f s, CPU time: %f s\n", step,step*dt,ms/1000.);
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
    free_BoundingBox(&box);
    cudaDeviceReset();
    return 0;
}
