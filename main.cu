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
#include "TriangleMesh.h"
#include "solver_output.h"
#include "Vec3.h"
#define DIM 3
#define OUTPUT 1
#define NONDIM 1
#define MAX_NEIGHBOR 50





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
    /* copy particle system for later swap */
    ParticleSystem tmpPs;

    BoundingBox box;

    /* =========== parameters ============= */
    double r = 0.01;
    double res = 0.3; //CoR
    double density = 1000;
    double m = density*3.14*r*r*r*4./3.;
    double k = 1e6;
    double mu = 0.3;

    ps.N = 1000;
    tmpPs.N = ps.N;

    /* read triangles */
    printf("\n Loading Triangles\n");
    const char* trianglesDir = "geometry/box.stl";
    TriangleMesh mesh;
    load_ascii_stl_double(trianglesDir,&mesh);
    printf("Loading Triangles done!\n");


    /*============ BoundingBox and walls ================== */
    ps.walls.N = 5;
    tmpPs.walls.N = ps.walls.N;

    double minx = mesh.gminx;
    double miny = mesh.gminy;
    double minz = mesh.gminz;

    double maxx = mesh.gmaxx;
    double maxy = 7.0;
    double maxz = mesh.gmaxz;

    printf("Bounding box min (x,y,z)= %f %f %f\n",minx ,miny, minz);
    printf("Bounding box max (x,y,z)= %f %f %f\n",maxx ,maxy, maxz);
       
    printf("allocating memory\n");
    ps.MAX_NEI=MAX_NEIGHBOR;
    tmpPs.MAX_NEI=ps.MAX_NEI;
    ps.mu = mu;

    allocateMemory(&ps);
    allocateMemory(&tmpPs);
    printf("allocating memory done\n");

    printf("initalizing particles\n");
    initializeParticles(&ps,r,m,k,res);
    initializeParticles(&tmpPs,r,m,k,res);
    printf("initalizing particles done\n");



    /* give gravity */
    ps.g[0] = 0.;
    ps.g[1] = -9.81;
    ps.g[2] = 0.;
    printf("%f %f %f\n", ps.g[0],ps.g[1],ps.g[2]);

    /* set time step */
    double dt = 1e-5;
    double out_time = 0.01;
    double end_time = 2.0;
    int outStep = round(out_time/dt);
    printf("Outstep = %d\n",outStep);

    ps.dt=dt;
    

    printf("\nInitializing the Bounding Box\n");
    initialize_BoundingBox(&ps, &box, &mesh, minx, maxx, miny, maxy, minz, maxz);
    printf("Initializing the Bounding Box Done!!\n");

    printf("\n Updating triangle list\n");
    update_tList(&box, &mesh);
    printf("Updating triangle list done!\n");



    /* initialize Walls */
    /*
    double walllist[ps.walls.N*DIM]={0.,1.,0.,
    -1.,0.,0.,
    1.,0.,0.,
    0.,0.,1.0,
    0.,0.,-1.};;;

    double wallPoint[ps.walls.N*DIM]={0.,0.,0.,
    0.2,0.,0.,
    -0.2,0.,0.,
    0.,0.,-0.2,
    0.,0.,0.2};

    
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
    */

    /* non dimensionalize */
    #if NONDIM
        printf("nondimensionalizing ...\n");
        nondimensionalize(&ps,&box,&mesh);
        printf("nondimensionalizing done \n");
        printf("\n");
        printf("g after nondim is %f %f %f \n",ps.g[0],ps.g[1],ps.g[2]);
        printf("m[0]:%f r[0]:%f dt:%f \n",ps.m[0],ps.r[0],ps.dt);
        printf("\n");
    #endif

    #if USE_GPU
        printf("copying memory to device\n");
        copyToDevice(&ps);
        copyToDevice(&tmpPs);
        copyToDeviceBox(&box,&ps);
        deviceMallocCopyTriangleMesh(&mesh);
        printf("copying memory to device done\n");
    #endif



    int steps = (int)(end_time/dt);

    #if USE_GPU
    int blockSize = 256;
    int gridSize = (ps.N + blockSize - 1) / blockSize;
    printf("grid=%d, block=%d\n", gridSize, blockSize);

    check_g_kernel<<<1, 1>>>(ps.d_groupPtr,mesh.d_meshPtr);
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
    for (int step = 1; step < steps; step++)
    {
        #if USE_GPU
        /* GPU */

        /* if want naive collision*/
            //integrateKernel<<<gridSize, blockSize>>>(ps.d_group);

            //device_dem(&ps, &box, gridSize, blockSize);
            device_dem_triangles(&ps, &box, &mesh,gridSize, blockSize);
            //device_dem_withSort(&ps, &tmpPs,&box, gridSize, blockSize,step);
            #if OUTPUT
                if (step % outStep == 0)
                {
                    cudaDeviceSynchronize();
                    copyFromDevice(&ps);
                    #if NONDIM
                        write_frame_bin(outdir,step,&ps,ps.length_factor);
                    #else
                        write_frame_bin(outdir,step,&ps,1.0);
                    #endif

                    cudaEventRecord(now);
                    cudaEventSynchronize(now);

                    cudaEventElapsedTime(&ms, start, now);
                    printf("Output step %d, current time: %f, GPU time: %f s\n", step, (step)*dt,ms/1000.0f);
                }
            #endif
        #else

                /* CPU */
                //cpu_dem_nosort(&ps, &tmpPs, &box);
                //cpu_dem_sort(&ps, &tmpPs, &box, step);
                cpu_dem_sort_triangles(&ps, &tmpPs, &box,&mesh, step);
                checkOoB(&ps,&box);

            #if OUTPUT
            if (step % outStep == 0)
            {
              //  writeParticlesVTKBinary(&ps, step);
                    #if NONDIM
                        //writeParticlesDimensionalizeVTK(&ps, step);
                        write_frame_bin(outdir,step,&ps,ps.length_factor);
                    #else
                        write_frame_bin(outdir,step,&ps,1.0);
                    #endif
                        end = clock();
                        ms = (double)(end-start)*1000./CLOCKS_PER_SEC;;
                        printf("Output step %d,current time: %f s, CPU time: %f s\n", step,(step)*dt,ms/1000.);
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
    freeMemory(&tmpPs);

    free_TriangleMesh(&mesh);
    free_BoundingBox(&box);

    #if USE_GPU
    cudaDeviceReset();
    #endif
    return 0;
}
