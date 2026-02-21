#include "ParticleSystem.h"
#include <stdio.h>


/*
============================================================
free
============================================================
*/
void freeMemory(ParticleSystem* ps)
{
    /* host*/
    free(ps->x);
    free(ps->v);
    free(ps->r);
    free(ps->rsq);
    free(ps->a);
    free(ps->f);
    free(ps->k);
    free(ps->m);
    free(ps->invm);
    free(ps->sqrtm);
    free(ps->invr);
    free(ps->g);
    
    free(ps->cellId);
    free(ps->cellx);

    free(ps->etaconst);
    free(ps->walls.n);
    free(ps->walls.d);

    #if USE_GPU
        cudaFree(ps->d_group.x);
        cudaFree(ps->d_group.v);
        cudaFree(ps->d_group.r);
        cudaFree(ps->d_group.rsq);
        cudaFree(ps->d_group.a);
        cudaFree(ps->d_group.f);
        cudaFree(ps->d_group.k);
        cudaFree(ps->d_group.m);
        cudaFree(ps->d_group.sqrtm);
        cudaFree(ps->d_group.invm);
        cudaFree(ps->d_group.invr);
        cudaFree(ps->d_group.etaconst);
        cudaFree(ps->d_group.g);

        cudaFree(ps->d_group.cellId);
        cudaFree(ps->d_group.cellx);

        cudaFree(ps->d_group.walls.n);
        cudaFree(ps->d_group.walls.d);
    #endif
}


/*
   ============================================================
   Host→Device転送
   ============================================================
 */
void copyToDevice(ParticleSystem *ps)
{
    size_t size = ps->N * sizeof(double);
    size_t size_walls = ps->walls.N * sizeof(double);

    ps->d_group.dt = ps->dt;
    ps->d_group.N = ps->N;
    ps->d_group.walls.N = ps->walls.N;

    cudaMemcpy(ps->d_group.x,  ps->x,  size*DIM, cudaMemcpyHostToDevice);
    cudaMemcpy(ps->d_group.v, ps->v, size*DIM, cudaMemcpyHostToDevice);
    cudaMemcpy(ps->d_group.r, ps->r, size, cudaMemcpyHostToDevice);
    cudaMemcpy(ps->d_group.rsq, ps->rsq, size, cudaMemcpyHostToDevice);
    cudaMemcpy(ps->d_group.a, ps->a, size*DIM, cudaMemcpyHostToDevice);
    cudaMemcpy(ps->d_group.f, ps->f, size*DIM, cudaMemcpyHostToDevice);
    cudaMemcpy(ps->d_group.m, ps->m, size, cudaMemcpyHostToDevice);
    cudaMemcpy(ps->d_group.invm, ps->invm, size, cudaMemcpyHostToDevice);
    cudaMemcpy(ps->d_group.sqrtm, ps->sqrtm, size, cudaMemcpyHostToDevice);
    cudaMemcpy(ps->d_group.invr, ps->invr, size, cudaMemcpyHostToDevice);
    cudaMemcpy(ps->d_group.etaconst, ps->etaconst, size, cudaMemcpyHostToDevice);
    cudaMemcpy(ps->d_group.k, ps->k, size, cudaMemcpyHostToDevice);
    cudaMemcpy(ps->d_group.g, ps->g, sizeof(double)*DIM, cudaMemcpyHostToDevice);

    cudaMemcpy(ps->d_group.walls.n,  ps->walls.n,  size_walls*DIM, cudaMemcpyHostToDevice);
    cudaMemcpy(ps->d_group.walls.d,  ps->walls.d,  size_walls, cudaMemcpyHostToDevice);
}


/*
   ============================================================
   Device→Host転送
   ============================================================
 */
void copyFromDevice(ParticleSystem* ps)
{
    size_t size = ps->N * sizeof(double);

    cudaMemcpy(ps->x, ps->d_group.x, size*DIM, cudaMemcpyDeviceToHost);
}

/*
   ============================================================
   メモリ確保
   ============================================================
 */
void allocateMemory(ParticleSystem* ps)
{
    size_t size = ps->N * sizeof(double);
    size_t size_walls = ps->walls.N * sizeof(double);

    // Host
    ps->x  = (double*)malloc(size*DIM);
    ps->v = (double*)malloc(size*DIM);
    ps->a = (double*)malloc(size*DIM);
    ps->f = (double*)malloc(size*DIM);
    ps->r = (double*)malloc(size);
    ps->rsq = (double*)malloc(size);
    ps->invr = (double*)malloc(size);
    ps->m = (double*)malloc(size);
    ps->sqrtm = (double*)malloc(size);
    ps->invm = (double*)malloc(size);
    ps->k = (double*)malloc(size);
    ps->etaconst = (double*)malloc(size);
    ps->g = (double*)malloc(sizeof(double)*DIM);

    ps->cellId = (int*)malloc(sizeof(int)*ps->N);
    ps->cellx = (int*)malloc(sizeof(int)*DIM*ps->N);

    ps->walls.n = (double*)malloc(size_walls*DIM);
    ps->walls.d = (double*)malloc(size_walls);

    // Device
    #if USE_GPU
    cudaMalloc(&ps->d_group.x, size*DIM);
    cudaMalloc(&ps->d_group.v, size*DIM);
    cudaMalloc(&ps->d_group.a, size*DIM);
    cudaMalloc(&ps->d_group.f, size*DIM);
    cudaMalloc(&ps->d_group.r, size);
    cudaMalloc(&ps->d_group.rsq, size);
    cudaMalloc(&ps->d_group.invr, size);
    cudaMalloc(&ps->d_group.m, size);
    cudaMalloc(&ps->d_group.sqrtm, size);
    cudaMalloc(&ps->d_group.invm, size);
    cudaMalloc(&ps->d_group.k, size);
    cudaMalloc(&ps->d_group.etaconst, size);


    cudaMalloc(&ps->d_group.cellId, sizeof(int)*ps->N);
    cudaMalloc(&ps->d_group.cellx, sizeof(int)*DIM*ps->N);

    cudaMalloc(&ps->d_group.walls.n, size_walls*DIM);
    cudaMalloc(&ps->d_group.walls.d, size_walls);
    cudaMalloc(&ps->d_group.g, sizeof(double)*DIM);
    #endif

}

/*
   ============================================================
   初期化
   ============================================================
 */
void initializeParticles(ParticleSystem* ps,double r,double m,double k,double res)
{
    int max_trials = 1000; // 1粒子あたりの再配置試行回数


    for (int i = 0; i < ps->N; i++){
        int trial = 0;
        int success = 0;
        while (trial < max_trials)
        {
            // ランダム配置
            double x = (double)rand() / RAND_MAX * 0.25 + 0.1;
            double y = (double)rand() / RAND_MAX * 2.0 + 0.5;
            double z = (double)rand() / RAND_MAX * 0.25 + 0.1;

            // 既存粒子との距離チェック
            int overlap = 0;
            for (int j = 0; j < i; j++)
            {
                double dx = x - ps->x[j*DIM+0];
                double dy = y - ps->x[j*DIM+1];
                double dz = z - ps->x[j*DIM+2];
                double d2 = dx*dx + dy*dy + dz*dz;
                if (d2 < 4*r*r) // 半径2倍以上離れているか
                {
                    overlap = 1;
                    break;
                }
            }

            if (!overlap)
            {
                ps->x[i*DIM+0] = x;
                ps->x[i*DIM+1] = y;
                ps->x[i*DIM+2] = z;
                success = 1;
                break;
            }

            trial++;
        }

        if (trial == max_trials){
            printf("!!!couldn't place them!!\n");
        };


        /*
           ps->x[i*DIM+0] = 0.25;
           ps->x[i*DIM+1] = (double)i*1.0+0.03;
           ps->x[i*DIM+2] = 0.25;
         */

        ps->r[i] = r;
        ps->rsq[i] = ps->r[i]*ps->r[i];
        ps->invr[i] = 1./ps->r[i];
        ps->k[i] = k;
        ps->m[i] = m;
        ps->sqrtm[i] = sqrt(m);
        ps->invm[i] = 1./m;
        ps->etaconst[i]=-2.*log(res)*sqrt(ps->k[i]/(3.1415*3.1415+log(res)*log(res)));
        ps->cellId[i] = -1;
    }


}

void nondimensionalize(ParticleSystem* ps, BoundingBox *box){
    /* get nondimensionalize factor */
    ps->time_factor = sqrt(ps->m[0]/ps->k[0]);
    ps->mass_factor = ps->m[0];
    ps->length_factor = ps->r[0];

    double inv_time_factor = 1./ps->time_factor;
    double inv_mass_factor = 1./ps->mass_factor;
    double inv_length_factor = 1./ps->length_factor;
    double inv_sqrtk_factor = 1./sqrt(ps->k[0]);

    ps->dt*=inv_time_factor;

    ps->g[0]*=inv_length_factor*ps->time_factor*ps->time_factor;
    ps->g[1]*=inv_length_factor*ps->time_factor*ps->time_factor;
    ps->g[2]*=inv_length_factor*ps->time_factor*ps->time_factor;

    for (int i=0; i<ps->N; i++){
        ps->x[i*DIM+0]*=inv_length_factor;
        ps->x[i*DIM+1]*=inv_length_factor;
        ps->x[i*DIM+2]*=inv_length_factor;
    }

    for (int i=0; i<ps->N; i++){
        ps->v[i*DIM+0]*=inv_length_factor*ps->time_factor;
        ps->v[i*DIM+1]*=inv_length_factor*ps->time_factor;
        ps->v[i*DIM+2]*=inv_length_factor*ps->time_factor;
    }

    for (int i=0; i<ps->N; i++){
        ps->r[i]*=inv_length_factor;
    }

    for (int i=0; i<ps->N; i++){
        ps->rsq[i]*=inv_length_factor*inv_length_factor;
    }

    for (int i=0; i<ps->N; i++){
        ps->invr[i]*=ps->length_factor;
    }

    for (int i=0; i<ps->N; i++){
        ps->m[i]*=inv_mass_factor;
    }


    for (int i=0; i<ps->N; i++){
        ps->sqrtm[i]*=sqrt(inv_mass_factor);
    }

    for (int i=0; i<ps->N; i++){
        ps->invm[i]*=ps->mass_factor;
    }

    for (int i=0; i<ps->N; i++){
        ps->k[i]*=inv_sqrtk_factor*inv_sqrtk_factor;
    }

    for (int i=0; i<ps->N; i++){
        ps->etaconst[i]*=inv_sqrtk_factor;
    }

    for (int i=0; i<ps->N; i++){
        ps->invr[i]*=ps->length_factor;
    }


    for (int i=0; i<ps->walls.N; i++){
        ps->walls.d[i]*=inv_length_factor;
    }

    /* bounding box*/

    box->minx*=inv_length_factor;
    box->miny*=inv_length_factor;
    box->minz*=inv_length_factor;
    box->maxx*=inv_length_factor;
    box->maxy*=inv_length_factor;
    box->maxz*=inv_length_factor;

    box->dx*=inv_length_factor;
    box->dy*=inv_length_factor;
    box->dz*=inv_length_factor;

    box->invdx*=ps->length_factor;
    box->invdy*=ps->length_factor;
    box->invdz*=ps->length_factor;

    box->rangex*=inv_length_factor;
    box->rangey*=inv_length_factor;
    box->rangez*=inv_length_factor;
}
