#include "ParticleSystem.h"
#include "hardCodedParameters.h"
#include <stdio.h>
#include <random>


/*
============================================================
free
============================================================
 */
void freeMemory(ParticleSystem* ps, int isGPUon){
    /* host*/
    free(ps->x);
    free(ps->v);
    free(ps->a);
    free(ps->f);
    free(ps->r);
    free(ps->rsq);
    free(ps->invr);


    free(ps->m);
    free(ps->sqrtm);
    free(ps->invm);
    free(ps->k);
    free(ps->etaconst);
    free(ps->g);

    free(ps->isActive);

    free(ps->angv);
    free(ps->anga);
    free(ps->moi);
    free(ps->invmoi);

    free(ps->mom);

    free(ps->deltHisx);
    free(ps->deltHisy);
    free(ps->deltHisz);

    free(ps->deltHisxWall);
    free(ps->deltHisyWall);
    free(ps->deltHiszWall);

    free(ps->indHis);
    free(ps->numCont);

    free(ps->isContact);
    free(ps->isContactWall);

    free(ps->indHisWall);
    free(ps->numContWall);

    free(ps->indHisWallNow);

    free(ps->indHisVorENow);

    free(ps->cellId);
    free(ps->cellx);

    free(ps->pId);
    free(ps->mortonKey);
    free(ps->mortonOrder);

    free(ps->tmpMortonKey);
    free(ps->tmpMortonOrder);

    /* ========== verlet list related ======= */
    free(ps->refx);
    free(ps->refy);
    free(ps->refz);
    free(ps->neiList);
    free(ps->numNei);
    free(ps->neiListWall);
    free(ps->numNeiWall);

    free(ps->walls.n);
    free(ps->walls.d);

    #if USE_GPU
    if (isGPUon ==1){
        cudaFree(ps->d_group.x);
        cudaFree(ps->d_group.v);
        cudaFree(ps->d_group.a);
        cudaFree(ps->d_group.f);
        cudaFree(ps->d_group.r);
        cudaFree(ps->d_group.rsq);
        cudaFree(ps->d_group.invr);


        cudaFree(ps->d_group.m);
        cudaFree(ps->d_group.sqrtm);
        cudaFree(ps->d_group.invm);
        cudaFree(ps->d_group.k);
        cudaFree(ps->d_group.etaconst);
        cudaFree(ps->d_group.g);

        cudaFree(ps->d_group.isActive);

        cudaFree(ps->d_group.angv);
        cudaFree(ps->d_group.anga);
        cudaFree(ps->d_group.moi);
        cudaFree(ps->d_group.invmoi);

        cudaFree(ps->d_group.mom);

        cudaFree(ps->d_group.deltHisx);
        cudaFree(ps->d_group.deltHisy);
        cudaFree(ps->d_group.deltHisz);

        cudaFree(ps->d_group.deltHisxWall);
        cudaFree(ps->d_group.deltHisyWall);
        cudaFree(ps->d_group.deltHiszWall);

        cudaFree(ps->d_group.indHis);
        cudaFree(ps->d_group.numCont);

        cudaFree(ps->d_group.isContact);
        cudaFree(ps->d_group.isContactWall);

        cudaFree(ps->d_group.indHisWall);
        cudaFree(ps->d_group.indHisWallNow);
        cudaFree(ps->d_group.indHisVorENow);

        cudaFree(ps->d_group.numContWall);


        cudaFree(ps->d_group.cellId);
        cudaFree(ps->d_group.cellx);

        cudaFree(ps->d_group.pId);
        cudaFree(ps->d_group.mortonKey);
        cudaFree(ps->d_group.mortonOrder);

        cudaFree(ps->d_group.tmpMortonKey);
        cudaFree(ps->d_group.tmpMortonOrder);

        cudaFree(ps->d_group.walls.n);
        cudaFree(ps->d_group.walls.d);

        cudaFree(ps->d_group.tmp_storage);


        /* ========== verlet list related ======= */
        cudaFree(ps->d_group.refx);
        cudaFree(ps->d_group.refy);
        cudaFree(ps->d_group.refz);

        cudaFree(ps->d_group.neiList);
        cudaFree(ps->d_group.numNei);

        cudaFree(ps->d_group.neiListWall);
        cudaFree(ps->d_group.numNeiWall);
        cudaFree(ps->d_group.refreshVerletFlag);

        /* ==== structs ====*/
        cudaFree(ps->d_groupPtr);
    }
    #endif
}


/* ==================== Memory Allocation =================== */
void allocate(ParticleSys<HostMemory> *ps){
    int N=ps->parameters.N;
    #define MEMBER(type,name,Np,SAVE_FLAG) ps->p.name = HostMemory::template allocate<type>(Np);
    #include "ParticleSystemMember_common.def"
    printf("Allocate Done\n");
    #undef MEMBER
}

void allocate(ParticleSys<DeviceMemory> *ps){
    int N=ps->parameters.N;
    printf("parameters N=%d\n",N);
    #define MEMBER(type,name,Np,SAVE_FLAG) ps->p.name = DeviceMemory::template allocate<type>((Np));
    #include "ParticleSystemMember_common.def"
    #undef MEMBER
    cudaMalloc((void**)&ps->deviceOnly.refreshVerletFlag,sizeof(int));
    
    /* ===== for sorting ===============*/
        ps->deviceOnly.tmp_storage=NULL;
        ps->deviceOnly.tmp_bytes=0;

        /* calculate temporarly size */

        cub::DeviceRadixSort::SortPairs(
                ps->deviceOnly.tmp_storage,
                ps->deviceOnly.tmp_bytes,
                ps->p.mortonKey,
                ps->p.tmpMortonKey,
                ps->p.mortonOrder,
                ps->p.tmpMortonOrder,
                N);

        cudaMalloc((void**)&ps->deviceOnly.tmp_storage,ps->deviceOnly.tmp_bytes);
}

void deallocate(ParticleSys<HostMemory> *ps){
    #define MEMBER(type,name,Np,SAVE_FLAG) HostMemory::deallocate(ps->p.name);
    #include "ParticleSystemMember_common.def"
    #undef MEMBER
}

void deallocate(ParticleSys<DeviceMemory> *ps){
    #define MEMBER(type,name,Np,SAVE_FLAG) DeviceMemory::deallocate(ps->p.name);
    #include "ParticleSystemMember_common.def"
    #undef MEMBER
    DeviceMemory::deallocate(ps->deviceOnly.refreshVerletFlag);
    DeviceMemory::deallocate(ps->deviceOnly.tmp_storage);
}

    


/*
   ============================================================
   Host→Device転送
   ============================================================
 */
void copyToDevice(ParticleSys<HostMemory> *ps, ParticleSys<DeviceMemory> *d_ps){
    size_t N = ps->parameters.N;

    #define MEMBER(type,name,Np,SAVE_FLAG) cudaMemcpy(d_ps->p.name,ps->p.name,Np*sizeof(type),cudaMemcpyHostToDevice); 
    #include "ParticleSystemMember_common.def"
    #undef MEMBER

    /* ========= parameters =========== */
    cudaMemcpy(&d_ps->parameters, &ps->parameters, sizeof(Parameters), cudaMemcpyHostToDevice);
    /* ========= commons  =========== */
    cudaMemcpy(&d_ps->p, &ps->p, sizeof(Common), cudaMemcpyHostToDevice);
}


/*
   ============================================================
   Device→Host転送
   ============================================================
 */
void copyFromDevice(ParticleSys<DeviceMemory>* d_ps, ParticleSys<HostMemory>* ps){
    size_t N = ps->parameters.N;

    #define MEMBER(type,name,Np,SAVE_FLAG) \
    if(SAVE_FLAG == SAVE_ON){ \
        cudaMemcpy(ps->p.name,d_ps->p.name, Np*sizeof(type), cudaMemcpyDeviceToHost);\
    }
    #include "ParticleSystemMember_common.def"
    #undef MEMBER
}

void copyFromDevice(ParticleSystem* ps)
{
    size_t size = ps->N * sizeof(double);

    cudaMemcpy(ps->x, ps->d_group.x, size*DIM, cudaMemcpyDeviceToHost);
    cudaMemcpy(ps->v, ps->d_group.v, size*DIM, cudaMemcpyDeviceToHost);
    cudaMemcpy(ps->a, ps->d_group.a, size*DIM, cudaMemcpyDeviceToHost);
    cudaMemcpy(ps->angv, ps->d_group.angv, size*DIM, cudaMemcpyDeviceToHost);
    cudaMemcpy(ps->r, ps->d_group.r, size, cudaMemcpyDeviceToHost);
    cudaMemcpy(ps->isActive, ps->d_group.isActive, ps->N*sizeof(int), cudaMemcpyDeviceToHost);
}

/*
   ============================================================
   メモリ確保
   ============================================================
 */
void allocateMemory(ParticleSystem* ps, int isGPUon){
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


    ps->isActive = (int*)malloc(sizeof(int)*ps->N);

    ps->angv = (double*)malloc(size*DIM);
    ps->anga = (double*)malloc(size*DIM);
    ps->moi = (double*)malloc(size);
    ps->invmoi = (double*)malloc(size);

    ps->mom = (double*)malloc(size*DIM);

    ps->deltHisx = (double*)malloc(sizeof(double)*ps->N*MAX_NEI);
    ps->deltHisy = (double*)malloc(sizeof(double)*ps->N*MAX_NEI);
    ps->deltHisz = (double*)malloc(sizeof(double)*ps->N*MAX_NEI);

    ps->deltHisxWall = (double*)malloc(sizeof(double)*ps->N*MAX_NEI);
    ps->deltHisyWall = (double*)malloc(sizeof(double)*ps->N*MAX_NEI);
    ps->deltHiszWall = (double*)malloc(sizeof(double)*ps->N*MAX_NEI);

    ps->indHis = (int*)malloc(sizeof(int)*ps->N*MAX_NEI);
    ps->numCont = (int*)malloc(sizeof(int)*ps->N);

    ps->isContact = (int*)malloc(sizeof(int)*ps->N*MAX_NEI);
    ps->isContactWall = (int*)malloc(sizeof(int)*ps->N*MAX_NEI);

    ps->indHisWall = (int*)malloc(sizeof(int)*ps->N*MAX_NEI);
    ps->numContWall = (int*)malloc(sizeof(int)*ps->N);

    ps->indHisWallNow = (int*)malloc(sizeof(int)*ps->N*MAX_NEI);

    ps->indHisVorENow = (int*)malloc(sizeof(int)*ps->N*MAX_NEI);

    ps->cellId = (int*)malloc(sizeof(int)*ps->N);
    ps->cellx = (int*)malloc(sizeof(int)*DIM*ps->N);

    ps->pId = (int*)malloc(sizeof(int)*ps->N);
    ps->mortonKey = (uint32_t*)malloc(sizeof(uint32_t)*ps->N);
    ps->mortonOrder = (int*)malloc(sizeof(int)*ps->N);
    ps->tmpMortonKey = (uint32_t*)malloc(sizeof(uint32_t)*ps->N);
    ps->tmpMortonOrder = (int*)malloc(sizeof(int)*ps->N);

    /* ========== verlet list related ======= */
    ps->refx = (double*)malloc(size);
    ps->refy = (double*)malloc(size);
    ps->refz = (double*)malloc(size);
    ps->neiList = (int*)malloc(sizeof(int)*ps->N*MAX_NEI);
    ps->numNei = (int*)malloc(sizeof(int)*ps->N);
    ps->neiListWall = (int*)malloc(sizeof(int)*ps->N*MAX_NEI);
    ps->numNeiWall = (int*)malloc(sizeof(int)*ps->N);

    ps->walls.n = (double*)malloc(size_walls*DIM);
    ps->walls.d = (double*)malloc(size_walls);

    // Device
    #if USE_GPU
    if (isGPUon == 1){
        cudaMalloc((void**)&ps->d_group.x,  size*DIM);
        cudaMalloc((void**)&ps->d_group.v,  size*DIM);
        cudaMalloc((void**)&ps->d_group.a,  size*DIM);
        cudaMalloc((void**)&ps->d_group.f,  size*DIM);
        cudaMalloc((void**)&ps->d_group.r,  size);
        cudaMalloc((void**)&ps->d_group.rsq,  size);
        cudaMalloc((void**)&ps->d_group.invr, size);

        cudaMalloc((void**)&ps->d_group.m, size);
        cudaMalloc((void**)&ps->d_group.sqrtm, size);
        cudaMalloc((void**)&ps->d_group.invm, size);
        cudaMalloc((void**)&ps->d_group.k, size);
        cudaMalloc((void**)&ps->d_group.etaconst, size);
        cudaMalloc((void**)&ps->d_group.g, sizeof(double)*DIM);

        cudaMalloc((void**)&ps->d_group.isActive, sizeof(int)*ps->N);

        cudaMalloc((void**)&ps->d_group.angv, size*DIM);
        cudaMalloc((void**)&ps->d_group.anga, size*DIM);
        cudaMalloc((void**)&ps->d_group.moi, size);
        cudaMalloc((void**)&ps->d_group.invmoi, size);

        cudaMalloc((void**)&ps->d_group.mom, size*DIM);

        cudaMalloc((void**)&ps->d_group.deltHisx, sizeof(double)*ps->N*MAX_NEI);
        cudaMalloc((void**)&ps->d_group.deltHisy, sizeof(double)*ps->N*MAX_NEI);
        cudaMalloc((void**)&ps->d_group.deltHisz, sizeof(double)*ps->N*MAX_NEI);

        cudaMalloc((void**)&ps->d_group.deltHisxWall, sizeof(double)*ps->N*MAX_NEI);
        cudaMalloc((void**)&ps->d_group.deltHisyWall, sizeof(double)*ps->N*MAX_NEI);
        cudaMalloc((void**)&ps->d_group.deltHiszWall, sizeof(double)*ps->N*MAX_NEI);

        cudaMalloc((void**)&ps->d_group.indHis, sizeof(int)*ps->N*MAX_NEI);
        cudaMalloc((void**)&ps->d_group.numCont, sizeof(int)*ps->N);

        cudaMalloc((void**)&ps->d_group.isContact, sizeof(int)*ps->N*MAX_NEI);
        cudaMalloc((void**)&ps->d_group.isContactWall, sizeof(int)*ps->N*MAX_NEI);

        cudaMalloc((void**)&ps->d_group.indHisWall, sizeof(int)*ps->N*MAX_NEI);
        cudaMalloc((void**)&ps->d_group.indHisWallNow, sizeof(int)*ps->N*MAX_NEI);
        cudaMalloc((void**)&ps->d_group.indHisVorENow, sizeof(int)*ps->N*MAX_NEI);

        cudaMalloc((void**)&ps->d_group.numContWall, sizeof(int)*ps->N);

        cudaMalloc((void**)&ps->d_group.cellId, sizeof(int)*ps->N);
        cudaMalloc((void**)&ps->d_group.cellx, sizeof(int)*DIM*ps->N);

        cudaMalloc((void**)&ps->d_group.pId, sizeof(int)*ps->N);
        cudaMalloc((void**)&ps->d_group.mortonKey, sizeof(uint32_t)*ps->N);
        cudaMalloc((void**)&ps->d_group.mortonOrder, sizeof(int)*ps->N);
        cudaMalloc((void**)&ps->d_group.tmpMortonKey, sizeof(uint32_t)*ps->N);

        /* ========== verlet list related ======= */
        cudaMalloc((void**)&ps->d_group.refx, size);
        cudaMalloc((void**)&ps->d_group.refy, size);
        cudaMalloc((void**)&ps->d_group.refz, size);

        cudaMalloc((void**)&ps->d_group.neiList, sizeof(int) * ps->N * MAX_NEI);
        cudaMalloc((void**)&ps->d_group.numNei, sizeof(int) * ps->N);

        cudaMalloc((void**)&ps->d_group.neiListWall, sizeof(int) * ps->N * MAX_NEI);
        cudaMalloc((void**)&ps->d_group.numNeiWall, sizeof(int) * ps->N);

        cudaMalloc((void**)&ps->d_group.tmpMortonOrder, sizeof(int)*ps->N);
        cudaMalloc((void**)&ps->d_group.walls.n, size_walls*DIM);
        cudaMalloc((void**)&ps->d_group.walls.d, size_walls);
        cudaMalloc((void**)&ps->d_group.refreshVerletFlag, sizeof(int));

        /* ======== for sorting ==============*/

        ps->d_group.tmp_storage=NULL;
        ps->d_group.tmp_bytes=0;

        /* calculate temporarly size */

        cub::DeviceRadixSort::SortPairs(
                ps->d_group.tmp_storage,
                ps->d_group.tmp_bytes,
                ps->d_group.mortonKey,
                ps->d_group.tmpMortonKey,
                ps->d_group.mortonOrder,
                ps->d_group.tmpMortonOrder,
                ps->N);

        cudaMalloc((void**)&ps->d_group.tmp_storage,ps->d_group.tmp_bytes);

        /* ======== for structs ============*/
        cudaMalloc((void**)&ps->d_groupPtr,sizeof(DeviceParticleGroup));

    }

    #endif

}

/*
   ============================================================
   初期化
   ============================================================
 */
void initializeTmpParticles(ParticleSys<HostMemory>* ps,cJSON *json_inlet, double r,double m,double k,double res)
{
    for (int i = 0; i < ps->parameters.N; i++){
        ps->p.x[i*DIM+0] = 0.;
        ps->p.x[i*DIM+1] = 0.;
        ps->p.x[i*DIM+2] = 0.;


        /*
           ps->p.x[i*DIM+0] = 0.25;
           ps->p.x[i*DIM+1] = (double)i*1.0+0.03;
           ps->p.x[i*DIM+2] = 0.25;
         */

        ps->p.r[i] = r;
        ps->p.rsq[i] = ps->p.r[i]*ps->p.r[i];
        ps->p.invr[i] = 1./ps->p.r[i];
        ps->p.k[i] = k;
        ps->p.m[i] = m;
        ps->p.sqrtm[i] = sqrt(m);
        ps->p.invm[i] = 1./m;
        ps->p.etaconst[i]=-2.*log(res)*sqrt(ps->p.k[i]/(3.1415*3.1415+log(res)*log(res)));
        ps->p.cellId[i] = -1;
        ps->p.isActive[i] = 1;


        ps->p.numCont[i] = 0;
        ps->p.numContWall[i] = 0;

        ps->p.v[i*DIM+0] = 0.;
        ps->p.v[i*DIM+1] = 0.;
        ps->p.v[i*DIM+2] = 0.;

        /* ======== for temporarly check========= */
        /*
           if (i==1){
           ps->p.v[i*DIM+0] = 1.;
           ps->p.v[i*DIM+1] = 0.;
           ps->p.v[i*DIM+2] = 0.;
           }
         */
        /* ======== for temporarly check========= */

        ps->p.angv[i*DIM+0] = 0.;
        ps->p.angv[i*DIM+1] = 0.;
        ps->p.angv[i*DIM+2] = 0.;

        ps->p.anga[i*DIM+0] = 0.;
        ps->p.anga[i*DIM+1] = 0.;
        ps->p.anga[i*DIM+2] = 0.;

        ps->p.moi[i] = 2./5. * ps->p.m[i]*ps->p.rsq[i];
        ps->p.invmoi[i] = 1./ps->p.moi[i];

        ps->p.pId[i] = i;


        for (int j=0; j<MAX_NEI; j++){
            ps->p.indHis[i*MAX_NEI+j] = -1;
            ps->p.indHisWall[i*MAX_NEI+j] = -1;
            ps->p.indHisWallNow[i*MAX_NEI+j] = -1;

            ps->p.isContact[i*MAX_NEI+j] = -1;
            ps->p.isContactWall[i*MAX_NEI+j] = -1;

            ps->p.deltHisx[i*MAX_NEI+j] = 0.;
            ps->p.deltHisy[i*MAX_NEI+j] = 0.;
            ps->p.deltHisz[i*MAX_NEI+j] = 0.;
            ps->p.deltHisxWall[i*MAX_NEI+j] = 0.;
            ps->p.deltHisyWall[i*MAX_NEI+j] = 0.;
            ps->p.deltHiszWall[i*MAX_NEI+j] = 0.;
        }
    }
}

void initializeParticles(ParticleSys<HostMemory>* ps,cJSON *json_inlet, double r,double m,double k,double res)
{
    int max_trials = 1000; // 1粒子あたりの再配置試行回数
    int seed = 1;
    std::mt19937 mt(seed);
    std::uniform_real_distribution<double> uni(0,1);


    if(strcmp(cJSON_GetObjectItem(json_inlet,"inputMode")->valuestring,"shuffle")==0){

        printf("shuffling particles\n");

        /* == get box for shuffling == */

        cJSON *json_shuffle = cJSON_GetObjectItem(json_inlet,"shuffle");
        double minx_shuffle = cJSON_GetObjectItem(json_shuffle,"minx")->valuedouble;
        double miny_shuffle = cJSON_GetObjectItem(json_shuffle,"miny")->valuedouble;
        double minz_shuffle = cJSON_GetObjectItem(json_shuffle,"minz")->valuedouble;
        double maxx_shuffle = cJSON_GetObjectItem(json_shuffle,"maxx")->valuedouble;
        double maxy_shuffle = cJSON_GetObjectItem(json_shuffle,"maxy")->valuedouble;
        double maxz_shuffle = cJSON_GetObjectItem(json_shuffle,"maxz")->valuedouble;

        for (int i = 0; i < ps->parameters.N; i++){
            int trial = 0;
            while (trial < max_trials)
            {
                // ランダム配置
                /*
                   double x = (double)rand() / RAND_MAX*(0.3)-(0.35+0.15);
                   double y = (double)rand() / RAND_MAX * 1.0 + 2.8;
                //double y = -0.97;
                double z = (double)rand() / RAND_MAX*(0.3)-(0.35+0.15);
                 */

                double x = uni(mt)*(maxx_shuffle-minx_shuffle) + minx_shuffle;
                double y = uni(mt)*(maxy_shuffle-miny_shuffle) + miny_shuffle;
                double z = uni(mt)*(maxz_shuffle-minz_shuffle) + minz_shuffle;


                /* ======== for temporarly check========= */
                /*
                   double x = 0.;
                   double y = 0.021;
                   double z = 0.0;

                   if (i==0){
                   x = 0.;
                   y = 0.0;
                   z = 0.0;
                   }else{
                   x = -0.2;
                   y = 0.005;
                   z = 0.0;

                   }
                 */
                /* ======== for temporarly check========= */


                // 既存粒子との距離チェック
                int overlap = 0;
                for (int j = 0; j < i; j++)
                {
                    double dx = x - ps->p.x[j*DIM+0];
                    double dy = y - ps->p.x[j*DIM+1];
                    double dz = z - ps->p.x[j*DIM+2];
                    double d2 = dx*dx + dy*dy + dz*dz;
                    if (d2 < 4*r*r) // 半径2倍以上離れているか
                    {
                        overlap = 1;
                        break;
                    }
                }

                if (!overlap)
                {
                    ps->p.x[i*DIM+0] = x;
                    ps->p.x[i*DIM+1] = y;
                    ps->p.x[i*DIM+2] = z;
                    break;
                }

                trial++;
            }

            if (trial == max_trials){
                printf("!!!couldn't place them!!\n");
                abort();
            };


            /*
               ps->p.x[i*DIM+0] = 0.25;
               ps->p.x[i*DIM+1] = (double)i*1.0+0.03;
               ps->p.x[i*DIM+2] = 0.25;
             */

            ps->p.r[i] = r;
            ps->p.rsq[i] = ps->p.r[i]*ps->p.r[i];
            ps->p.invr[i] = 1./ps->p.r[i];
            ps->p.k[i] = k;
            ps->p.m[i] = m;
            ps->p.sqrtm[i] = sqrt(m);
            ps->p.invm[i] = 1./m;
            ps->p.etaconst[i]=-2.*log(res)*sqrt(ps->p.k[i]/(3.1415*3.1415+log(res)*log(res)));
            ps->p.cellId[i] = -1;
            ps->p.isActive[i] = 1;


            ps->p.numCont[i] = 0;
            ps->p.numContWall[i] = 0;

            ps->p.v[i*DIM+0] = 0.;
            ps->p.v[i*DIM+1] = 0.;
            ps->p.v[i*DIM+2] = 0.;

            /* ======== for temporarly check========= */
            /*
               if (i==1){
               ps->p.v[i*DIM+0] = 1.;
               ps->p.v[i*DIM+1] = 0.;
               ps->p.v[i*DIM+2] = 0.;
               }
             */
            /* ======== for temporarly check========= */

            ps->p.angv[i*DIM+0] = 0.;
            ps->p.angv[i*DIM+1] = 0.;
            ps->p.angv[i*DIM+2] = 0.;

            ps->p.anga[i*DIM+0] = 0.;
            ps->p.anga[i*DIM+1] = 0.;
            ps->p.anga[i*DIM+2] = 0.;

            ps->p.moi[i] = 2./5. * ps->p.m[i]*ps->p.rsq[i];
            ps->p.invmoi[i] = 1./ps->p.moi[i];

            ps->p.pId[i] = i;


            for (int j=0; j<MAX_NEI; j++){
                ps->p.indHis[i*MAX_NEI+j] = -1;
                ps->p.indHisWall[i*MAX_NEI+j] = -1;
                ps->p.indHisWallNow[i*MAX_NEI+j] = -1;

                ps->p.isContact[i*MAX_NEI+j] = -1;
                ps->p.isContactWall[i*MAX_NEI+j] = -1;

                ps->p.deltHisx[i*MAX_NEI+j] = 0.;
                ps->p.deltHisy[i*MAX_NEI+j] = 0.;
                ps->p.deltHisz[i*MAX_NEI+j] = 0.;
                ps->p.deltHisxWall[i*MAX_NEI+j] = 0.;
                ps->p.deltHisyWall[i*MAX_NEI+j] = 0.;
                ps->p.deltHiszWall[i*MAX_NEI+j] = 0.;
            }
        }
    }else if(strcmp(cJSON_GetObjectItem(json_inlet,"inputMode")->valuestring,"file")==0){
        printf("to be created...");
        abort();

    }else{
        printf("no proper input mode defined!\n");
        abort();
    }

}

void nondimensionalize(ParticleSys<HostMemory>* ps, BoundingBox *box, TriangleMesh* mesh){
    /* get nondimensionalize factor */
    ps->parameters.time_factor = sqrt(ps->p.m[0]/ps->p.k[0]);
    ps->parameters.mass_factor = ps->p.m[0];
    ps->parameters.length_factor = ps->p.r[0];
    double time_factor = ps->parameters.time_factor;
    double mass_factor = ps->parameters.mass_factor;
    double length_factor = ps->parameters.length_factor;

    double inv_time_factor = 1./time_factor;
    double inv_mass_factor = 1./mass_factor;
    double inv_length_factor = 1./length_factor;
    double inv_sqrtk_factor = 1./sqrt(ps->p.k[0]);

    ps->parameters.dt*=inv_time_factor;

    ps->p.g[0]*=inv_length_factor*time_factor*time_factor;
    ps->p.g[1]*=inv_length_factor*time_factor*time_factor;
    ps->p.g[2]*=inv_length_factor*time_factor*time_factor;



    for (int i=0; i<ps->parameters.N; i++){
        ps->p.x[i*DIM+0]*=inv_length_factor;
        ps->p.x[i*DIM+1]*=inv_length_factor;
        ps->p.x[i*DIM+2]*=inv_length_factor;

        ps->p.v[i*DIM+0]*=inv_length_factor*time_factor;
        ps->p.v[i*DIM+1]*=inv_length_factor*time_factor;
        ps->p.v[i*DIM+2]*=inv_length_factor*time_factor;

        ps->p.r[i]*=inv_length_factor;
        ps->p.rsq[i]*=inv_length_factor*inv_length_factor;
        ps->p.invr[i]*=length_factor;

        ps->p.m[i]*=inv_mass_factor;
        ps->p.sqrtm[i]*=sqrt(inv_mass_factor);

        ps->p.invm[i]*=mass_factor;

        ps->p.moi[i]*=inv_mass_factor*inv_length_factor*inv_length_factor;
        ps->p.invmoi[i]*=mass_factor*length_factor*length_factor;

        ps->p.k[i]*=inv_sqrtk_factor*inv_sqrtk_factor;
        ps->p.etaconst[i]*=inv_sqrtk_factor;
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

    box->invdx*=length_factor;
    box->invdy*=length_factor;
    box->invdz*=length_factor;

    box->rangex*=inv_length_factor;
    box->rangey*=inv_length_factor;
    box->rangez*=inv_length_factor;
    box->skinR*=inv_length_factor;
    box->refreshThresh*=inv_length_factor;
    box->refreshThreshSq*=inv_length_factor*inv_length_factor;

    /* triangles */
    for (int i=0; i<mesh->nTri; i++){
        mesh->cx[i]*= inv_length_factor;
        mesh->cy[i]*= inv_length_factor;
        mesh->cz[i]*= inv_length_factor;

        mesh->d[i]*= inv_length_factor;

        mesh->e01x[i]*= inv_length_factor;
        mesh->e01y[i]*= inv_length_factor;
        mesh->e01z[i]*= inv_length_factor;

        mesh->e02x[i]*= inv_length_factor;
        mesh->e02y[i]*= inv_length_factor;
        mesh->e02z[i]*= inv_length_factor;

        mesh->e12x[i]*= inv_length_factor;
        mesh->e12y[i]*= inv_length_factor;
        mesh->e12z[i]*= inv_length_factor;

        mesh->d00[i]*= inv_length_factor*inv_length_factor;
        mesh->d00inv[i]*= length_factor*length_factor;
        mesh->d01[i]*= inv_length_factor*inv_length_factor;
        mesh->d11[i]*= inv_length_factor*inv_length_factor;
        mesh->d11inv[i]*=length_factor*length_factor;
        mesh->d22[i]*= inv_length_factor*inv_length_factor;
        mesh->d22inv[i]*=length_factor*length_factor;

        mesh->denom[i]*= length_factor*length_factor*length_factor*length_factor;

        mesh->minx[i]*= inv_length_factor;
        mesh->miny[i]*= inv_length_factor;
        mesh->minz[i]*= inv_length_factor;
        mesh->maxx[i]*= inv_length_factor;
        mesh->maxy[i]*= inv_length_factor;
        mesh->maxz[i]*= inv_length_factor;
    }
    for (int i=0; i<mesh->nVert; i++){
        mesh->mx[i]*= inv_length_factor;
        mesh->my[i]*= inv_length_factor;
        mesh->mz[i]*= inv_length_factor;
    }

}
