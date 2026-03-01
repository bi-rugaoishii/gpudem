#include "BoundingBox.h"
#include "ParticleSystem.h"

/* =========== used for sort =========== */


int compare_int(const void *a, const void *b){
    int ia = *(const int*)a;
    int ib = *(const int*)b;

    if (ia < ib) return -1;
    if (ia > ib) return  1;
    return 0;
}

/* ============ Morton key related ================ */

static inline uint32_t expandBits(uint32_t v){
    v &= 0x000003ff;                 // 10bit
    v = (v | (v << 16)) & 0x030000FF;
    v = (v | (v <<  8)) & 0x0300F00F;
    v = (v | (v <<  4)) & 0x030C30C3;
    v = (v | (v <<  2)) & 0x09249249;
    return v;
}
static inline uint32_t morton3D(uint32_t ix, uint32_t iy,uint32_t iz){
    return (expandBits(ix) << 2)| (expandBits(iy) << 1)| (expandBits(iz));
}

// ======================================================
// ポインタスワップ版 radix sort（最速設計）
// keyPtr/indexPtr は「ポインタへのポインタ」
// ======================================================
void radixSortUint32(
    uint32_t** keyPtr,
    int**      indexPtr,
    int N,
    uint32_t*  workKey,
    int*       workIndex)
{
    const int RADIX = 256;
    const int PASS  = 4;   // 30bit Mortonなら3passでOK

    uint32_t* srcKey = *keyPtr;
    uint32_t* dstKey = workKey;

    int* srcIdx = *indexPtr;
    int* dstIdx = workIndex;

    int pass;

    for(pass=0; pass<PASS; pass++)
    {
        int count[RADIX];
        int i;

        memset(count,0,sizeof(count));

        int shift = pass*8;

        // --- histogram ---
        for(i=0;i<N;i++)
        {
            uint32_t digit = (srcKey[i] >> shift) & 0xFF;
            count[digit]++;
        }

        // --- prefix sum ---
        int sum = 0;
        for(i=0;i<RADIX;i++)
        {
            int c = count[i];
            count[i] = sum;
            sum += c;
        }

        // --- scatter ---
        for(i=0;i<N;i++)
        {
            uint32_t k = srcKey[i];
            uint32_t digit = (k >> shift) & 0xFF;

            int pos = count[digit]++;

            dstKey[pos] = k;
            dstIdx[pos] = srcIdx[i];
        }

        // --- swap pointer only ---
        {
            uint32_t* tk = srcKey;
            srcKey = dstKey;
            dstKey = tk;
        }
        {
            int* ti = srcIdx;
            srcIdx = dstIdx;
            dstIdx = ti;
        }
    }

    // 最終ポインタを返す（コピー無し）
    *keyPtr   = srcKey;
    *indexPtr = srcIdx;
}

void swap_ps(ParticleSystem *p, ParticleSystem *tmp){
    int N = p->N;


    for (int i=0; i<N; i++){
        int src = p->mortonOrder[i];

        // --- position ---
        tmp->x[i*3+0] = p->x[src*3+0];
        tmp->x[i*3+1] = p->x[src*3+1];
        tmp->x[i*3+2] = p->x[src*3+2];

        // --- velocity ---
        tmp->v[i*3+0] = p->v[src*3+0];
        tmp->v[i*3+1] = p->v[src*3+1];
        tmp->v[i*3+2] = p->v[src*3+2];

        // --- ang velocity ---
        tmp->angv[i*3+0] = p->angv[src*3+0];
        tmp->angv[i*3+1] = p->angv[src*3+1];
        tmp->angv[i*3+2] = p->angv[src*3+2];

        // --- scalar ---
        tmp->r[i] = p->r[src];
        tmp->rsq[i] = p->rsq[src];
        tmp->invr[i] = p->invr[src];

        tmp->m[i] = p->m[src];
        tmp->sqrtm[i] = p->sqrtm[src];
        tmp->invm[i] = p->invm[src];
        tmp->moi[i] = p->moi[src];
        tmp->invmoi[i] = p->invmoi[src];
        tmp->etaconst[i] = p->etaconst[src];

        tmp->isActive[i] = p->isActive[src];

        int bi = i*p->MAX_NEI;
        int bsrc = src*p->MAX_NEI;
        for (int j=0; j<p->MAX_NEI; j++){
            // ---- contact history (particle) ----
            tmp->deltHisx[bi+j] = p->deltHisx[bsrc+j];
            tmp->deltHisy[bi+j] = p->deltHisy[bsrc+j];
            tmp->deltHisz[bi+j] = p->deltHisz[bsrc+j];

            tmp->isContact[bi+j] = p->isContact[bsrc+j];

            /* get the after order index of contact particle */
            tmp->indHis[bi+j] = p->indHis[bsrc+j];

            // ---- contact history (wall) ----
            tmp->deltHisxWall[bi+j] = p->deltHisxWall[bsrc+j];
            tmp->deltHisyWall[bi+j] = p->deltHisyWall[bsrc+j];
            tmp->deltHiszWall[bi+j] = p->deltHiszWall[bsrc+j];
            tmp->isContactWall[bi+j] = p->isContactWall[bsrc+j];
            tmp->indHisWall[bi+j] = p->indHisWall[bsrc+j];

        }

        // ---- number of contact ----
        tmp->numCont[i]     = p->numCont[src];
        tmp->numContWall[i] = p->numContWall[src];

        // ---- cell info ----
        tmp->cellId[i] = p->cellId[src];
        tmp->cellx[i*DIM+0]  = p->cellx[src*DIM+0];
        tmp->cellx[i*DIM+1]  = p->cellx[src*DIM+1];
        tmp->cellx[i*DIM+2]  = p->cellx[src*DIM+2];


        // ---- morton key ----
        tmp->mortonKey[i] = p->mortonKey[src];
        tmp->pId[i] = p->pId[src];
    }


    /* swap pointers*/
    double* td;
    int*    ti;
    uint32_t* tu;

    td=p->x; p->x=tmp->x; tmp->x=td;
    td=p->v; p->v=tmp->v; tmp->v=td;
    td=p->r; p->r=tmp->r; tmp->r=td;
    td=p->rsq; p->rsq=tmp->rsq; tmp->rsq=td;
    td=p->invr; p->invr=tmp->invr; tmp->invr=td;

    td=p->m; p->m=tmp->m; tmp->m=td;
    td=p->invm; p->invm=tmp->invm; tmp->invm=td;
    td=p->sqrtm; p->sqrtm=tmp->sqrtm; tmp->sqrtm=td;
    td=p->moi; p->moi=tmp->moi; tmp->moi=td;
    td=p->invmoi; p->invmoi=tmp->invmoi; tmp->invmoi=td;
    td=p->etaconst; p->etaconst=tmp->etaconst; tmp->etaconst=td;

    td=p->angv; p->angv=tmp->angv; tmp->angv=td;
    // ---- contact history ----
    td=p->deltHisx;     p->deltHisx=tmp->deltHisx;     tmp->deltHisx=td;
    td=p->deltHisy;     p->deltHisy=tmp->deltHisy;     tmp->deltHisy=td;
    td=p->deltHisz;     p->deltHisz=tmp->deltHisz;     tmp->deltHisz=td;

    td=p->deltHisxWall; p->deltHisxWall=tmp->deltHisxWall; tmp->deltHisxWall=td;
    td=p->deltHisyWall; p->deltHisyWall=tmp->deltHisyWall; tmp->deltHisyWall=td;
    td=p->deltHiszWall; p->deltHiszWall=tmp->deltHiszWall; tmp->deltHiszWall=td;

    // ---- contact number ----
    ti=p->numCont;      p->numCont=tmp->numCont;      tmp->numCont=ti;
    ti=p->numContWall;  p->numContWall=tmp->numContWall; tmp->numContWall=ti;

    // ---- cell info ----
    ti=p->cellId;       p->cellId=tmp->cellId;       tmp->cellId=ti;
    ti=p->cellx;        p->cellx=tmp->cellx;         tmp->cellx=ti;

    // ---- particle id ----
    ti=p->pId;          p->pId=tmp->pId;             tmp->pId=ti;

    ti=p->isActive;       p->isActive=tmp->isActive;       tmp->isActive=ti;
    // ---- morton key ----
    tu=p->mortonKey;    p->mortonKey=tmp->mortonKey; tmp->mortonKey=tu;

}

void copyToDeviceBox(BoundingBox *box, ParticleSystem *ps){

    size_t size = box->N;

    box->d_box.minx = box->minx;
    box->d_box.miny = box->miny;
    box->d_box.minz = box->minz;
    box->d_box.maxx = box->maxx;
    box->d_box.maxy = box->maxy;
    box->d_box.maxz = box->maxz;

    box->d_box.dx = box->dx;
    box->d_box.dy = box->dy;
    box->d_box.dz = box->dz;

    box->d_box.invdx = box->invdx;
    box->d_box.invdy = box->invdy;
    box->d_box.invdz = box->invdz;

    box->d_box.rangex = box->rangex;
    box->d_box.rangey = box->rangey;
    box->d_box.rangez = box->rangez;

    box->d_box.sizex = box->sizex;
    box->d_box.sizey = box->sizey;
    box->d_box.sizez = box->sizez;

    box->d_box.N = box->N;

    cudaMemcpy(box->d_box.pList,  box->pList,  sizeof(int)*ps->N, cudaMemcpyHostToDevice);
    cudaMemcpy(box->d_box.pNum,  box->pNum,  sizeof(int)*size, cudaMemcpyHostToDevice);
    cudaMemcpy(box->d_box.pStart,  box->pStart,  sizeof(int)*size, cudaMemcpyHostToDevice);
    cudaMemcpy(box->d_box.cellOffset,  box->cellOffset,  sizeof(int)*size, cudaMemcpyHostToDevice);
    cudaMemcpy(box->d_box.tList,  box->tList,  sizeof(int)*size*box->MAX_TRI, cudaMemcpyHostToDevice);
    cudaMemcpy(box->d_box.tNum,  box->tNum,  sizeof(int)*size, cudaMemcpyHostToDevice);
    cudaMemcpy(box->d_boxPtr,  &box->d_box,  sizeof(DeviceBoundingBox), cudaMemcpyHostToDevice);
}

void free_BoundingBox(BoundingBox *box){

    free(box->pList);
    free(box->pNum);
    free(box->pStart);
    free(box->cellOffset);
    free(box->usedCells);

    #if USE_GPU
    cudaFree(box->d_box.pList);
    cudaFree(box->d_box.tmpExSum);
    cudaFree(box->d_box.pNum);
    cudaFree(box->d_box.pStart);
    cudaFree(box->d_box.cellOffset);
    cudaFree(box->d_box.tList);
    cudaFree(box->d_box.tNum);
    cudaFree(box->d_boxPtr);
    #endif
}

void initialize_BoundingBox(ParticleSystem *p, BoundingBox *box,TriangleMesh *mesh, double minx, double maxx, double miny, double maxy, double minz, double maxz){
    /* find maximum radius */

    double maxr=0.;
    for (int i=0; i<p->N; i++){
        if (maxr<p->r[i]){
            maxr=p->r[i];
        }
    }

    box->minx = minx;
    box->miny = miny;
    box->minz = minz;
    box->maxx = maxx;
    box->maxy = maxy;
    box->maxz = maxz;
    box->rangex = maxx - minx;
    box->rangey = maxy - miny;
    box->rangez = maxz - minz;

    box->dx = maxr*4. ;
    box->dy = maxr*4. ;
    box->dz = maxr*4. ;

    box->invdx = 1./box->dx;
    box->invdy = 1./box->dy;
    box->invdz = 1./box->dz;

    box->sizex = ceil(box->rangex/box->dx) +2;
    box->sizey = ceil(box->rangey/box->dy) +2;
    box->sizez = ceil(box->rangez/box->dz) +2;

    int sizeBox= box->sizex*box->sizey*box->sizez;
    box->N= sizeBox;

    box->pList = (int*)malloc(sizeof(int)*p->N);
    box->pNum = (int*)calloc(sizeBox,sizeof(int));
    box->pStart = (int*)calloc(sizeBox,sizeof(int));
    box->cellOffset = (int*)calloc(sizeBox,sizeof(int));

    box->usedCells = (int*)malloc(sizeof(int)*sizeBox);
    box->numUsedCells = 0;

    /* ==== for triangles ==== */

    box->MAX_TRI = 100;
    box->tList = (int*)calloc(sizeBox*box->MAX_TRI,sizeof(int));
    box->tNum = (int*)calloc(sizeBox,sizeof(int));

    #if USE_GPU
    box->d_box.N= sizeBox;
    box->d_box.MAX_TRI = box->MAX_TRI;
    cudaMalloc(&box->d_box.pList, sizeof(int)*p->N);
    cudaMalloc(&box->d_box.pNum, sizeof(int)*sizeBox);
    cudaMalloc(&box->d_box.pStart, sizeof(int)*sizeBox);
    cudaMalloc(&box->d_box.tList, sizeof(int)*sizeBox*box->MAX_TRI);
    cudaMalloc(&box->d_box.tNum, sizeof(int)*sizeBox);
    cudaMalloc(&box->d_box.cellOffset, sizeof(int)*sizeBox);

    /* for cub */

    box->d_box.tmpExSum =NULL;
    box->d_box.scanTmpBytes = 0;

    cub::DeviceScan::ExclusiveSum(
            box->d_box.tmpExSum,
            box->d_box.scanTmpBytes,
            box->d_box.pNum,
            box->d_box.pStart,
            box->d_box.N);
    cudaMalloc(&box->d_box.tmpExSum, box->d_box.scanTmpBytes);

    /* for structs*/
    cudaMalloc(&box->d_boxPtr,sizeof(DeviceBoundingBox));
    #endif

}

void update_pList_withSort_fast(ParticleSystem *p, ParticleSystem *tmpPs,BoundingBox *box){
    /* initialize */
    for (int i=0; i<box->numUsedCells; i++){
        int cid = box->usedCells[i];
        box->pNum[cid] = 0;
        box->pStart[cid] = 0;
        box->cellOffset[cid] = 0;
    }

    for(int i=0;i<p->N;i++){
        p->mortonOrder[i] = i;
    }
    /* get cellId and get mortonkey*/
    for (int i=0; i<p->N; i++){
        int dx = floor((p->x[i*DIM+0]-box->minx)*box->invdx)+1; //+1 for ghost cell
        int dy = floor((p->x[i*DIM+1]-box->miny)*box->invdy)+1; //+1 for ghost cell
        int dz = floor((p->x[i*DIM+2]-box->minz)*box->invdz)+1; //+1 for ghost cell
        p->cellx[i*DIM+0] = dx;
        p->cellx[i*DIM+1] = dy;
        p->cellx[i*DIM+2] = dz;
        p->cellId[i] = (box->sizey*dz+dy)*box->sizex+dx;
        p->mortonKey[i]=morton3D((uint32_t)dx,(uint32_t)dy,(uint32_t)dz);
    }

    /* ========= sort and reorder ========= */
    radixSortUint32(&p->mortonKey,&p->mortonOrder,p->N,p->tmpMortonKey,p->tmpMortonOrder);

    swap_ps(p, tmpPs);



    /* count number of particle in cell */

    box->numUsedCells = 0;

    for (int i=0; i<p->N; i++){
        if (p->isActive[i]!=1){
            continue;
        }
        int cellId = p->cellId[i];

        if (box->pNum[cellId]==0){
            box->usedCells[box->numUsedCells] = cellId;
            box->numUsedCells += 1;
        }

        box->pNum[cellId] +=1;
    }


    int runningSum = 0;
    for (int i=0; i<box->numUsedCells; i++){
        int cid = box->usedCells[i];
        box->pStart[cid]=runningSum;
        runningSum += box->pNum[cid];
    }


    for (int i=0; i<p->N; i++){
        if (p->isActive[i]!=1){
            continue;
        }

        int cellId = p->cellId[i];
        int startId = box->pStart[cellId];
        int offset = box->cellOffset[cellId];

        box->pList[startId+offset] = i;  
        box->cellOffset[cellId] +=1;
    }
}

void update_pList_withSort(ParticleSystem *p, ParticleSystem *tmpPs,BoundingBox *box){
    /* initialize */
    for (int i=0; i<box->N; i++){
        box->pNum[i] = 0;
    }

    for (int i=0; i<box->N; i++){
        box->pStart[i] = 0;
    }

    for(int i=0;i<p->N;i++){
        p->mortonOrder[i] = i;
    }

    /* get cellId and get mortonkey*/
    for (int i=0; i<p->N; i++){
        int dx = floor((p->x[i*DIM+0]-box->minx)*box->invdx)+1; //+1 for ghost cell
        int dy = floor((p->x[i*DIM+1]-box->miny)*box->invdy)+1; //+1 for ghost cell
        int dz = floor((p->x[i*DIM+2]-box->minz)*box->invdz)+1; //+1 for ghost cell
        p->cellx[i*DIM+0] = dx;
        p->cellx[i*DIM+1] = dy;
        p->cellx[i*DIM+2] = dz;
        p->cellId[i] = (box->sizey*dz+dy)*box->sizex+dx;
        p->mortonKey[i]=morton3D((uint32_t)dx,(uint32_t)dy,(uint32_t)dz);
    }

    /* ========= sort and reorder ========= */
    radixSortUint32(&p->mortonKey,&p->mortonOrder,p->N,p->tmpMortonKey,p->tmpMortonOrder);

    swap_ps(p, tmpPs);



    /* count number of particle in cell */

    for (int i=0; i<p->N; i++){
        if (p->isActive[i]!=1){
            continue;
        }
        int cellId = p->cellId[i];
        box->pNum[cellId] +=1;
    }

    box->pStart[0]=0;
    for (int i=1; i<box->N; i++){
        box->pStart[i]=box->pNum[i-1]+box->pStart[i-1];
    }

    for (int i=0; i<box->N; i++){
        box->cellOffset[i] = 0;
    }

    for (int i=0; i<p->N; i++){
        if (p->isActive[i]!=1){
            continue;
        }
        int cellId = p->cellId[i];
        int startId = box->pStart[cellId];
        int offset = box->cellOffset[cellId];


        box->pList[startId+offset] = i;  
        box->cellOffset[cellId] +=1;
    }

}

void update_tList(BoundingBox *box, TriangleMesh *mesh){
    for (int i=0; i<mesh->nTri; i++){
        int sx = floor((mesh->minx[i]-box->minx)*box->invdx)+1;
        int sy = floor((mesh->miny[i]-box->miny)*box->invdy)+1;
        int sz = floor((mesh->minz[i]-box->minz)*box->invdz)+1;

        int ex = ceil((mesh->maxx[i]-box->minx)*box->invdx)+1;
        int ey = ceil((mesh->maxy[i]-box->miny)*box->invdy)+1;
        int ez = ceil((mesh->maxz[i]-box->minz)*box->invdz)+1;
        printf("%d\n",i);

        for (int dz=sz; dz<=ez; dz++){
            for (int dy=sy; dy<=ey; dy++){
                for (int dx=sx; dx<=ex; dx++){
                    int cellId= (box->sizey*dz+dy)*box->sizex+dx;
                    if(cellId>=box->N){
                        printf("Cell Id is FUCKED BRO!");
                    }

                    if (box->tNum[cellId]>= box->MAX_TRI){
                        printf("!!!!!!!Too many triangles in a cell!!!!\n");
                    }

                    int tNum = box->tNum[cellId];
                    box->tList[cellId*box->MAX_TRI+tNum]=i;
                    box->tNum[cellId]++;
                }
            }
        }
    }
}

void update_pList_fast(ParticleSystem *p, BoundingBox *box){

    /* initialize */
    for (int i=0; i<box->numUsedCells; i++){
        int cid = box->usedCells[i];
        box->pNum[cid] = 0;
        box->pStart[cid] = 0;
        box->cellOffset[cid] = 0;
    }


    /* get cellId*/
    for (int i=0; i<p->N; i++){
        if (p->isActive[i]!=1){
            continue;
        }
        int dx = floor((p->x[i*DIM+0]-box->minx)*box->invdx)+1; //+1 for ghost cell
        int dy = floor((p->x[i*DIM+1]-box->miny)*box->invdy)+1; //+1 for ghost cell
        int dz = floor((p->x[i*DIM+2]-box->minz)*box->invdz)+1; //+1 for ghost cell
        p->cellx[i*DIM+0] = dx;
        p->cellx[i*DIM+1] = dy;
        p->cellx[i*DIM+2] = dz;
        p->cellId[i] = (box->sizey*dz+dy)*box->sizex+dx;

    }

    /* count number of particle in cell */

    box->numUsedCells = 0;

    for (int i=0; i<p->N; i++){
        if (p->isActive[i]!=1){
            continue;
        }
        int cellId = p->cellId[i];

        if (box->pNum[cellId]==0){
            box->usedCells[box->numUsedCells] = cellId;
            box->numUsedCells += 1;
        }

        box->pNum[cellId] +=1;
    }

    //qsort(box->usedCells, box->numUsedCells,sizeof(int), compare_int );

    int runningSum = 0;
    for (int i=0; i<box->numUsedCells; i++){
        int cid = box->usedCells[i];
        box->pStart[cid]=runningSum;
        runningSum += box->pNum[cid];
    }


    for (int i=0; i<p->N; i++){
        if (p->isActive[i]!=1){
            continue;
        }
        int cellId = p->cellId[i];
        int startId = box->pStart[cellId];
        int offset = box->cellOffset[cellId];

        box->pList[startId+offset] = i;  
        box->cellOffset[cellId] +=1;
    }

}

void update_pList(ParticleSystem *p, BoundingBox *box){
    /* initialize */
    for (int i=0; i<box->N; i++){
        box->pNum[i] = 0;
    }

    for (int i=0; i<box->N; i++){
        box->pStart[i] = 0;
    }

    /* get cellId*/
    for (int i=0; i<p->N; i++){
        if (p->isActive[i]!=1){
            continue;
        }
        int dx = floor((p->x[i*DIM+0]-box->minx)*box->invdx)+1; //+1 for ghost cell
        int dy = floor((p->x[i*DIM+1]-box->miny)*box->invdy)+1; //+1 for ghost cell
        int dz = floor((p->x[i*DIM+2]-box->minz)*box->invdz)+1; //+1 for ghost cell
        p->cellx[i*DIM+0] = dx;
        p->cellx[i*DIM+1] = dy;
        p->cellx[i*DIM+2] = dz;
        p->cellId[i] = (box->sizey*dz+dy)*box->sizex+dx;
    }

    /* count number of particle in cell */

    for (int i=0; i<p->N; i++){
        if (p->isActive[i]!=1){
            continue;
        }
        int cellId = p->cellId[i];
        box->pNum[cellId] +=1;
    }

    box->pStart[0]=0;
    for (int i=1; i<box->N; i++){
        box->pStart[i]=box->pNum[i-1]+box->pStart[i-1];
    }

    for (int i=0; i<box->N; i++){
        box->cellOffset[i] = 0;
    }

    for (int i=0; i<p->N; i++){
        if (p->isActive[i]!=1){
            continue;
        }
        int cellId = p->cellId[i];
        int startId = box->pStart[cellId];
        int offset = box->cellOffset[cellId];


        box->pList[startId+offset] = i;  
        box->cellOffset[cellId] +=1;
    }
}

/*
   ===================================
   device  related functions
   ===================================
 */

__global__ void dk_build_cellCount(DeviceParticleGroup* p, DeviceBoundingBox* box){

    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if(i>=p->N || p->isActive[i]!=1){
        return;
    }

    int cellId = d_calcCellId(p,i,box);
    p->cellId[i] = cellId;

    atomicAdd(&box->pNum[cellId],1);

}

__global__ void dk_build_pList(DeviceParticleGroup* p, DeviceBoundingBox* box){

    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if(i>=p->N || p->isActive[i]!=1){
        return;
    }

    int cellId = p->cellId[i];
    int offset = atomicAdd(&box->cellOffset[cellId],1);
    int index = box->pStart[cellId] + offset;
    box->pList[index] = i;
}

__device__ int d_calcCellId(DeviceParticleGroup* p,int i, DeviceBoundingBox* box){
    int dx = floor((p->x[i*DIM+0]-box->minx)*box->invdx)+1; //+1 for ghost cell
    int dy = floor((p->x[i*DIM+1]-box->miny)*box->invdy)+1; //+1 for ghost cell
    int dz = floor((p->x[i*DIM+2]-box->minz)*box->invdz)+1; //+1 for ghost cell
    p->cellx[i*DIM+0] = dx;
    p->cellx[i*DIM+1] = dy;
    p->cellx[i*DIM+2] = dz;

    return (box->sizey*dz+dy)*box->sizex+dx;
}

void d_update_pList(ParticleSystem *p, BoundingBox *box,int gridSize, int blockSize){
    // printf("N=%d MAX_NEI=%d numCell=%d\n", p->d_group.N, p->d_group.MAX_NEI, box->N);

    /* initialize */
    cudaMemset(box->d_box.pNum, 0, sizeof(int)*box->N);
    cudaMemset(box->d_box.pStart, 0, sizeof(int)*box->N);
    cudaMemset(box->d_box.cellOffset, 0, sizeof(int)*box->N);

    dk_build_cellCount<<<gridSize,blockSize>>>(p->d_groupPtr,box->d_boxPtr);

    cub::DeviceScan::ExclusiveSum(box->d_box.tmpExSum,
            box->d_box.scanTmpBytes,
            box->d_box.pNum,
            box->d_box.pStart,
            box->d_box.N);

    dk_build_pList<<<gridSize,blockSize>>>(p->d_groupPtr,box->d_boxPtr);

}
