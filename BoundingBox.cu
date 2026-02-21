#include "ParticleSystem.h"
#include "BoundingBox.h"



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
}

void free_BoundingBox(BoundingBox *box){
    free(box->pList);
    free(box->pNum);
    free(box->pStart);
    free(box->cellOffset);

    #if USE_GPU
        cudaFree(box->d_box.pList);
        cudaFree(box->d_box.tmpExSum);
        cudaFree(box->d_box.pNum);
        cudaFree(box->d_box.pStart);
        cudaFree(box->d_box.cellOffset);
    #endif
}

void initialize_BoundingBox(ParticleSystem *p, BoundingBox *box, double minx, double maxx, double miny, double maxy, double minz, double maxz){
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
    box->pNum = (int*)malloc(sizeof(int)*sizeBox);
    box->pStart = (int*)malloc(sizeof(int)*sizeBox);
    box->cellOffset = (int*)malloc(sizeof(int)*sizeBox);

    #if USE_GPU
        box->d_box.N= sizeBox;
        cudaMalloc(&box->d_box.pList, sizeof(int)*p->N);
        cudaMalloc(&box->d_box.pNum, sizeof(int)*sizeBox);
        cudaMalloc(&box->d_box.pStart, sizeof(int)*sizeBox);
        cudaMalloc(&box->d_box.cellOffset, sizeof(int)*sizeBox);
        cudaMalloc(&box->d_box.tmpExSum, sizeof(int)*sizeBox);

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

    #endif


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
        int cellId = p->cellId[i];
        box->pNum[cellId] +=1;
    }

    for (int i=1; i<box->N; i++){
        box->pStart[i]=box->pNum[i-1]+box->pStart[i-1];
    }

    for (int i=0; i<box->N; i++){
        box->cellOffset[i] = 0;
    }

    for (int i=0; i<p->N; i++){
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

__global__ void dk_build_cellCount(DeviceParticleGroup p, DeviceBoundingBox box){

    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if(i>=p.N){
        return;
    }

    int cellId = d_calcCellId(p,i,box);
    p.cellId[i] = cellId;

    atomicAdd(&box.pNum[cellId],1);

}

__global__ void dk_build_pList(DeviceParticleGroup p, DeviceBoundingBox box){

    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if(i>=p.N){
        return;
    }

    int cellId = p.cellId[i];
    int offset = atomicAdd(&box.cellOffset[cellId],1);
    int index = box.pStart[cellId] + offset;
    box.pList[index] = i;
}

__device__ int d_calcCellId(DeviceParticleGroup p,int i, DeviceBoundingBox box){
        int dx = floor((p.x[i*DIM+0]-box.minx)*box.invdx)+1; //+1 for ghost cell
        int dy = floor((p.x[i*DIM+1]-box.miny)*box.invdy)+1; //+1 for ghost cell
        int dz = floor((p.x[i*DIM+2]-box.minz)*box.invdz)+1; //+1 for ghost cell
        p.cellx[i*DIM+0] = dx;
        p.cellx[i*DIM+1] = dy;
        p.cellx[i*DIM+2] = dz;

        return (box.sizey*dz+dy)*box.sizex+dx;
}

void d_update_pList(ParticleSystem *p, BoundingBox *box,int gridSize, int blockSize){
    /* initialize */
    cudaMemset(box->d_box.pNum, 0, sizeof(int)*box->N);
    cudaMemset(box->d_box.pStart, 0, sizeof(int)*box->N);
    cudaMemset(box->d_box.cellOffset, 0, sizeof(int)*box->N);

    dk_build_cellCount<<<gridSize,blockSize>>>(p->d_group,box->d_box);

    cub::DeviceScan::ExclusiveSum(box->d_box.tmpExSum,
            box->d_box.scanTmpBytes,
            box->d_box.pNum,
            box->d_box.pStart,
            box->d_box.N);

    dk_build_pList<<<gridSize,blockSize>>>(p->d_group,box->d_box);

}
