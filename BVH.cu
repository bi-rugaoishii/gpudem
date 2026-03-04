#include "BVH.h"

typedef struct {
    int start;
    int end;   // [start,end)
    int node;
} BuildTask;

void initializeBVH(BVH *bvh, int numTriangles){
    int maxNodes = 2*numTriangles-1;
    bvh->minx = (double*)malloc(sizeof(double)*maxNodes);
    bvh->miny = (double*)malloc(sizeof(double)*maxNodes);
    bvh->minz = (double*)malloc(sizeof(double)*maxNodes);
    bvh->maxx = (double*)malloc(sizeof(double)*maxNodes);
    bvh->maxy = (double*)malloc(sizeof(double)*maxNodes);
    bvh->maxz = (double*)malloc(sizeof(double)*maxNodes);

    bvh->left = (int*)malloc(sizeof(int)*maxNodes);
    bvh->right = (int*)malloc(sizeof(int)*maxNodes);
    bvh->tri = (int*)malloc(sizeof(int)*maxNodes);
    bvh->nodeCount = 0;
}

void free_BVH(BVH *bvh){
    free(bvh->minx) ;
    free(bvh->miny) ;
    free(bvh->minz) ;
    free(bvh->maxx) ;
    free(bvh->maxy) ;
    free(bvh->maxz) ;

    free(bvh->left) ;
    free(bvh->right);
    free(bvh->tri) ;
}

void computeNodeAABB(int start, int end,TriangleMesh *mesh, BVH *bvh,int k){
    int startTri = mesh->sortedIndex[start];
    double minx = mesh->minx[startTri];
    double miny = mesh->miny[startTri];
    double minz = mesh->minz[startTri];
    double maxx = mesh->maxx[startTri];
    double maxy = mesh->maxy[startTri];
    double maxz = mesh->maxz[startTri];

    for(int i = start+1; i < end; i++){
        int sortedi = mesh->sortedIndex[i];
        if(mesh->minx[sortedi] < minx) minx = mesh->minx[sortedi];
        if(mesh->miny[sortedi] < miny) miny = mesh->miny[sortedi];
        if(mesh->minz[sortedi] < minz) minz = mesh->minz[sortedi];

        if(mesh->maxx[sortedi] > maxx) maxx = mesh->maxx[sortedi];
        if(mesh->maxy[sortedi] > maxy) maxy = mesh->maxy[sortedi];
        if(mesh->maxz[sortedi] > maxz) maxz = mesh->maxz[sortedi];
    }

    bvh->minx[k] = minx;
    bvh->miny[k] = miny;
    bvh->minz[k] = minz;

    bvh->maxx[k] = maxx;
    bvh->maxy[k] = maxy;
    bvh->maxz[k] = maxz;
}

int buildBVH(BVH* bvh, TriangleMesh *mesh){

    BuildTask* stack = (BuildTask*)malloc(sizeof(BuildTask) * (2*mesh->nTri));
    int sp = 0;

    int nodeCount = 0;

    int root = nodeCount++;
    stack[sp++] = (BuildTask){0, mesh->nTri, root};

    while(sp > 0){
        BuildTask t = stack[--sp];

        int start = t.start;
        int end   = t.end;
        int node  = t.node;

        computeNodeAABB(start, end,mesh, bvh, node);

        int count = end - start;

        if(count == 1)
        {
            bvh->left[node]  = -1;
            bvh->right[node] = -1;
            bvh->tri[node]   = mesh->sortedIndex[start];   // Morton順
            continue;
        }

        int mid = (start + end) / 2;

        int leftChild  = nodeCount++;
        int rightChild = nodeCount++;

        bvh->left[node]  = leftChild;
        bvh->right[node] = rightChild;
        bvh->tri[node]   = -1;

        stack[sp++] = (BuildTask){mid, end, rightChild};
        stack[sp++] = (BuildTask){start, mid, leftChild};
    }

    bvh->nodeCount = nodeCount;

    free(stack);
    return root;
}
