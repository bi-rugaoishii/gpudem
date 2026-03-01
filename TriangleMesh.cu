#include "TriangleMesh.h"

/* ========== hash key related functions used for merging ===========*/
static inline unsigned long long hash_func(long long x,long long y,long long z){
    unsigned long long h = 1469598103934665603ULL;
    h ^= (unsigned long long)x; h *= 1099511628211ULL;
    h ^= (unsigned long long)y; h *= 1099511628211ULL;
    h ^= (unsigned long long)z; h *= 1099511628211ULL;
    return h;
}

void vertex_hash_init(VertexHash* h, size_t size){
    h->size = size;
    h->table = (HashEntry*)calloc(size, sizeof(HashEntry));
}

int vertex_hash_get_or_insert(VertexHash* h,
                              double x, double y, double z,
                              double eps,
                              int vIndex){
    long long ix = llround(x / eps);
    long long iy = llround(y / eps);
    long long iz = llround(z / eps);

    unsigned long long hash = hash_func(ix,iy,iz);
    size_t pos = hash % h->size;

    while(1){
        if(!h->table[pos].used){
            /* 新規登録 */
            h->table[pos].used  = 1;
            h->table[pos].ix    = ix;
            h->table[pos].iy    = iy;
            h->table[pos].iz    = iz;
            h->table[pos].index = vIndex;

            return h->table[pos].index;
        }

        /* 一致チェック */
        if(h->table[pos].ix == ix &&
           h->table[pos].iy == iy &&
           h->table[pos].iz == iz)
        {
            return h->table[pos].index;
        }

        /* 線形探索 */
        pos = (pos + 1) % h->size;
    }
}

int get_edge_index(int v0, int v1, int* edge,int eIndex){

    int tv0 = (v0<v1) ? v0:v1;
    int tv1 = (v0>v1) ? v0:v1;
    for (int k=0; k<eIndex; k++){
        if (edge[k*2+0]==v0 && edge[k*2+1]==v1){
            return k;
            break;
        }
    }
    edge[eIndex*2+0]= tv0;
    edge[eIndex*2+1]= tv1;
    return eIndex;

}

void sort3(int *x, int *y, int *z){
    int tmp;

    if (*x > *y){
        tmp = *x; *x = *y; *y = tmp;
    }

    if (*y > *z){
        tmp = *y; *y = *z; *z = tmp;
    }

    if (*x > *y){
        tmp = *x; *x = *y; *y = tmp;
    }
}

/* ========== Device Related===============*/
#if USE_GPU
void deviceMallocCopyTriangleMesh(TriangleMesh *mesh){
    int Nv = mesh->nVert;
    int Nt = mesh->nTri;

    mesh->d_mesh.nVert = mesh->nVert;
    mesh->d_mesh.nTri = mesh->nTri;
    mesh->d_mesh.nShift = mesh->nShift;


    cudaMalloc(&mesh->d_mesh.mx, sizeof(double)*Nv);
    cudaMalloc(&mesh->d_mesh.my, sizeof(double)*Nv);
    cudaMalloc(&mesh->d_mesh.mz, sizeof(double)*Nv);

    cudaMalloc(&mesh->d_mesh.nx, sizeof(double)*Nt);
    cudaMalloc(&mesh->d_mesh.ny, sizeof(double)*Nt);
    cudaMalloc(&mesh->d_mesh.nz, sizeof(double)*Nt);
    cudaMalloc(&mesh->d_mesh.d,  sizeof(double)*Nt);

    cudaMalloc(&mesh->d_mesh.edge, sizeof(int)*Nv*2);

    cudaMalloc(&mesh->d_mesh.e01x, sizeof(double)*Nt);
    cudaMalloc(&mesh->d_mesh.e01y, sizeof(double)*Nt);
    cudaMalloc(&mesh->d_mesh.e01z, sizeof(double)*Nt);

    cudaMalloc(&mesh->d_mesh.e02x, sizeof(double)*Nt);
    cudaMalloc(&mesh->d_mesh.e02y, sizeof(double)*Nt);
    cudaMalloc(&mesh->d_mesh.e02z, sizeof(double)*Nt);

    cudaMalloc(&mesh->d_mesh.e12x, sizeof(double)*Nt);
    cudaMalloc(&mesh->d_mesh.e12y, sizeof(double)*Nt);
    cudaMalloc(&mesh->d_mesh.e12z, sizeof(double)*Nt);

    cudaMalloc(&mesh->d_mesh.d00,    sizeof(double)*Nt);
    cudaMalloc(&mesh->d_mesh.d00inv, sizeof(double)*Nt);
    cudaMalloc(&mesh->d_mesh.d01,    sizeof(double)*Nt);
    cudaMalloc(&mesh->d_mesh.d11,    sizeof(double)*Nt);
    cudaMalloc(&mesh->d_mesh.d11inv, sizeof(double)*Nt);
    cudaMalloc(&mesh->d_mesh.d22,    sizeof(double)*Nt);
    cudaMalloc(&mesh->d_mesh.d22inv, sizeof(double)*Nt);
    cudaMalloc(&mesh->d_mesh.denom,  sizeof(double)*Nt);

    cudaMalloc(&mesh->d_mesh.tri_i0, sizeof(int)*Nt);
    cudaMalloc(&mesh->d_mesh.tri_i1, sizeof(int)*Nt);
    cudaMalloc(&mesh->d_mesh.tri_i2, sizeof(int)*Nt);

    cudaMalloc(&mesh->d_mesh.tri_e0, sizeof(int)*Nt);
    cudaMalloc(&mesh->d_mesh.tri_e1, sizeof(int)*Nt);
    cudaMalloc(&mesh->d_mesh.tri_e2, sizeof(int)*Nt);

    cudaMalloc(&mesh->d_mesh.minx, sizeof(double)*Nt);
    cudaMalloc(&mesh->d_mesh.maxx, sizeof(double)*Nt);

    cudaMalloc(&mesh->d_mesh.miny, sizeof(double)*Nt);
    cudaMalloc(&mesh->d_mesh.maxy, sizeof(double)*Nt);

    cudaMalloc(&mesh->d_mesh.minz, sizeof(double)*Nt);
    cudaMalloc(&mesh->d_mesh.maxz, sizeof(double)*Nt);

    /* ================= malloc struct =============*/
    cudaMalloc(&mesh->d_meshPtr, sizeof(DeviceTriangleMesh));


    /* ================ memcpy ================= */

    /* ---- vertex ---- */
    cudaMemcpy(mesh->d_mesh.mx, mesh->mx, sizeof(double)*Nv, cudaMemcpyHostToDevice);
    cudaMemcpy(mesh->d_mesh.my, mesh->my, sizeof(double)*Nv, cudaMemcpyHostToDevice);
    cudaMemcpy(mesh->d_mesh.mz, mesh->mz, sizeof(double)*Nv, cudaMemcpyHostToDevice);

    /* ---- normal & plane ---- */
    cudaMemcpy(mesh->d_mesh.nx, mesh->nx, sizeof(double)*Nt, cudaMemcpyHostToDevice);
    cudaMemcpy(mesh->d_mesh.ny, mesh->ny, sizeof(double)*Nt, cudaMemcpyHostToDevice);
    cudaMemcpy(mesh->d_mesh.nz, mesh->nz, sizeof(double)*Nt, cudaMemcpyHostToDevice);
    cudaMemcpy(mesh->d_mesh.d,  mesh->d,  sizeof(double)*Nt, cudaMemcpyHostToDevice);

    /* ---- edges (index) ---- */
    cudaMemcpy(mesh->d_mesh.edge, mesh->edge, sizeof(int)*Nv*2, cudaMemcpyHostToDevice);

    /* ---- triangle edge vectors ---- */
    cudaMemcpy(mesh->d_mesh.e01x, mesh->e01x, sizeof(double)*Nt, cudaMemcpyHostToDevice);
    cudaMemcpy(mesh->d_mesh.e01y, mesh->e01y, sizeof(double)*Nt, cudaMemcpyHostToDevice);
    cudaMemcpy(mesh->d_mesh.e01z, mesh->e01z, sizeof(double)*Nt, cudaMemcpyHostToDevice);

    cudaMemcpy(mesh->d_mesh.e02x, mesh->e02x, sizeof(double)*Nt, cudaMemcpyHostToDevice);
    cudaMemcpy(mesh->d_mesh.e02y, mesh->e02y, sizeof(double)*Nt, cudaMemcpyHostToDevice);
    cudaMemcpy(mesh->d_mesh.e02z, mesh->e02z, sizeof(double)*Nt, cudaMemcpyHostToDevice);

    cudaMemcpy(mesh->d_mesh.e12x, mesh->e12x, sizeof(double)*Nt, cudaMemcpyHostToDevice);
    cudaMemcpy(mesh->d_mesh.e12y, mesh->e12y, sizeof(double)*Nt, cudaMemcpyHostToDevice);
    cudaMemcpy(mesh->d_mesh.e12z, mesh->e12z, sizeof(double)*Nt, cudaMemcpyHostToDevice);

    /* ---- barycentric precompute ---- */
    cudaMemcpy(mesh->d_mesh.d00,    mesh->d00,    sizeof(double)*Nt, cudaMemcpyHostToDevice);
    cudaMemcpy(mesh->d_mesh.d00inv, mesh->d00inv, sizeof(double)*Nt, cudaMemcpyHostToDevice);
    cudaMemcpy(mesh->d_mesh.d01,    mesh->d01,    sizeof(double)*Nt, cudaMemcpyHostToDevice);
    cudaMemcpy(mesh->d_mesh.d11,    mesh->d11,    sizeof(double)*Nt, cudaMemcpyHostToDevice);
    cudaMemcpy(mesh->d_mesh.d11inv, mesh->d11inv, sizeof(double)*Nt, cudaMemcpyHostToDevice);
    cudaMemcpy(mesh->d_mesh.d22,    mesh->d22,    sizeof(double)*Nt, cudaMemcpyHostToDevice);
    cudaMemcpy(mesh->d_mesh.d22inv, mesh->d22inv, sizeof(double)*Nt, cudaMemcpyHostToDevice);
    cudaMemcpy(mesh->d_mesh.denom,  mesh->denom,  sizeof(double)*Nt, cudaMemcpyHostToDevice);

    /* ---- triangle vertex index ---- */
    cudaMemcpy(mesh->d_mesh.tri_i0, mesh->tri_i0, sizeof(int)*Nt, cudaMemcpyHostToDevice);
    cudaMemcpy(mesh->d_mesh.tri_i1, mesh->tri_i1, sizeof(int)*Nt, cudaMemcpyHostToDevice);
    cudaMemcpy(mesh->d_mesh.tri_i2, mesh->tri_i2, sizeof(int)*Nt, cudaMemcpyHostToDevice);

    /* ---- triangle edge index ---- */
    cudaMemcpy(mesh->d_mesh.tri_e0, mesh->tri_e0, sizeof(int)*Nt, cudaMemcpyHostToDevice);
    cudaMemcpy(mesh->d_mesh.tri_e1, mesh->tri_e1, sizeof(int)*Nt, cudaMemcpyHostToDevice);
    cudaMemcpy(mesh->d_mesh.tri_e2, mesh->tri_e2, sizeof(int)*Nt, cudaMemcpyHostToDevice);

    /* ---- triangle AABB ---- */
    cudaMemcpy(mesh->d_mesh.minx, mesh->minx, sizeof(double)*Nt, cudaMemcpyHostToDevice);
    cudaMemcpy(mesh->d_mesh.maxx, mesh->maxx, sizeof(double)*Nt, cudaMemcpyHostToDevice);
    cudaMemcpy(mesh->d_mesh.miny, mesh->miny, sizeof(double)*Nt, cudaMemcpyHostToDevice);
    cudaMemcpy(mesh->d_mesh.maxy, mesh->maxy, sizeof(double)*Nt, cudaMemcpyHostToDevice);
    cudaMemcpy(mesh->d_mesh.minz, mesh->minz, sizeof(double)*Nt, cudaMemcpyHostToDevice);
    cudaMemcpy(mesh->d_mesh.maxz, mesh->maxz, sizeof(double)*Nt, cudaMemcpyHostToDevice);

    /* --- copy struct ----*/
    cudaMemcpy(mesh->d_meshPtr,  &mesh->d_mesh,  sizeof(DeviceTriangleMesh), cudaMemcpyHostToDevice);

}


#endif

/* ========== triangle mesh ===============*/
void free_TriangleMesh(TriangleMesh* mesh){
    free(mesh->mx);
    free(mesh->my);
    free(mesh->mz);

    free(mesh->nx);
    free(mesh->ny);
    free(mesh->nz);
    free(mesh->d);

    free(mesh->edge);


    free(mesh->e01x);
    free(mesh->e01y);
    free(mesh->e01z);

    free(mesh->e02x);
    free(mesh->e02y);
    free(mesh->e02z);

    free(mesh->e12x);
    free(mesh->e12y);
    free(mesh->e12z);

    free(mesh->d00);
    free(mesh->d00inv);
    free(mesh->d01);
    free(mesh->d11);
    free(mesh->d11inv);
    free(mesh->d22);
    free(mesh->d22inv);
    free(mesh->denom);

    free(mesh->minx);
    free(mesh->miny);
    free(mesh->minz);

    free(mesh->maxx);
    free(mesh->maxy);
    free(mesh->maxz);


    free(mesh->tri_i0);
    free(mesh->tri_i1);
    free(mesh->tri_i2);

    free(mesh->tri_e0);
    free(mesh->tri_e1);
    free(mesh->tri_e2);

    #if USE_GPU
    cudaFree(mesh->d_mesh.mx);
    cudaFree(mesh->d_mesh.my);
    cudaFree(mesh->d_mesh.mz);

    cudaFree(mesh->d_mesh.nx);
    cudaFree(mesh->d_mesh.ny);
    cudaFree(mesh->d_mesh.nz);
    cudaFree(mesh->d_mesh.d);

    cudaFree(mesh->d_mesh.edge);

    cudaFree(mesh->d_mesh.e01x);
    cudaFree(mesh->d_mesh.e01y);
    cudaFree(mesh->d_mesh.e01z);

    cudaFree(mesh->d_mesh.e02x);
    cudaFree(mesh->d_mesh.e02y);
    cudaFree(mesh->d_mesh.e02z);

    cudaFree(mesh->d_mesh.e12x);
    cudaFree(mesh->d_mesh.e12y);
    cudaFree(mesh->d_mesh.e12z);

    cudaFree(mesh->d_mesh.d00);
    cudaFree(mesh->d_mesh.d00inv);
    cudaFree(mesh->d_mesh.d01);
    cudaFree(mesh->d_mesh.d11);
    cudaFree(mesh->d_mesh.d11inv);
    cudaFree(mesh->d_mesh.d22);
    cudaFree(mesh->d_mesh.d22inv);
    cudaFree(mesh->d_mesh.denom);

    cudaFree(mesh->d_mesh.tri_i0);
    cudaFree(mesh->d_mesh.tri_i1);
    cudaFree(mesh->d_mesh.tri_i2);

    cudaFree(mesh->d_mesh.tri_e0);
    cudaFree(mesh->d_mesh.tri_e1);
    cudaFree(mesh->d_mesh.tri_e2);

    cudaFree(mesh->d_mesh.minx);
    cudaFree(mesh->d_mesh.maxx);
    cudaFree(mesh->d_mesh.miny);
    cudaFree(mesh->d_mesh.maxy);
    cudaFree(mesh->d_mesh.minz);
    cudaFree(mesh->d_mesh.maxz);

    cudaFree(mesh->d_meshPtr);
    #endif
}

static int count_ascii_stl_triangles(FILE* fp){
    char line[256];
    int count = 0;

    while(fgets(line,256,fp))
    {
        if(strncmp(line,"facet normal",12)==0 ||
                strncmp(line," facet normal",13)==0)
        {
            count++;
        }
    }
    return count;
}


int load_ascii_stl_double(const char* filename, TriangleMesh* mesh){
    FILE* fp = fopen(filename,"r");
    if(!fp){
        printf("stl open error\n");
        return 0;
    }

    int triCount = count_ascii_stl_triangles(fp);
    rewind(fp);

    mesh->nTri = triCount;
    mesh->nVert  = triCount * 3;

    int Nv = mesh->nVert;
    int Nt = mesh->nTri;
    /* ===== global mesh AABB init ===== */

    mesh->gminx =  1e300;
    mesh->gminy =  1e300;
    mesh->gminz =  1e300;

    mesh->gmaxx = -1e300;
    mesh->gmaxy = -1e300;
    mesh->gmaxz = -1e300;

    /* allocate SoA */
    mesh->mx = (double*)malloc(sizeof(double)*Nv);
    mesh->my = (double*)malloc(sizeof(double)*Nv);
    mesh->mz = (double*)malloc(sizeof(double)*Nv);

    mesh->nx = (double*)malloc(sizeof(double)*Nt);
    mesh->ny = (double*)malloc(sizeof(double)*Nt);
    mesh->nz = (double*)malloc(sizeof(double)*Nt);
    mesh->d = (double*)malloc(sizeof(double)*Nt);

    mesh->edge = (int*)malloc(sizeof(int)*Nv*2);

    mesh->e01x = (double*)malloc(sizeof(double)*Nt);
    mesh->e01y = (double*)malloc(sizeof(double)*Nt);
    mesh->e01z = (double*)malloc(sizeof(double)*Nt);

    mesh->e02x = (double*)malloc(sizeof(double)*Nt);
    mesh->e02y = (double*)malloc(sizeof(double)*Nt);
    mesh->e02z = (double*)malloc(sizeof(double)*Nt);

    mesh->e12x = (double*)malloc(sizeof(double)*Nt);
    mesh->e12y = (double*)malloc(sizeof(double)*Nt);
    mesh->e12z = (double*)malloc(sizeof(double)*Nt);

    mesh->d00 = (double*)malloc(sizeof(double)*Nt);
    mesh->d00inv = (double*)malloc(sizeof(double)*Nt);
    mesh->d01 = (double*)malloc(sizeof(double)*Nt);
    mesh->d11 = (double*)malloc(sizeof(double)*Nt);
    mesh->d11inv = (double*)malloc(sizeof(double)*Nt);
    mesh->d22 = (double*)malloc(sizeof(double)*Nt);
    mesh->d22inv = (double*)malloc(sizeof(double)*Nt);
    mesh->denom = (double*)malloc(sizeof(double)*Nt);


    mesh->tri_i0 = (int*)malloc(sizeof(int)*Nt);
    mesh->tri_i1 = (int*)malloc(sizeof(int)*Nt);
    mesh->tri_i2 = (int*)malloc(sizeof(int)*Nt);

    mesh->tri_e0 = (int*)malloc(sizeof(int)*Nt);
    mesh->tri_e1 = (int*)malloc(sizeof(int)*Nt);
    mesh->tri_e2 = (int*)malloc(sizeof(int)*Nt);

    mesh->minx = (double*)malloc(sizeof(double)*Nt);
    mesh->maxx = (double*)malloc(sizeof(double)*Nt);

    mesh->miny = (double*)malloc(sizeof(double)*Nt);
    mesh->maxy = (double*)malloc(sizeof(double)*Nt);

    mesh->minz = (double*)malloc(sizeof(double)*Nt);
    mesh->maxz = (double*)malloc(sizeof(double)*Nt);

    char line[256];

    int triIndex = 0;
    int vIndex   = 0;
    int eIndex = 0;
    mesh->nShift   = Nv;

    double nx,ny,nz;
    double x,y,z;

    /* ========= for merging ============= */
    VertexHash vh;
    vertex_hash_init(&vh,mesh->nVert*3);

    while(fgets(line,256,fp))
    {
        //if(sscanf(line," facet normal %lf %lf %lf",&nx,&ny,&nz)==3 ||       sscanf(line,"facet normal %lf %lf %lf",&nx,&ny,&nz)==3)
        /* need better scan */
        if(sscanf(line," facet normal %lf %lf %lf",&nx,&ny,&nz)==3){


            /* outer loop */
            if(!fgets(line,256,fp)) break;
            int v0,v1,v2;

            /* v0 */
            if(!fgets(line,256,fp)) break;
            sscanf(line," vertex %lf %lf %lf",&x,&y,&z);
            mesh->mx[vIndex]=x;
            mesh->my[vIndex]=y;
            mesh->mz[vIndex]=z;
            v0=vertex_hash_get_or_insert(&vh,x,y,z,1e-9,vIndex);
            if (v0==vIndex){
                vIndex   += 1; //was new  
            }



            /* v1 */
            if(!fgets(line,256,fp)) break;
            sscanf(line," vertex %lf %lf %lf",&x,&y,&z);
            mesh->mx[vIndex]=x;
            mesh->my[vIndex]=y;
            mesh->mz[vIndex]=z;
            v1=vertex_hash_get_or_insert(&vh,x,y,z,1e-9,vIndex);
            if (v1==vIndex){
                vIndex   += 1; //was new  
            }

            /* v2 */
            if(!fgets(line,256,fp)) break;
            sscanf(line," vertex %lf %lf %lf",&x,&y,&z);
            mesh->mx[vIndex]=x;
            mesh->my[vIndex]=y;
            mesh->mz[vIndex]=z;
            v2=vertex_hash_get_or_insert(&vh,x,y,z,1e-9,vIndex);
            if (v2==vIndex){
                vIndex   += 1; //was new  
            }

            /* reorder vertices in ascending order */
            sort3(&v0,&v1,&v2);

            /* make edge */
            int tmpEdge;
            tmpEdge = get_edge_index(v0,v1,mesh->edge,eIndex);
            if(tmpEdge==eIndex){ //new index
                eIndex += 1;
            }
            mesh->tri_e0[triIndex] = tmpEdge+mesh->nShift;

            tmpEdge = get_edge_index(v1,v2,mesh->edge,eIndex);
            if(tmpEdge==eIndex){ //new index
                eIndex += 1;
            }
            mesh->tri_e1[triIndex] = tmpEdge+mesh->nShift;

            tmpEdge = get_edge_index(v0,v2,mesh->edge,eIndex);
            if(tmpEdge==eIndex){ //new index
                eIndex += 1;
            }
            mesh->tri_e2[triIndex] = tmpEdge+mesh->nShift;



            /* calculate edge vectors */
            Vec3 e01;

            e01.x = mesh->mx[v1]-mesh->mx[v0];
            e01.y = mesh->my[v1]-mesh->my[v0];
            e01.z = mesh->mz[v1]-mesh->mz[v0];

            Vec3 e02;
            e02.x = mesh->mx[v2]-mesh->mx[v0];
            e02.y = mesh->my[v2]-mesh->my[v0];
            e02.z = mesh->mz[v2]-mesh->mz[v0];

            Vec3 e12;
            e12.x = mesh->mx[v2]-mesh->mx[v1];
            e12.y = mesh->my[v2]-mesh->my[v1];
            e12.z = mesh->mz[v2]-mesh->mz[v1];

            mesh->e01x[triIndex] = e01.x;
            mesh->e01y[triIndex] = e01.y;
            mesh->e01z[triIndex] = e01.z;

            mesh->e02x[triIndex] = e02.x;
            mesh->e02y[triIndex] = e02.y;
            mesh->e02z[triIndex] = e02.z;

            mesh->e12x[triIndex] = e12.x;
            mesh->e12y[triIndex] = e12.y;
            mesh->e12z[triIndex] = e12.z;

            mesh->d00[triIndex] = vdot(e01,e01);
            mesh->d00inv[triIndex] = 1./mesh->d00[triIndex];
            mesh->d01[triIndex] = vdot(e01,e02);
            mesh->d11[triIndex] = vdot(e02,e02);
            mesh->d11inv[triIndex] = 1./mesh->d11[triIndex];
            mesh->d22[triIndex] = vdot(e12,e12);
            mesh->d22inv[triIndex] = 1./mesh->d22[triIndex];

            mesh->denom[triIndex] = mesh->d00[triIndex]*mesh->d11[triIndex]-mesh->d01[triIndex]*mesh->d01[triIndex];
            mesh->denom[triIndex] = 1./mesh->denom[triIndex];

            /* calculate normals */
            Vec3 n;
            n = vcross(e01,e02);
            n = vscalar(1./sqrt(vdot(n,n)),n);
            mesh->nx[triIndex] = n.x;
            mesh->ny[triIndex] = n.y;
            mesh->nz[triIndex] = n.z;
            mesh->d[triIndex]=-n.x*x-n.y*y-n.z*z;

            mesh->tri_i0[triIndex]=v0;
            mesh->tri_i1[triIndex]=v1;
            mesh->tri_i2[triIndex]=v2;

            /* ===== compute triangle AABB ===== */

            double x0 = mesh->mx[v0];
            double y0 = mesh->my[v0];
            double z0 = mesh->mz[v0];

            double x1 = mesh->mx[v1];
            double y1 = mesh->my[v1];
            double z1 = mesh->mz[v1];

            double x2 = mesh->mx[v2];
            double y2 = mesh->my[v2];
            double z2 = mesh->mz[v2];

            /* min */
            double minx = x0;
            if(x1 < minx) minx = x1;
            if(x2 < minx) minx = x2;

            double miny = y0;
            if(y1 < miny) miny = y1;
            if(y2 < miny) miny = y2;

            double minz = z0;
            if(z1 < minz) minz = z1;
            if(z2 < minz) minz = z2;

            /* max */
            double maxx = x0;
            if(x1 > maxx) maxx = x1;
            if(x2 > maxx) maxx = x2;

            double maxy = y0;
            if(y1 > maxy) maxy = y1;
            if(y2 > maxy) maxy = y2;

            double maxz = z0;
            if(z1 > maxz) maxz = z1;
            if(z2 > maxz) maxz = z2;

            mesh->minx[triIndex] = minx;
            mesh->maxx[triIndex] = maxx;

            mesh->miny[triIndex] = miny;
            mesh->maxy[triIndex] = maxy;

            mesh->minz[triIndex] = minz;
            mesh->maxz[triIndex] = maxz;

            /* ===== update global mesh AABB ===== */

            if(minx < mesh->gminx) mesh->gminx = minx;
            if(miny < mesh->gminy) mesh->gminy = miny;
            if(minz < mesh->gminz) mesh->gminz = minz;

            if(maxx > mesh->gmaxx) mesh->gmaxx = maxx;
            if(maxy > mesh->gmaxy) mesh->gmaxy = maxy;
            if(maxz > mesh->gmaxz) mesh->gmaxz = maxz;

            triIndex += 1;

            if(!fgets(line,256,fp)) break;
            if(!fgets(line,256,fp)) break;
        }
    }

    printf(" final number of edges are %d \n", eIndex);

    fclose(fp);

    printf("ASCII STL (double) loaded : %d triangles\n",triCount);
    return 1;
}

