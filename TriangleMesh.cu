#include "TriangleMesh.h"

void free_TriangleMesh(TriangleMesh* mesh){
    free(mesh->mx);
    free(mesh->my);
    free(mesh->mz);

    free(mesh->nx);
    free(mesh->ny);
    free(mesh->nz);

    free(mesh->e01x);
    free(mesh->e01y);
    free(mesh->e01z);

    free(mesh->e02x);
    free(mesh->e02y);
    free(mesh->e02z);

    free(mesh->d00);
    free(mesh->d01);
    free(mesh->d11);
    free(mesh->denom);

    free(mesh->minx);
    free(mesh->miny);
    free(mesh->minz);

    free(mesh->maxx);
    free(mesh->maxy);
    free(mesh->maxz);
    free(mesh->d);

    
    free(mesh->tri_i0);
    free(mesh->tri_i1);
    free(mesh->tri_i2);
}

static int count_ascii_stl_triangles(FILE* fp)
{
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

    
    mesh->e01x = (double*)malloc(sizeof(double)*Nt);
    mesh->e01y = (double*)malloc(sizeof(double)*Nt);
    mesh->e01z = (double*)malloc(sizeof(double)*Nt);

    mesh->e02x = (double*)malloc(sizeof(double)*Nt);
    mesh->e02y = (double*)malloc(sizeof(double)*Nt);
    mesh->e02z = (double*)malloc(sizeof(double)*Nt);

    mesh->d00 = (double*)malloc(sizeof(double)*Nt);
    mesh->d01 = (double*)malloc(sizeof(double)*Nt);
    mesh->d11 = (double*)malloc(sizeof(double)*Nt);
    mesh->denom = (double*)malloc(sizeof(double)*Nt);

    mesh->tri_i0 = (int*)malloc(sizeof(int)*Nt);
    mesh->tri_i1 = (int*)malloc(sizeof(int)*Nt);
    mesh->tri_i2 = (int*)malloc(sizeof(int)*Nt);

    mesh->minx = (double*)malloc(sizeof(double)*Nt);
    mesh->maxx = (double*)malloc(sizeof(double)*Nt);

    mesh->miny = (double*)malloc(sizeof(double)*Nt);
    mesh->maxy = (double*)malloc(sizeof(double)*Nt);

    mesh->minz = (double*)malloc(sizeof(double)*Nt);
    mesh->maxz = (double*)malloc(sizeof(double)*Nt);

    char line[256];

    int triIndex = 0;
    int vIndex   = 0;

    double nx,ny,nz;
    double x,y,z;

    while(fgets(line,256,fp))
    {
        if(sscanf(line," facet normal %lf %lf %lf",&nx,&ny,&nz)==3 ||
                sscanf(line,"facet normal %lf %lf %lf",&nx,&ny,&nz)==3)
        {
            /* normalize */
            double norm = nx*nx+ny*ny+nz*nz;
            norm=sqrt(norm);

            mesh->nx[triIndex] = nx/norm;
            mesh->ny[triIndex] = ny/norm;
            mesh->nz[triIndex] = nz/norm;


            /* outer loop */
            fgets(line,256,fp);

            /* v0 */
            fgets(line,256,fp);
            sscanf(line," vertex %lf %lf %lf",&x,&y,&z);
            mesh->mx[vIndex+0]=x;
            mesh->my[vIndex+0]=y;
            mesh->mz[vIndex+0]=z;

            /* v1 */
            fgets(line,256,fp);
            sscanf(line," vertex %lf %lf %lf",&x,&y,&z);
            mesh->mx[vIndex+1]=x;
            mesh->my[vIndex+1]=y;
            mesh->mz[vIndex+1]=z;

            /* v2 */
            fgets(line,256,fp);
            sscanf(line," vertex %lf %lf %lf",&x,&y,&z);
            mesh->mx[vIndex+2]=x;
            mesh->my[vIndex+2]=y;
            mesh->mz[vIndex+2]=z;

            mesh->d[triIndex]=-nx*x-ny*y-nz*z;

            /* calculate edges */
            mesh->e01x[triIndex] = mesh->mx[vIndex+1]-mesh->mx[vIndex+0];
            mesh->e01y[triIndex] = mesh->my[vIndex+1]-mesh->my[vIndex+0];
            mesh->e01z[triIndex] = mesh->mz[vIndex+1]-mesh->mz[vIndex+0];

            mesh->e02x[triIndex] = mesh->mx[vIndex+2]-mesh->mx[vIndex+0];
            mesh->e02y[triIndex] = mesh->my[vIndex+2]-mesh->my[vIndex+0];
            mesh->e02z[triIndex] = mesh->mz[vIndex+2]-mesh->mz[vIndex+0];

            mesh->d00[triIndex] = mesh->e01x[triIndex]*mesh->e01x[triIndex]+mesh->e01y[triIndex]*mesh->e01y[triIndex]+mesh->e01z[triIndex]*mesh->e01z[triIndex];
            mesh->d01[triIndex] = mesh->e01x[triIndex]*mesh->e02x[triIndex]+mesh->e01y[triIndex]*mesh->e02y[triIndex]+mesh->e01z[triIndex]*mesh->e02z[triIndex];
            mesh->d11[triIndex] = mesh->e02x[triIndex]*mesh->e02x[triIndex]+mesh->e02y[triIndex]*mesh->e02y[triIndex]+mesh->e02z[triIndex]*mesh->e02z[triIndex];

            mesh->denom[triIndex] = mesh->d00[triIndex]*mesh->d11[triIndex]-mesh->d01[triIndex]*mesh->d01[triIndex];
            mesh->denom[triIndex] = 1./mesh->denom[triIndex];

            mesh->tri_i0[triIndex]=vIndex+0;
            mesh->tri_i1[triIndex]=vIndex+1;
            mesh->tri_i2[triIndex]=vIndex+2;

            /* ===== compute triangle AABB ===== */

            double x0 = mesh->mx[vIndex+0];
            double y0 = mesh->my[vIndex+0];
            double z0 = mesh->mz[vIndex+0];

            double x1 = mesh->mx[vIndex+1];
            double y1 = mesh->my[vIndex+1];
            double z1 = mesh->mz[vIndex+1];

            double x2 = mesh->mx[vIndex+2];
            double y2 = mesh->my[vIndex+2];
            double z2 = mesh->mz[vIndex+2];

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

            vIndex   += 3;
            triIndex += 1;

            fgets(line,256,fp);
            fgets(line,256,fp);
        }
    }

    fclose(fp);

    printf("ASCII STL (double) loaded : %d triangles\n",triCount);
    return 1;
}
