#include "TriangleMesh.h"

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

    /* allocate SoA */
    mesh->mx = (double*)malloc(sizeof(double)*Nv);
    mesh->my = (double*)malloc(sizeof(double)*Nv);
    mesh->mz = (double*)malloc(sizeof(double)*Nv);

    mesh->nx = (double*)malloc(sizeof(double)*Nt);
    mesh->ny = (double*)malloc(sizeof(double)*Nt);
    mesh->nz = (double*)malloc(sizeof(double)*Nt);

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
            mesh->nx[triIndex] = nx;
            mesh->ny[triIndex] = ny;
            mesh->nz[triIndex] = nz;

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
