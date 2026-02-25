#include "device_dem.h"
#define DIM 3
#define SMALL_NUM 1e-15


/*
============================================================
__device__ 関数群
============================================================
*/
__device__ __forceinline__
ContactCache d_calc_normal_force_wall(DeviceParticleGroup *p,int i,int j,Vec3 n,double delMag,double dist){

    int bi = i*DIM;
    int bj = j*DIM;

    /* get deltas */
    Vec3 del;
    del = vscalar(delMag,n);

    /* get relative velocity (for now assumed not moving)*/ 
    Vec3 v_rel;
    v_rel.x = p->v[bi+0];
    v_rel.y = p->v[bi+1];
    v_rel.z = p->v[bi+2];

    double v_reldotn = vdot(v_rel,n);

    /* get normal relative velocity*/
    ContactCache result;
    result.vn_rel = vscalar(v_reldotn,n);


    double eta = p->etaconst[i]*p->sqrtm[i];
    result.eta = eta;

    result.fn.x = p->k[i]*del.x - eta*result.vn_rel.x;
    result.fn.y = p->k[i]*del.y - eta*result.vn_rel.y;
    result.fn.z = p->k[i]*del.z - eta*result.vn_rel.z;


    /* calculate relative tangential velocity */
    Vec3 vrot;
    vrot.x = p->r[i]*p->angv[i*DIM+0];
    vrot.y = p->r[i]*p->angv[i*DIM+1];
    vrot.z = p->r[i]*p->angv[i*DIM+2];
    vrot = vcross(vrot,n);

    /* calculate relative tangential velocity */
    Vec3 vt;
    vt=vsub(v_rel,result.vn_rel);
    result.vt=vadd(vt,vrot);

    /* add force */

    p->f[bi+0] += result.fn.x;
    p->f[bi+1] += result.fn.y;
    p->f[bi+2] += result.fn.z;

    result.n = n;

    return result;
}

__device__ __forceinline__
void d_calc_tangential_force_wall(DeviceParticleGroup *p,int i,int j,ContactCache c){

    int bi = i*DIM;
    int ci = i*p->MAX_NEI;

    /* check history */

    int isInHis = 0;
    int neiInd=0;

    for (int k=0; k<p->numContWall[i]; k++){
        if (j == p->indHisWall[ci+k]){ /* if the contact particle is in the history */
            isInHis =1;
            neiInd=k;
            p->isContactWall[ci+neiInd]=1;
            break;
        }
    }

    if(isInHis == 0){
        neiInd=p->numContWall[i];
        p->isContactWall[ci+neiInd]=1;
        p->indHisWall[ci+neiInd]=j;
        p->numContWall[i] +=1;
    }


    double vtmag = vdot(c.vt,c.vt);
    vtmag = sqrt(vtmag);

    /* get the magnitude of delta t history */
    Vec3 delt_old;
    delt_old.x = p->deltHisxWall[ci+neiInd];
    delt_old.y = p->deltHisyWall[ci+neiInd];
    delt_old.z = p->deltHiszWall[ci+neiInd];

    double deltMag = vdot(delt_old,delt_old);
    deltMag = sqrt(deltMag);

    Vec3 t; /* tangential normal vector */

    if(vtmag>SMALL_NUM){
        t = vscalar(1./vtmag,c.vt);
    }else{
        t = vscalar(1./(deltMag+SMALL_NUM),delt_old); 
    }



    double dt = p->dt;

    /* get new delta t*/
    Vec3 delt_new;
    delt_new.x =deltMag*t.x + c.vt.x*dt;
    delt_new.y =deltMag*t.y + c.vt.y*dt;
    delt_new.z =deltMag*t.z + c.vt.z*dt;

    Vec3 ft;

    ft.x = p->k[i]*delt_new.x - c.eta*c.vt.x;
    ft.y = p->k[i]*delt_new.y - c.eta*c.vt.y;
    ft.z = p->k[i]*delt_new.z - c.eta*c.vt.z;

    /* ========== Friction ============ */
    double ftsq = vdot(ft,ft);
    double fnsq = vdot(c.fn,c.fn);

    if(ftsq>(p->mu*p->mu)*fnsq){ /* slip */
        double fnnorm = sqrt(fnsq);
        ft = vscalar(-p->mu*fnnorm,t);
        delt_new = delt_old;
    }

    /* ======== add angular acceleration ========== */
    Vec3 mom;
    mom = vcross(c.n,ft);
    mom = vscalar(p->r[i],mom);

    p->mom[bi+0] += mom.x;
    p->mom[bi+1] += mom.y;
    p->mom[bi+2] += mom.z;

    /* ======== add force ========== */
    p->f[bi+0] += ft.x;
    p->f[bi+1] += ft.y;
    p->f[bi+2] += ft.z;


    /* ====== refresh history ======= */
    p->deltHisxWall[ci+neiInd] = delt_new.x;
    p->deltHisyWall[ci+neiInd] = delt_new.y;
    p->deltHiszWall[ci+neiInd] = delt_new.z;
}


__device__ __forceinline__
void d_calc_tangential_force(DeviceParticleGroup *p,int i,int j,ContactCache c){

    int bi = i*DIM;
    int ci = i*p->MAX_NEI;

    /* check history */

    int isInHis = 0;
    int neiInd=0;
    int pId = p->pId[j];

    for (int k=0; k<p->numCont[i]; k++){
        if (pId == p->indHis[ci+k]){ /* if the contact particle is in the history */
            isInHis =1;
            neiInd=k;
            p->isContact[ci+neiInd]=1;
            break;
        }
    }

    if(isInHis == 0){
        neiInd=p->numCont[i];
    }

    p->isContact[ci+neiInd]=1;
    p->indHis[ci+neiInd]=pId;
    p->numCont[i] +=1;




    double vtmag = vdot(c.vt,c.vt);
    vtmag = sqrt(vtmag);

    /* get the magnitude of delta t history */
    Vec3 delt_old;
    delt_old.x = p->deltHisx[ci+neiInd];
    delt_old.y = p->deltHisy[ci+neiInd];
    delt_old.z = p->deltHisz[ci+neiInd];

    double deltMag = vdot(delt_old,delt_old);
    deltMag = sqrt(deltMag);

    Vec3 t; /* tangential normal vector */

    if(vtmag>SMALL_NUM){
        t = vscalar(1./vtmag,c.vt);
    }else{
        t = vscalar(1./(deltMag+SMALL_NUM),delt_old); 
    }



    double dt = p->dt;

    /* get new delta t*/
    Vec3 delt_new;
    delt_new.x =deltMag*t.x + c.vt.x*dt;
    delt_new.y =deltMag*t.y + c.vt.y*dt;
    delt_new.z =deltMag*t.z + c.vt.z*dt;

    Vec3 ft;

    ft.x = p->k[i]*delt_new.x - c.eta*c.vt.x;
    ft.y = p->k[i]*delt_new.y - c.eta*c.vt.y;
    ft.z = p->k[i]*delt_new.z - c.eta*c.vt.z;


    /* ========== Friction ============ */
    double ftsq = vdot(ft,ft);
    double fnsq = vdot(c.fn,c.fn);

    if(ftsq>(p->mu*p->mu)*fnsq){ /* slip */
        double fnnorm = sqrt(fnsq);
        ft = vscalar(-p->mu*fnnorm,t);
        delt_new = delt_old;
    }

    /* ======== add force ========== */
    p->f[bi+0] += ft.x;
    p->f[bi+1] += ft.y;
    p->f[bi+2] += ft.z;

    /* ======== add angular acceleration ========== */
    Vec3 mom;
    mom = vcross(c.n,ft);
    mom = vscalar(p->r[i],mom);

    p->mom[bi+0] += mom.x;
    p->mom[bi+1] += mom.y;
    p->mom[bi+2] += mom.z;

    /* ====== refresh history ======= */
    p->deltHisx[ci+neiInd] = delt_new.x;
    p->deltHisy[ci+neiInd] = delt_new.y;
    p->deltHisz[ci+neiInd] = delt_new.z;

}

__device__ __forceinline__
void d_update_history_wall(DeviceParticleGroup *p,int i){

    int ci = i*p->MAX_NEI;

    for(int k=0; k<p->numContWall[i]; k++){
        if (p->isContactWall[ci+k]==0){ 
            /* particle contact lost */
            int last = p->numContWall[i]-1;
            p->deltHisxWall[ci+k] = p->deltHisxWall[ci+last];
            p->deltHisyWall[ci+k] = p->deltHisyWall[ci+last];
            p->deltHiszWall[ci+k] = p->deltHiszWall[ci+last];

            p->deltHisxWall[ci+last]=0.;
            p->deltHisyWall[ci+last]=0.;
            p->deltHiszWall[ci+last]=0.;

            p->numContWall[i] -= 1;
            k -= 1; /* check if swapped particle is in contact */
        }
    }
}

__device__ __forceinline__
void d_update_history(DeviceParticleGroup *p,int i){

    int ci = i*p->MAX_NEI;
    for(int k=0; k<p->numCont[i]; k++){
        if (p->isContact[ci+k]==0){ 
            /* particle contact lost */
            int last = p->numCont[i]-1;
            p->deltHisx[ci+k] = p->deltHisx[ci+last];
            p->deltHisy[ci+k] = p->deltHisy[ci+last];
            p->deltHisz[ci+k] = p->deltHisz[ci+last];

            p->deltHisx[ci+last]=0.;
            p->deltHisy[ci+last]=0.;
            p->deltHisz[ci+last]=0.;

            p->numCont[i] -= 1;
            k -= 1; /* check if swapped particle is in contact */
        }
    }
}

__device__ __forceinline__
ContactCache d_calc_normal_force(DeviceParticleGroup* p,int i,int j,Vec3 n,double delMag,double dist){

    int bi = i*DIM;
    int bj = j*DIM;

    /* get deltas */
    Vec3 del;
    del = vscalar(delMag,n);

    /* get relative velocity */
    Vec3 v_rel;
    v_rel.x = p->v[bi+0] - p->v[bj+0];
    v_rel.y = p->v[bi+1] - p->v[bj+1];
    v_rel.z = p->v[bi+2] - p->v[bj+2];

    double v_reldotn = vdot(v_rel,n);

    /* get relative normal velocity*/
    ContactCache result;
    result.vn_rel = vscalar(v_reldotn,n);


    double m_eff = (p->m[i]*p->m[j])/(p->m[i]+p->m[j]);
    double eta = p->etaconst[i]*sqrt(m_eff);
    result.eta = eta;

    result.fn.x = p->k[i]*del.x - eta*result.vn_rel.x;
    result.fn.y = p->k[i]*del.y - eta*result.vn_rel.y;
    result.fn.z = p->k[i]*del.z - eta*result.vn_rel.z;

    /* calculate relative tangential velocity */
    Vec3 vrot;
    vrot.x = p->r[i]*p->angv[i*DIM+0]+p->r[j]*p->angv[j*DIM+0];
    vrot.y = p->r[i]*p->angv[i*DIM+1]+p->r[j]*p->angv[j*DIM+1];
    vrot.z = p->r[i]*p->angv[i*DIM+2]+p->r[j]*p->angv[j*DIM+2];
    vrot = vcross(vrot,n);

    Vec3 vt;
    vt=vsub(v_rel,result.vn_rel);
    result.vt=vadd(vt,vrot);

    /* add force */

    p->f[bi+0] += result.fn.x;
    p->f[bi+1] += result.fn.y;
    p->f[bi+2] += result.fn.z;

    result.n = n;

    return result;
}

__device__ __forceinline__
void d_wall_collision_naive(DeviceParticleGroup* ps,int i){
    //particle-wall
    for (int j=0; j<ps->walls.N; j++){
        double distsq = 0.;
        int bi = i*DIM;
        int bj = j*DIM;
        distsq += ps->walls.n[bj+0]*ps->x[bi+0];
        distsq += ps->walls.n[bj+1]*ps->x[bi+1];
        distsq += ps->walls.n[bj+2]*ps->x[bi+2];
        distsq += ps->walls.d[j];
        distsq = distsq*distsq;

        if (distsq < ps->rsq[i]){
            double dist = sqrt(distsq);
            double delmag = ps->r[i]-dist;
            if (delmag*ps->invr[i]*0.5>0.1){
                printf("overlap over 10%% with wall!!!!\n");
            }

            Vec3 n;
            n.x = ps->walls.n[bj+0];
            n.y = ps->walls.n[bj+1];
            n.z = ps->walls.n[bj+2];

            ContactCache c;
            c = d_calc_normal_force_wall(ps,i,j,n,delmag,dist);
            d_calc_tangential_force_wall(ps,i,j,c);

        }
    }
    d_update_history_wall(ps,i);
}

__device__ __forceinline__
void d_particle_collision_cell_linked_noVec3(DeviceParticleGroup* p, int i, DeviceBoundingBox* box){
    //particle-particle

    /* cycle through neighbor cells */
    int bi = i*DIM;
    int x=p->cellx[bi+0];
    int y=p->cellx[bi+1];
    int z=p->cellx[bi+2];

    for (int sx=-1; sx<=1; sx++){
        for (int sy=-1; sy<=1; sy++){
            for (int sz=-1; sz<=1; sz++){
                int cellId = (box->sizey*(z+sz)+y+sy)*box->sizex+x+sx;
                int start = box->pStart[cellId];
                int end = start+box->pNum[cellId];
                for (int k=start; k<end; k++){
                    int j = box->pList[k];
                    if (i==j){
                        continue;
                    }else{
                        int bj = j*DIM;
                        /* normal points toward particle i */
                        double dx = p->x[bi+0]- p->x[bj+0];
                        double dy = p->x[bi+1]- p->x[bj+1];
                        double dz = p->x[bi+2]- p->x[bj+2];
                        double distsq =dx*dx+dy*dy+dz*dz;
                        double R = p->r[i]+p->r[j];

                        if (distsq<R*R){
                            double dist = sqrt(distsq);
                            double delta = R-dist;
                            if(delta*p->invr[i]*0.5>0.05){
                                printf("overlap over 5%%!!!!\n");
                                printf("delta is %f %%, dist is %f\n",delta*p->invr[i]+0.5*100.,dist);
                            }

                            /* get normal direction */
                            double nx,ny,nz;
                            nx = dx/dist;
                            ny = dy/dist;
                            nz = dz/dist;

                            /* get deltas */
                            double delx,dely,delz;
                            delx = nx*delta;
                            dely = ny*delta;
                            delz = nz*delta;

                            /* get relative velocity */
                            double v_relx,v_rely,v_relz;
                            v_relx = p->v[bi+0] - p->v[bj+0];
                            v_rely = p->v[bi+1] - p->v[bj+1];
                            v_relz = p->v[bi+2] - p->v[bj+2];

                            double v_reldotn = v_relx*nx+v_rely*ny+v_relz*nz;

                            /* get normal relative velocity*/
                            double vn_relx,vn_rely,vn_relz;
                            vn_relx = v_reldotn*nx;
                            vn_rely = v_reldotn*ny;
                            vn_relz = v_reldotn*nz;

                            double m_eff = 1./(p->invm[i]+p->invm[j]);
                            double eta = p->etaconst[i]*sqrt(m_eff);

                            p->f[bi+0]+= p->k[i]*delx - eta*vn_relx;
                            p->f[bi+1] += p->k[i]*dely - eta*vn_rely;
                            p->f[bi+2] += p->k[i]*delz - eta*vn_relz;
                        }
                    }
                }

            }
        }
    }
}

__device__ __forceinline__
void d_particle_collision_cell_linked(DeviceParticleGroup* p, int i, DeviceBoundingBox* box){
    //particle-particle

    /* cycle through neighbor cells */
    int bi = i*DIM;
    int x=p->cellx[bi+0];
    int y=p->cellx[bi+1];
    int z=p->cellx[bi+2];

    for (int sx=-1; sx<=1; sx++){
        for (int sy=-1; sy<=1; sy++){
            for (int sz=-1; sz<=1; sz++){
                int cellId = (box->sizey*(z+sz)+y+sy)*box->sizex+x+sx;


                int start = box->pStart[cellId];
                int end = start+box->pNum[cellId];
                for (int k=start; k<end; k++){
                    int j = box->pList[k];
                    if (i==j){
                        continue;
                    }else{
                        int bj = j*DIM;

                        Vec3 del;
                        /* normal points toward particle i */
                        del.x = p->x[bi+0]- p->x[bj+0];
                        del.y = p->x[bi+1]- p->x[bj+1];
                        del.z = p->x[bi+2]- p->x[bj+2];
                        double distsq = vdot(del,del);
                        double R = p->r[i]+p->r[j];

                        if (distsq<R*R){
                            double dist = sqrt(distsq);
                            double delMag = R-dist;
                            if(delMag*p->invr[i]*0.5>0.05){
                                printf("overlap over 5%%!!!!\n");
                                printf("delta is %f %%, dist is %f\n",delMag*p->invr[i]+0.5*100.,dist);
                            }

                            /* ======================================================
                               Force Calculation
                               ======================================================*/
                            /* get normal direction */
                            Vec3 n;
                            n.x = del.x/dist;
                            n.y = del.y/dist;
                            n.z = del.z/dist;
                            ContactCache c;
                            c = d_calc_normal_force(p,i,j,n,delMag,dist);


                            d_calc_tangential_force(p,i,j,c);

                        }
                    }
                }
            }
        }
    }/* neighbor cell search done */

    /* update contact history */
    d_update_history(p,i);
}

__device__ __forceinline__
void particle_collision_naive(DeviceParticleGroup* ps, int i){
    //particle-particle
    for (int j=0; j<ps->N; j++){
        if (i==j){
            continue;
        }else{
            int bi = i*DIM;
            int bj = j*DIM;
            /* normal points toward particle i */
            double dx = ps->x[bi+0]- ps->x[bj+0];
            double dy = ps->x[bi+1]- ps->x[bj+1];
            double dz = ps->x[bi+2]- ps->x[bj+2];
            double distsq =dx*dx+dy*dy+dz*dz;
            double R = ps->r[i]+ps->r[j];

            if (distsq<R*R){
                double dist = sqrt(distsq);
                double delta = R-dist;

                /* get normal direction */
                double nx,ny,nz;
                nx = dx/dist;
                ny = dy/dist;
                nz = dz/dist;

                /* get deltas */
                double delx,dely,delz;
                delx = nx*delta;
                dely = ny*delta;
                delz = nz*delta;

                /* get relative velocity */
                double v_relx,v_rely,v_relz;
                v_relx = ps->v[bi+0] - ps->v[bj+0];
                v_rely = ps->v[bi+1] - ps->v[bj+1];
                v_relz = ps->v[bi+2] - ps->v[bj+2];

                double v_reldotn = v_relx*nx+v_rely*ny+v_relz*nz;

                /* get normal relative velocity*/
                double vn_relx,vn_rely,vn_relz;
                vn_relx = v_reldotn*nx;
                vn_rely = v_reldotn*ny;
                vn_relz = v_reldotn*nz;

                double m_eff = (ps->m[i]*ps->m[j])/(ps->m[i]+ps->m[j]);
                double eta = ps->etaconst[i]*sqrt(m_eff);

                ps->f[bi+0]+= ps->k[i]*delx - eta*vn_relx;
                ps->f[bi+1] += ps->k[i]*dely - eta*vn_rely;
                ps->f[bi+2] += ps->k[i]*delz - eta*vn_relz;


            }
        }
    }
}
    __device__ __forceinline__
void updateAcceleration(DeviceParticleGroup* p,
        int i
        )
{
    int bi = i*DIM;
    p->a[bi+0] = (p->g[0]+p->f[bi+0]*p->invm[i]) ;
    p->a[bi+1] = (p->g[1]+p->f[bi+1]*p->invm[i]);
    p->a[bi+2] = (p->g[2]+p->f[bi+2]*p->invm[i]);
}

/* 速度更新（重力適用） */
    __device__ __forceinline__
void updateVelocity(DeviceParticleGroup* p,
        int i)
{
    int bi = i*DIM;
    p->v[bi+0] += p->a[bi+0] * p->dt;
    p->v[bi+1] += p->a[bi+1] * p->dt;
    p->v[bi+2] += p->a[bi+2] * p->dt;
}


/* 位置更新（オイラー法） */
    __device__ __forceinline__
void updatePosition(DeviceParticleGroup* p,
        int i)
{
    int bi = i*DIM;
    p->x[bi+0] += p->v[bi+0] * p->dt;
    p->x[bi+1] += p->v[bi+1] * p->dt;
    p->x[bi+2] += p->v[bi+2] * p->dt;
}




/*
   ============================================================
   カーネル
   ============================================================
 */
__global__ void check_g_kernel(DeviceParticleGroup* ps){
    printf("g in kernel: %f %f %f\n",ps->g[0],ps->g[1],ps->g[2]);
    printf("wall 0 in kernel: %f %f %f\n",ps->walls.n[0],ps->walls.n[1],ps->walls.n[2]);
    printf("wall 1 in kernel: %f %f %f\n",ps->walls.n[1*DIM+0],ps->walls.n[1*DIM+1],ps->walls.n[1*DIM+2]);
    printf("value in kernel: %f %f %f\n",ps->r[0],ps->m[0],ps->invm[0]);
    printf("value in kernel: dt = %f, etaconst = %f\n",ps->dt,ps->etaconst[0]);
}

__global__ void k_collision(DeviceParticleGroup* p, DeviceBoundingBox* box){
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= p->N) return;

    d_particle_collision_cell_linked(p,i,box);
    // d_particle_collision_cell_linked_noVec3(p,i,box);
    d_wall_collision_naive(p,i);

}

__global__ void integrateKernel(DeviceParticleGroup* p){
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= p->N) return;

    int bi = i*DIM;
    /* initialize */
    p->f[bi+0]=0.;
    p->f[bi+1]=0.;
    p->f[bi+2]=0.;

    particle_collision_naive(p,i);
    d_wall_collision_naive(p,i);
    updateAcceleration(p,i);
    // 処理A
    updateVelocity(p, i);

    // 処理B
    updatePosition(p, i);

}

__global__ void k_integrate(DeviceParticleGroup* p){
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= p->N) return;

    updateAcceleration(p,i);
    // 処理A
    updateVelocity(p, i);

    // 処理B
    updatePosition(p, i);

    int bi = i*DIM;


    // angular acceleration
    p->anga[bi+0] = p->mom[bi+0]*p->invmoi[i];
    p->anga[bi+1] = p->mom[bi+1]*p->invmoi[i];
    p->anga[bi+2] = p->mom[bi+2]*p->invmoi[i];

    // angular velocity
    p->angv[bi+0] += p->anga[bi+0] * p->dt;
    p->angv[bi+1] += p->anga[bi+1] * p->dt;
    p->angv[bi+2] += p->anga[bi+2] * p->dt;


}


/*
   ======================================================
   main routine 
   ======================================================
 */

void device_dem(ParticleSystem *p, BoundingBox *box, int gridSize, int blockSize){

    /* initialize force */
    cudaMemset(p->d_group.f, 0, sizeof(double)*DIM*p->d_group.N);
    cudaMemset(p->d_group.mom, 0, sizeof(double)*DIM*p->d_group.N);
    cudaMemset(p->d_group.isContact, 0, sizeof(int)*p->d_group.N*p->d_group.MAX_NEI);
    cudaMemset(p->d_group.isContactWall, 0, sizeof(int)*p->d_group.N*p->d_group.MAX_NEI);

    d_update_pList(p,box,gridSize, blockSize);
    k_collision<<<gridSize,blockSize>>>(p->d_groupPtr,box->d_boxPtr);
    k_integrate<<<gridSize, blockSize>>>(p->d_groupPtr);
    // cudaError_t err = cudaGetLastError();
    //printf("CUDA error = %s\n", cudaGetErrorString(err));
    //    cudaDeviceSynchronize();

}

