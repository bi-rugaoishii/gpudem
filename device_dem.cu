#include "device_dem.h"
#define DIM 3


/*
============================================================
__device__ 関数群
============================================================
*/

__device__ __forceinline__
void wall_collision_naive(DeviceParticleGroup* ps,int i){
    //particle-wall
    for (int j=0; j<ps->walls.N; j++){
        double distsq = 0.;
        distsq += ps->walls.n[j*DIM+0]*ps->x[i*DIM+0];
        distsq += ps->walls.n[j*DIM+1]*ps->x[i*DIM+1];
        distsq += ps->walls.n[j*DIM+2]*ps->x[i*DIM+2];
        distsq += ps->walls.d[j];
        distsq = distsq*distsq;
        
        if (distsq < ps->rsq[i]){
            double dist = sqrt(distsq);
            double delta = ps->r[i]-dist;
                if (delta*ps->invr[i]>0.01){
                    printf("overlap over 10% with wall!!!!\n");
                }


            /* get deltas */
            double delx,dely,delz;
            delx = ps->walls.n[j*DIM+0]*delta;
            dely = ps->walls.n[j*DIM+1]*delta;
            delz = ps->walls.n[j*DIM+2]*delta;

            /* get relative velocity (here we assumed walls are stationary) */
            double v_relx,v_rely,v_relz;
            v_relx = ps->v[i*DIM+0];
            v_rely = ps->v[i*DIM+1];
            v_relz = ps->v[i*DIM+2];

            double v_reldotn = v_relx*ps->walls.n[j*DIM+0]+v_rely*ps->walls.n[j*DIM+1]+v_relz*ps->walls.n[j*DIM+2];

            /* get normal relative velocity*/
            double vn_relx,vn_rely,vn_relz;
            vn_relx = v_reldotn*ps->walls.n[j*DIM+0];
            vn_rely = v_reldotn*ps->walls.n[j*DIM+1];
            vn_relz = v_reldotn*ps->walls.n[j*DIM+2];

            double eta = ps->etaconst[i]*ps->sqrtm[i];

            
            ps->f[i*DIM+0] += ps->k[i]*delx - eta*vn_relx;
            ps->f[i*DIM+1] += ps->k[i]*dely - eta*vn_rely;
            ps->f[i*DIM+2] += ps->k[i]*delz - eta*vn_relz;
            
            /*
            ps->f[i*DIM+0] += ps->k[i]*delx;
            ps->f[i*DIM+1] += ps->k[i]*dely;
            ps->f[i*DIM+2] += ps->k[i]*delz;
            */

        }
    }
}

__device__ __forceinline__
void particle_collision_naive(DeviceParticleGroup* ps, int i){
    //particle-particle
    for (int j=0; j<ps->N; j++){
        if (i==j){
            continue;
        }else{
            /* normal points toward particle i */
            double dx = ps->x[i*DIM+0]- ps->x[j*DIM+0];
            double dy = ps->x[i*DIM+1]- ps->x[j*DIM+1];
            double dz = ps->x[i*DIM+2]- ps->x[j*DIM+2];
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
                v_relx = ps->v[i*DIM+0] - ps->v[j*DIM+0];
                v_rely = ps->v[i*DIM+1] - ps->v[j*DIM+1];
                v_relz = ps->v[i*DIM+2] - ps->v[j*DIM+2];

                double v_reldotn = v_relx*nx+v_rely*ny+v_relz*nz;

                /* get normal relative velocity*/
                double vn_relx,vn_rely,vn_relz;
                vn_relx = v_reldotn*nx;
                vn_rely = v_reldotn*ny;
                vn_relz = v_reldotn*nz;

                double m_eff = (ps->m[i]*ps->m[j])/(ps->m[i]+ps->m[j]);
                double eta = ps->etaconst[i]*sqrt(m_eff);

                ps->f[i*DIM+0]+= ps->k[i]*delx - eta*vn_relx;
                ps->f[i*DIM+1] += ps->k[i]*dely - eta*vn_rely;
                ps->f[i*DIM+2] += ps->k[i]*delz - eta*vn_relz;

                /*
                   ps->f[i*DIM+0]+= ps->k[i]*delx;
                   ps->f[i*DIM+1]+= ps->k[i]*dely;
                   ps->f[i*DIM+2]+= ps->k[i]*delz;
                 */

            }
        }
    }
}
__device__ __forceinline__
void updateAcceleration(DeviceParticleGroup* p,
        int i
        )
{
    p->a[i*DIM+0] = (p->g[0]+p->f[i*DIM+0]*p->invm[i]) ;
    p->a[i*DIM+1] = (p->g[1]+p->f[i*DIM+1]*p->invm[i]);
    p->a[i*DIM+2] = (p->g[2]+p->f[i*DIM+2]*p->invm[i]);
}
/* 速度更新（重力適用） */
    __device__ __forceinline__
void updateVelocity(DeviceParticleGroup* p,
        int i)
{
    p->v[i*DIM+0] += p->a[i*DIM+0] * p->dt;
    p->v[i*DIM+1] += p->a[i*DIM+1] * p->dt;
    p->v[i*DIM+2] += p->a[i*DIM+2] * p->dt;
}


/* 位置更新（オイラー法） */
    __device__ __forceinline__
void updatePosition(DeviceParticleGroup* p,
        int i)
{
    p->x[i*DIM+0] += p->v[i*DIM+0] * p->dt;
    p->x[i*DIM+1] += p->v[i*DIM+1] * p->dt;
    p->x[i*DIM+2] += p->v[i*DIM+2] * p->dt;
}




/*
   ============================================================
   カーネル
   ============================================================
 */
__global__ void check_g_kernel(DeviceParticleGroup ps){
    printf("g in kernel: %f %f %f\n",ps.g[0],ps.g[1],ps.g[2]);
    printf("wall 0 in kernel: %f %f %f\n",ps.walls.n[0],ps.walls.n[1],ps.walls.n[2]);
    printf("wall 1 in kernel: %f %f %f\n",ps.walls.n[1*DIM+0],ps.walls.n[1*DIM+1],ps.walls.n[2*DIM+2]);
    printf("value in kernel: %f %f %f\n",ps.r[0],ps.m[0],ps.invm[0]);
    printf("value in kernel: dt = %f, etaconst = %f\n",ps.dt,ps.etaconst[0]);
}
__global__ void integrateKernel(DeviceParticleGroup p)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= p.N) return;

    /* initialize */
    p.f[i*DIM+0]=0.;
    p.f[i*DIM+1]=0.;
    p.f[i*DIM+2]=0.;

    particle_collision_naive(&p,i);
    wall_collision_naive(&p,i);
    updateAcceleration(&p,i);
    // 処理A
    updateVelocity(&p, i);

    // 処理B
    updatePosition(&p, i);

}


