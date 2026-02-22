#include "cpu_dem.h"

/*
============================================================
cpu functions
============================================================
*/

void wall_collision_naive(ParticleSystem* ps){
    for (int i=0; i<ps->N; i++){
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
                double delta = ps->r[i]-dist;
                if (delta*ps->invr[i]>0.01){
                    printf("overlap over 10% with wall!!!!\n");
                }


                /* get deltas */
                double delx,dely,delz;
                delx = ps->walls.n[bj+0]*delta;
                dely = ps->walls.n[bj+1]*delta;
                delz = ps->walls.n[bj+2]*delta;

                /* get relative velocity (here we assumed walls are stationary) */
                double v_relx,v_rely,v_relz;
                v_relx = ps->v[bi+0];
                v_rely = ps->v[bi+1];
                v_relz = ps->v[bi+2];

                double v_reldotn = v_relx*ps->walls.n[bj+0]+v_rely*ps->walls.n[bj+1]+v_relz*ps->walls.n[bj+2];

                /* get normal relative velocity*/
                double vn_relx,vn_rely,vn_relz;
                vn_relx = v_reldotn*ps->walls.n[bj+0];
                vn_rely = v_reldotn*ps->walls.n[bj+1];
                vn_relz = v_reldotn*ps->walls.n[bj+2];

                double eta = ps->etaconst[i]*ps->sqrtm[i];

                ps->f[bi+0] += ps->k[i]*delx - eta*vn_relx;
                ps->f[bi+1] += ps->k[i]*dely - eta*vn_rely;
                ps->f[bi+2] += ps->k[i]*delz - eta*vn_relz;

            }
        }
    }
}

inline ContactCache calc_normal_force(ParticleSystem *p,int i,int j,Vec3 n,double delMag,double dist){

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

    /* get normal relative velocity*/
    ContactCache result;
    result.vn_rel = vscalar(v_reldotn,n);


    double m_eff = (p->m[i]*p->m[j])/(p->m[i]+p->m[j]);
    double eta = p->etaconst[i]*sqrt(m_eff);

    result.fn.x = p->k[i]*del.x - eta*result.vn_rel.x;
    result.fn.y = p->k[i]*del.y - eta*result.vn_rel.y;
    result.fn.z = p->k[i]*del.z - eta*result.vn_rel.z;

    return result;
}

void particle_collision_cell_linked(ParticleSystem* p, BoundingBox *box){
    update_pList(p,box);

    for (int i=0; i<p->N; i++){
        //particle-particle
        /* cycle through neighbor cells */
        int bi=i*DIM;
        int x=p->cellx[bi+0];
        int y=p->cellx[bi+1];
        int z=p->cellx[bi+2];

        for (int sx=-1; sx<=1; sx++){
            for (int sy=-1; sy<=1; sy++){
                for (int sz=-1; sz<=1; sz++){
                    int cellId = (box->sizey*(z+sz)+y+sy)*box->sizex+x+sx;

                    int start = box->pStart[cellId];
                    int end = start+box->pNum[cellId];
                    for (int k=box->pStart[cellId]; k<end; k++){
                        int j = box->pList[k];
                        if (i==j){
                            continue;
                        }else{
                            int bj=j*DIM;

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
                                if (delMag*p->invr[i]*0.5>0.01){
                                    printf("overlap over 10%!!!!\n");
                                    printf("delta is %f, dist is %f\n",delMag,dist);
                                }
                                /* ======================================================
                                   Normal Force Calculation
                                   ======================================================*/
                                /* get normal direction */
                                Vec3 n;
                                n.x = del.x/dist;
                                n.y = del.y/dist;
                                n.z = del.z/dist;
                                ContactCache c;
                                c = calc_normal_force(p,i,j,n,delMag,dist);
                                p->f[bi+0] += c.fn.x;
                                p->f[bi+1] += c.fn.y;
                                p->f[bi+2] += c.fn.z;

                            }
                        }
                    }
                }
            }
        }
    }
}

void particle_collision_cell_linked_noVec3(ParticleSystem* ps, BoundingBox *box){
    update_pList(ps,box);

    for (int i=0; i<ps->N; i++){
        //particle-particle
        /* cycle through neighbor cells */
        int bi=i*DIM;
        int x=ps->cellx[bi+0];
        int y=ps->cellx[bi+1];
        int z=ps->cellx[bi+2];

        for (int sx=-1; sx<=1; sx++){
            for (int sy=-1; sy<=1; sy++){
                for (int sz=-1; sz<=1; sz++){
                    int cellId = (box->sizey*(z+sz)+y+sy)*box->sizex+x+sx;

                    int start = box->pStart[cellId];
                    int end = start+box->pNum[cellId];
                    for (int k=box->pStart[cellId]; k<end; k++){
                        int j = box->pList[k];
                        if (i==j){
                            continue;
                        }else{
                            int bj=j*DIM;
                            /* normal points toward particle i */
                            double dx = ps->x[bi+0]- ps->x[bj+0];
                            double dy = ps->x[bi+1]- ps->x[bj+1];
                            double dz = ps->x[bi+2]- ps->x[bj+2];
                            double distsq =dx*dx+dy*dy+dz*dz;
                            double R = ps->r[i]+ps->r[j];

                            if (distsq<R*R){
                                double dist = sqrt(distsq);
                                double delta = R-dist;
                                if (delta*ps->invr[i]*0.5>0.01){
                                    printf("overlap over 10%!!!!\n");
                                    printf("delta is %f, dist is %f\n",delta,dist);
                                }
                                /* ======================================================
                                   Normal Force Calculation
                                   ======================================================*/
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

                            }
                        }
                    }
                }
            }
        }
    }
}

void particle_collision_naive(ParticleSystem* ps){
    for (int i=0; i<ps->N; i++){
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
                    if (delta*ps->invr[i]*0.5>0.01){
                        printf("overlap over 10%!!!!\n");
                        printf("delta is %f, dist is %f\n",delta,dist);
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
}

void integrateCPU(ParticleSystem* ps, BoundingBox* box)
{
    /* initialize */
    for (int i=0; i<ps->N*DIM; i++){
        ps->f[i]=0.;
    }

    //particle_collision_naive(ps);
    particle_collision_cell_linked(ps,box);
    //particle_collision_cell_linked_noVec3(ps,box);
    wall_collision_naive(ps);

    for (int i = 0; i < ps->N; i++)
    {
        int bi = i*DIM;
        // acceleration
        ps->a[bi+0] = ps->f[bi+0]*ps->invm[i]+ps->g[0];
        ps->a[bi+1] = ps->f[bi+1]*ps->invm[i]+ps->g[1];
        ps->a[bi+2] = ps->f[bi+2]*ps->invm[i]+ps->g[2];

        // velocity
        ps->v[bi+0] += ps->a[bi+0] * ps->dt;
        ps->v[bi+1] += ps->a[bi+1] * ps->dt;
        ps->v[bi+2] += ps->a[bi+2] * ps->dt;

        // position
        ps->x[bi+0] += ps->v[bi+0] * ps->dt;
        ps->x[bi+1] += ps->v[bi+1] * ps->dt;
        ps->x[bi+2] += ps->v[bi+2] * ps->dt;

    }
}
