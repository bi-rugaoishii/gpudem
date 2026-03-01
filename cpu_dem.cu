#include "cpu_dem.h"
#define SMALL_NUM 1e-15

/*
============================================================
cpu functions
============================================================
*/

/* === triangle related === */

TriangleContactCache dist_triangle(ParticleSystem* ps, int i, TriangleMesh* mesh, int j){
    /* check with plane first */
    double dist = 0.;
    int bi = i*DIM;
    TriangleContactCache res;
    Vec3 ntri;
    ntri.x = mesh->nx[j];
    ntri.y = mesh->ny[j];
    ntri.z = mesh->nz[j];

    res.n.x=0.;
    res.n.y=0.;
    res.n.z=0.;
    res.dist=1.e10;

    Vec3 pos;
    pos.x = ps->x[bi+0];
    pos.y = ps->x[bi+1];
    pos.z = ps->x[bi+2];

    dist += ntri.x*pos.x;
    dist += ntri.y*pos.y;
    dist += ntri.z*pos.z;
    dist += mesh->d[j];
    double distsq = dist*dist;
    if(distsq > ps->rsq[i]){
        /* no collision */
        return res;
    }

    /* check if inside the triangle */

    /* get vertex indices */
    int v0 = mesh->tri_i0[j];
    int v1 = mesh->tri_i1[j];
    int v2 = mesh->tri_i2[j];

    Vec3 ap;
    ap.x = pos.x - mesh->mx[v0];
    ap.y = pos.y - mesh->my[v0];
    ap.z = pos.z - mesh->mz[v0];

    double d00 = mesh->d00[j];
    double d01 = mesh->d01[j];
    double d11 = mesh->d11[j];
    double denom = mesh->denom[j];

    double d20 = ap.x*mesh->e01x[j] + ap.y*mesh->e01y[j] + ap.z*mesh->e01z[j];
    double d21 = ap.x*mesh->e02x[j] + ap.y*mesh->e02y[j] + ap.z*mesh->e02z[j];
    
    double v = (d11*d20 - d01*d21)*denom;
    double w = (d00*d21 - d01*d20)*denom;
    double u = 1.0 - v - w;
    /* 1) face inside */
    if (u>0. &&  v >0. && w>0.){
        double sign= dist<0 ? -1.0 : 1.0;
        res.n = vscalar(sign,ntri);
        res.dist = sign*dist;
        res.hitAt = -1;
        return res;
    }

    /* ========== vertices ========== */

    if (v<=0. && w<=0.0){
        Vec3 tx;
        tx.x = mesh->mx[v0];
        tx.y = mesh->my[v0];
        tx.z = mesh->mz[v0];

        res.n = vsub(pos,tx);
       distsq = vdot(res.n,res.n);
       if(distsq < ps->rsq[i]){
           res.dist = sqrt(distsq);
           res.n = vscalar(1./res.dist,res.n);
       }
       res.hitAt = v0; 
       return res;
    }

    if (u<=0. && w<=0.0){
        Vec3 tx;
        tx.x = mesh->mx[v1];
        tx.y = mesh->my[v1];
        tx.z = mesh->mz[v1];

        res.n = vsub(pos,tx);
        distsq = vdot(res.n,res.n);
       if(distsq < ps->rsq[i]){
           res.dist = sqrt(distsq);
           res.n = vscalar(1./res.dist,res.n);
       }

       res.hitAt = v1; 
       return res;
    }

    if (u<=0. && v<=0.0){
        Vec3 tx;
        tx.x = mesh->mx[v2];
        tx.y = mesh->my[v2];
        tx.z = mesh->mz[v2];

        res.n = vsub(pos,tx);
        distsq = vdot(res.n,res.n);
        if(distsq < ps->rsq[i]){
            res.dist = sqrt(distsq); 
            res.n = vscalar(1./res.dist,res.n);
        }

        res.hitAt = v2; 
        return res;
    }

    /* ========== edge ========== */
    if (u<=0.0){
        Vec3 bp;

        Vec3 tmpv;
        tmpv.x = mesh->mx[v1];
        tmpv.y = mesh->my[v1];
        tmpv.z = mesh->mz[v1];

        Vec3 tmpE;
        tmpE.x = mesh->e12x[j];
        tmpE.y = mesh->e12y[j];
        tmpE.z = mesh->e12z[j];

        bp = vsub(pos,tmpv);

        double t = vdot(bp,tmpE)*mesh->d22inv[j];
        t = fmax(0., fmin(1.0,t));

        Vec3 hitPoint = vadd(tmpv,vscalar(t,tmpE));
        res.n = vsub(pos,hitPoint);
        distsq = vdot(res.n,res.n);
        if(distsq < ps->rsq[i]){
            res.dist = sqrt(distsq);
            res.n = vscalar(1./res.dist,res.n);
        }
        res.hitAt = mesh->tri_e1[j]; 
        return res;
    }

    if (v<=0.0){
        Vec3 tmpE;
        tmpE.x = mesh->e02x[j];
        tmpE.y = mesh->e02y[j];
        tmpE.z = mesh->e02z[j];

        double t = d21*mesh->d11inv[j];
        t = fmax(0., fmin(1.0,t));
        Vec3 tmpv;
        tmpv.x = mesh->mx[v0];
        tmpv.y = mesh->my[v0];
        tmpv.z = mesh->mz[v0];

        Vec3 hitPoint = vadd(tmpv,vscalar(t,tmpE));
        res.n = vsub(pos,hitPoint);
        distsq = vdot(res.n,res.n);

        if(distsq < ps->rsq[i]){
            res.dist = sqrt(distsq);
            res.n = vscalar(1./res.dist,res.n);
        }
        res.hitAt = mesh->tri_e2[j]; 
        return res;
    }

    /* w=<0.0 or 0.,0.,0. */
    Vec3 tmpE;
    tmpE.x = mesh->e01x[j];
    tmpE.y = mesh->e01y[j];
    tmpE.z = mesh->e01z[j];


    double t = d20*mesh->d00inv[j];
    t = fmax(0., fmin(1.0,t));
    Vec3 tmpv;
    tmpv.x = mesh->mx[v0];
    tmpv.y = mesh->my[v0];
    tmpv.z = mesh->mz[v0];

    Vec3 hitPoint = vadd(tmpv,vscalar(t,tmpE));
    res.n = vsub(pos,hitPoint);
    distsq = vdot(res.n,res.n);
    if(distsq < ps->rsq[i]){
        res.dist = sqrt(distsq);
        res.n = vscalar(1./res.dist,res.n);
    }
    res.hitAt = mesh->tri_e0[j]; 
    return res;
}

void wall_collision_triangles(ParticleSystem* p,BoundingBox *box, TriangleMesh* mesh){
    for (int i=0; i<p->N; i++){
        if(p->isActive[i]!=1){
            continue;
        }
        //particle-particle
        /* cycle through neighbor cells */
        int bi=i*DIM;
        int x=p->cellx[bi+0];
        int y=p->cellx[bi+1];
        int z=p->cellx[bi+2];

        int numContVorENow = 0;
        int numContWallNow = 0;
        for (int sx=-1; sx<=1; sx++){
            for (int sy=-1; sy<=1; sy++){
                for (int sz=-1; sz<=1; sz++){
                    int cellId = (box->sizey*(z+sz)+y+sy)*box->sizex+x+sx;
                    for(int j=0; j<box->tNum[cellId]; j++){
                        int indTri = box->tList[cellId*box->MAX_TRI+j];

                        int wasHitBefore = 0;
                        /* check if the triangle was already hit */
                        for (int k=0; k<numContWallNow; k++){
                            if (indTri == p->indHisWallNow[i*p->MAX_NEI+k]){
                                wasHitBefore =1;
                                continue;
                            }
                        }
                        if (wasHitBefore == 1){
                            continue;
                        }


                        TriangleContactCache tc;
                        tc = dist_triangle(p,i,mesh,indTri); 

                        if(tc.dist<p->r[i]){
                            double delmag;
                            if(tc.hitAt==-1){ //hit at face
                                delmag = p->r[i]-tc.dist;
                            }else{ //hit as face or edge
                                int end = numContVorENow;
                                int hadDuplicate=0;

                                for (int k=0; k<end; k++){
                                    if(p->indHisVorENow[i*p->MAX_NEI+k]==tc.hitAt){
                                        printf("duplicate collision of vertex or edge!!\n");
                                        hadDuplicate=1;
                                        break;
                                    }
                                }

                                if(hadDuplicate ==1){
                                    continue;

                                }

                                delmag = p->r[i]-tc.dist;


                                p->indHisVorENow[i*p->MAX_NEI+end] = tc.hitAt;
                                numContVorENow+=1;
                            }
                            if (delmag*p->invr[i]*0.5>0.1){
                                printf("overlap over 10%% with wall!!!!\n");
                            }

                            ContactCache c;
                            c = calc_normal_force_wall(p,i,j,tc.n,delmag);
                            calc_tangential_force_wall(p,i,j,c);
                            int numWall = numContWallNow;
                            p->indHisWallNow[i*p->MAX_NEI+numWall] = indTri;
                            numContWallNow+=1;
                        }

                    }
                }
            }
        }
        update_history_wall(p,i);

    }
}

inline ContactCache calc_normal_force_wall(ParticleSystem *p,int i,int j,Vec3 n,double delMag){

    int bi = i*DIM;

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

void wall_collision_naive(ParticleSystem* ps){
    for (int i=0; i<ps->N; i++){
        if(ps->isActive[i]!=1){
            continue;
        }
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
                c = calc_normal_force_wall(ps,i,j,n,delmag);
                calc_tangential_force_wall(ps,i,j,c);

            }
        }
        update_history_wall(ps,i);
    }


}

inline void update_history(ParticleSystem *p,int i){

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

inline void update_history_wall(ParticleSystem *p,int i){

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

inline void calc_tangential_force_wall(ParticleSystem *p,int i,int j,ContactCache c){

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
        if(p->numContWall[i] > p->MAX_NEI){
            printf("OH MY OH MY\n");
        }
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


inline void calc_tangential_force(ParticleSystem *p,int i,int j,ContactCache c){

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
        p->isContact[ci+neiInd]=1;
        p->indHis[ci+neiInd]=pId;
        p->numCont[i] +=1;
    }


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

inline ContactCache calc_normal_force(ParticleSystem *p,int i,int j,Vec3 n,double delMag){

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


    double m_eff = 1./(p->invm[i]+p->invm[j]);
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

void particle_collision_cell_linked_withSort_fastUpdate(ParticleSystem* p,ParticleSystem* tmpPs, BoundingBox *box){
    update_pList_withSort_fast(p, tmpPs, box);
    for (int i=0; i<p->N; i++){
        if(p->isActive[i]!=1){
            continue;
        }
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
                                if (delMag*p->invr[i]*0.5>0.05){
                                    printf("overlap over 5%% with pair %d %d!!!!\n",i,j);
                                    printf("overlap is %f %% with pair %d %d\n",delMag*p->invr[i]*0.5*100.,i,j);
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
                                c = calc_normal_force(p,i,j,n,delMag);


                                calc_tangential_force(p,i,j,c);

                            }
                        }
                    }
                }
            }
        }/* neighbor cell search done */

        /* update contact history */
        update_history(p,i);

    }
}

void particle_collision_verlet(ParticleSystem* p, BoundingBox *box){
    for (int i=0; i<p->N; i++){
        if(p->isActive[i]!=1){
            continue;
        }
        //particle-particle
        /* cycle through neighbor list */

        int bi=i*DIM;
        int end = p->numNei[i];
        int ci = i*p->MAX_NEI;
        for (int k=0; k<end; k++){
            int j = p->neiList[ci+k];
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
                if (delMag*p->invr[i]*0.5>0.05){
                    printf("overlap over 5%% with pair %d %d!!!!\n",i,j);
                    printf("overlap is %f %% with pair %d %d\n",delMag*p->invr[i]*0.5*100.,i,j);
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
                c = calc_normal_force(p,i,j,n,delMag);


                 calc_tangential_force(p,i,j,c);

            }
        }/* neighbor search done */

        /* update contact history */
        update_history(p,i);
    }
}

void particle_collision_cell_linked_withSort(ParticleSystem* p,ParticleSystem* tmpPs, BoundingBox *box){
    update_pList_withSort(p, tmpPs, box);
    for (int i=0; i<p->N; i++){
        if(p->isActive[i]!=1){
            continue;
        }
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
                                if (delMag*p->invr[i]*0.5>0.05){
                                    printf("overlap over 5%% with pair %d %d!!!!\n",i,j);
                                    printf("overlap is %f %% with pair %d %d\n",delMag*p->invr[i]*0.5*100.,i,j);
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
                                c = calc_normal_force(p,i,j,n,delMag);


                                calc_tangential_force(p,i,j,c);

                            }
                        }
                    }
                }
            }
        }/* neighbor cell search done */

        /* update contact history */
        update_history(p,i);

    }
}

void particle_collision_cell_linked_fastUpdate(ParticleSystem* p, BoundingBox *box){
    //update_pList(p,box);
    update_pList_fast(p,box);

    for (int i=0; i<p->N; i++){
        if(p->isActive[i]!=1){
            continue;
        }
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
                                if (delMag*p->invr[i]*0.5>0.05){
                                    printf("overlap over 5%%!!!!\n");
                                    printf("overlap is %f %%\n",delMag*p->invr[i]*0.5*100.);
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
                                c = calc_normal_force(p,i,j,n,delMag);


                                calc_tangential_force(p,i,j,c);

                            }
                        }
                    }
                }
            }
        }/* neighbor cell search done */

        /* update contact history */
        update_history(p,i);

    }
}

void particle_collision_cell_linked(ParticleSystem* p, BoundingBox *box){
    update_pList(p,box);

    for (int i=0; i<p->N; i++){
        if(p->isActive[i]!=1){
            continue;
        }
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
                                if (delMag*p->invr[i]*0.5>0.05){
                                    printf("overlap over %%!!!!\n");
                                    printf("overlap is %f %%\n",delMag*p->invr[i]*0.5*100.);
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
                                c = calc_normal_force(p,i,j,n,delMag);


                                calc_tangential_force(p,i,j,c);

                            }
                        }
                    }
                }
            }
        }/* neighbor cell search done */

        /* update contact history */
        update_history(p,i);

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
                                    printf("overlap over 10%%!!!!\n");
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

void particle_collision_naive(ParticleSystem* p){
    for (int i=0; i<p->N; i++){
        if(p->isActive[i]!=1){
            continue;
        }
        //particle-particle
        int bi=i*DIM;
        for (int j=0; j<p->N; j++){
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
                    if (delMag*p->invr[i]*0.5>0.05){
                        printf("overlap over 5%%!!!!\n");
                        printf("overlap is %f %%\n",delMag*p->invr[i]*0.5*100.);
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
                    c = calc_normal_force(p,i,j,n,delMag);
                    calc_tangential_force(p,i,j,c);

                }
            }
        }
        update_history(p,i);
    }
}

/* =========== verlet list related =============== */

int shouldRefresh(ParticleSystem *p, BoundingBox* box){
    int flag = 0;
    double threshSq = box->refreshThreshSq;
    for(int i=0; i<p->N; i++){
        if(p->isActive[i]!=1){
            continue;
        }
        int bi = i*DIM;
        Vec3 del;

        del.x = p->x[bi+0]- p->refx[i];
        del.y = p->x[bi+1]- p->refy[i];
        del.z = p->x[bi+2]- p->refz[i];

        double distsq=vdot(del,del);
        if(distsq>threshSq){
            flag=1;
            break;
        }
    }
    return flag;
}

/* ============== check Out of Bounds ============  */

void checkOoB(ParticleSystem *p, BoundingBox* box){
    for(int i=0; i<p->N; i++){
        double x=p->x[i*DIM+0];
        double y=p->x[i*DIM+1];
        double z=p->x[i*DIM+2];
        if( x<box->minx || x>box->maxx ||y<box->miny || y>box->maxy || z<box->minz || z>box->maxz){
            p->isActive[i]=0;
            continue;
        }       
    }
}

/* ============== dem mains ================= */

void cpu_dem_verlet_triangles(ParticleSystem* p, ParticleSystem *tmpP, BoundingBox* box,TriangleMesh *mesh, int step){
    /* initialize */
    for (int i=0; i<p->N*DIM; i++){
        p->f[i]=0.;
    }

    for (int i=0; i<p->N*DIM; i++){
        p->mom[i]=0.;
    }

    for (int i=0; i<p->N*p->MAX_NEI; i++){
        p->isContact[i]=0;
        p->isContactWall[i]=0;
    }



    particle_collision_verlet(p,box);

    wall_collision_triangles(p,box,mesh);

    /* update */
    for (int i = 0; i < p->N; i++)
    {
        if (p->isActive[i]!=1){
            continue;
        }
        int bi = i*DIM;
        // acceleration
        p->a[bi+0] = p->f[bi+0]*p->invm[i]+p->g[0];
        p->a[bi+1] = p->f[bi+1]*p->invm[i]+p->g[1];
        p->a[bi+2] = p->f[bi+2]*p->invm[i]+p->g[2];

        // velocity
        p->v[bi+0] += p->a[bi+0] * p->dt;
        p->v[bi+1] += p->a[bi+1] * p->dt;
        p->v[bi+2] += p->a[bi+2] * p->dt;

        // position
        p->x[bi+0] += p->v[bi+0] * p->dt;
        p->x[bi+1] += p->v[bi+1] * p->dt;
        p->x[bi+2] += p->v[bi+2] * p->dt;

        // angular acceleration
        p->anga[bi+0] = p->mom[bi+0]*p->invmoi[i];
        p->anga[bi+1] = p->mom[bi+1]*p->invmoi[i];
        p->anga[bi+2] = p->mom[bi+2]*p->invmoi[i];

        // angular velocity
        p->angv[bi+0] += p->anga[bi+0] * p->dt;
        p->angv[bi+1] += p->anga[bi+1] * p->dt;
        p->angv[bi+2] += p->anga[bi+2] * p->dt;
    }
    int flag = shouldRefresh(p,box);
    if (flag ==1){
        update_neighborlist(p,tmpP,box);
    }
}

void cpu_dem_sort_triangles(ParticleSystem* ps, ParticleSystem *tmpPs, BoundingBox* box,TriangleMesh *mesh, int step){
    /* initialize */
    for (int i=0; i<ps->N*DIM; i++){
        ps->f[i]=0.;
    }

    for (int i=0; i<ps->N*DIM; i++){
        ps->mom[i]=0.;
    }

    for (int i=0; i<ps->N*ps->MAX_NEI; i++){
        ps->isContact[i]=0;
        ps->isContactWall[i]=0;
    }


    int reorderFreq=100;

    if(step%(reorderFreq)==0){
        particle_collision_cell_linked_withSort_fastUpdate(ps,tmpPs,box);
    }else{
        particle_collision_cell_linked_fastUpdate(ps,box);
    }

    wall_collision_triangles(ps,box,mesh);

    /* update */
    for (int i = 0; i < ps->N; i++)
    {
        if (ps->isActive[i]!=1){
            continue;
        }
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

        // angular acceleration
        ps->anga[bi+0] = ps->mom[bi+0]*ps->invmoi[i];
        ps->anga[bi+1] = ps->mom[bi+1]*ps->invmoi[i];
        ps->anga[bi+2] = ps->mom[bi+2]*ps->invmoi[i];

        // angular velocity
        ps->angv[bi+0] += ps->anga[bi+0] * ps->dt;
        ps->angv[bi+1] += ps->anga[bi+1] * ps->dt;
        ps->angv[bi+2] += ps->anga[bi+2] * ps->dt;

    }
}

void cpu_dem_naive_triangle(ParticleSystem* ps, BoundingBox* box, TriangleMesh *mesh){
    /* initialize */
    for (int i=0; i<ps->N*DIM; i++){
        ps->f[i]=0.;
    }

    for (int i=0; i<ps->N*DIM; i++){
        ps->mom[i]=0.;
    }

    for (int i=0; i<ps->N*ps->MAX_NEI; i++){
        ps->isContact[i]=0;
        ps->isContactWall[i]=0;
    }



    update_pList_fast(ps,box);
    particle_collision_naive(ps);

    wall_collision_triangles(ps,box,mesh);

    /* update */
    for (int i = 0; i < ps->N; i++)
    {
        if (ps->isActive[i]!=1){
            continue;
        }
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

        // angular acceleration
        ps->anga[bi+0] = ps->mom[bi+0]*ps->invmoi[i];
        ps->anga[bi+1] = ps->mom[bi+1]*ps->invmoi[i];
        ps->anga[bi+2] = ps->mom[bi+2]*ps->invmoi[i];

        // angular velocity
        ps->angv[bi+0] += ps->anga[bi+0] * ps->dt;
        ps->angv[bi+1] += ps->anga[bi+1] * ps->dt;
        ps->angv[bi+2] += ps->anga[bi+2] * ps->dt;

    }
}

void cpu_dem_nosort_triangle(ParticleSystem* ps, ParticleSystem *tmpPs, BoundingBox* box,TriangleMesh *mesh){
    /* initialize */
    for (int i=0; i<ps->N*DIM; i++){
        ps->f[i]=0.;
    }

    for (int i=0; i<ps->N*DIM; i++){
        ps->mom[i]=0.;
    }

    for (int i=0; i<ps->N*ps->MAX_NEI; i++){
        ps->isContact[i]=0;
        ps->isContactWall[i]=0;
    }



    particle_collision_cell_linked_fastUpdate(ps,box);

    wall_collision_triangles(ps,box,mesh);

    /* update */
    for (int i = 0; i < ps->N; i++)
    {
        if (ps->isActive[i]!=1){
            continue;
        }
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

        // angular acceleration
        ps->anga[bi+0] = ps->mom[bi+0]*ps->invmoi[i];
        ps->anga[bi+1] = ps->mom[bi+1]*ps->invmoi[i];
        ps->anga[bi+2] = ps->mom[bi+2]*ps->invmoi[i];

        // angular velocity
        ps->angv[bi+0] += ps->anga[bi+0] * ps->dt;
        ps->angv[bi+1] += ps->anga[bi+1] * ps->dt;
        ps->angv[bi+2] += ps->anga[bi+2] * ps->dt;

    }
}

void cpu_dem_sort(ParticleSystem* ps, ParticleSystem *tmpPs, BoundingBox* box, int step){
    /* initialize */
    for (int i=0; i<ps->N*DIM; i++){
        ps->f[i]=0.;
    }

    for (int i=0; i<ps->N*DIM; i++){
        ps->mom[i]=0.;
    }

    for (int i=0; i<ps->N*ps->MAX_NEI; i++){
        ps->isContact[i]=0;
        ps->isContactWall[i]=0;
    }


    int reorderFreq=100;

    if(step%(reorderFreq)==0){
        particle_collision_cell_linked_withSort_fastUpdate(ps,tmpPs,box);
    }else{
        particle_collision_cell_linked_fastUpdate(ps,box);
    }

    wall_collision_naive(ps);

    /* update */
    for (int i = 0; i < ps->N; i++)
    {
        if (ps->isActive[i]!=1){
            continue;
        }
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

        // angular acceleration
        ps->anga[bi+0] = ps->mom[bi+0]*ps->invmoi[i];
        ps->anga[bi+1] = ps->mom[bi+1]*ps->invmoi[i];
        ps->anga[bi+2] = ps->mom[bi+2]*ps->invmoi[i];

        // angular velocity
        ps->angv[bi+0] += ps->anga[bi+0] * ps->dt;
        ps->angv[bi+1] += ps->anga[bi+1] * ps->dt;
        ps->angv[bi+2] += ps->anga[bi+2] * ps->dt;

    }
}

