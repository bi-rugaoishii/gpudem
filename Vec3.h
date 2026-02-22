#ifndef _VEC3_H_
#define _VEC3_H_

typedef struct Vec3{
    double x;
    double y;
    double z;

#if USE_GPU
    __host__ __device__
        Vec3(){}
#endif

#if USE_GPU
    __host__ __device__
        Vec3(double X,double Y,double Z)
        : x(X),y(Y),z(Z){}
#endif

}Vec3;

#if USE_GPU
__host__ __device__
#endif
inline double vdot(Vec3 a,Vec3 b){
    return a.x*b.x+a.y*b.y+a.z*b.z;
}

#if USE_GPU
__host__ __device__
#endif
inline Vec3 vscalar(double a,Vec3 b){
    Vec3 c;
    c.x = a*b.x;
    c.y = a*b.y;
    c.z = a*b.z;
    return c;
}

#if USE_GPU
__host__ __device__
#endif
inline Vec3 vadd(Vec3 a,Vec3 b){
    Vec3 c;
    c.x = a.x+b.x;
    c.y = a.y+b.y;
    c.z = a.z+b.z;
    return c;
}

#if USE_GPU
__host__ __device__
#endif
inline Vec3 vsub(Vec3 a,Vec3 b){
    Vec3 c;
    c.x = a.x-b.x;
    c.y = a.y-b.y;
    c.z = a.z-b.z;
    return c;
}

#if USE_GPU
__host__ __device__
#endif
inline Vec3 vcross(Vec3 a,Vec3 b){
    Vec3 c;
    c.x = a.y*b.z-a.z*b.y;
    c.y = a.z*b.x-a.x*b.z;
    c.z = a.x*b.y-a.y*b.x;

    return c;
}


#endif
