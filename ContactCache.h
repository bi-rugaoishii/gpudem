
#ifndef _CONTACTCACHE_H_
#define _CONTACTCACHE_H_

#include "Vec3.h"

typedef struct ContactCache{
    Vec3 fn;
    Vec3 vn_rel;
    Vec3 vt;

    #if USE_GPU
    __host__ __device__
    ContactCache(){}

    #endif

}ContactCache;

#endif
