#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#include <math.h>
#include <chrono>
#include "ParticleSystem.h"
#ifndef _OUTPUT_H_
#define _OUTPUT_H_
#define DIM 3




/*
============================================================
出力
============================================================
*/
void writeParticles(ParticleSystem* ps, int step);

/*
============================================================
VTK Binary出力（ParaView用・高速）
============================================================
*/

/* ============================================================
   エンディアン変換（little → big）
============================================================ */
void swapBytes(void* data, size_t size);
void writeParticlesVTKBinary(ParticleSystem* ps, int step);

/*
============================================================
VTK出力（ParaView用）
============================================================
*/
void writeParticlesVTK(ParticleSystem* ps, int step);
void writeParticlesDimensionalizeVTK(ParticleSystem* ps, int step);

#endif
