#ifndef HOTSPOT3D_H
#define HOTSPOT3D_H

#define NX 514
#define NY 514
#define NZ 8

#define NUM_ITERATIONS 100

#define TILE_COLS 128

typedef enum {
	FRONT_TO_BACK,
	BACK_TO_FRONT,	
} Direction;

void computeTempFPGA(float *pIn, float* tIn, float *tOut, float Cap, float Rx, float Ry, float Rz, float dt, int numiter);
void computeTempCPU(float *pIn, float* tIn, float *tOut, int nx, int ny, int nz, float Cap, float Rx, float Ry, float Rz, float dt, int numiter);

#endif
