#include <math.h>
#include <stdio.h>
#include <stdlib.h> 
#include <string.h>

#include "hotspot3D.h"
#include "my_timer.h"

#define MAX_PD	(3.0e6)
// required precision in degrees
#define PRECISION	0.001
#define SPEC_HEAT_SI 1.75e6
#define K_SI 100
// capacitance fitting factor
#define FACTOR_CHIP	0.5

/* chip parameters	*/
float t_chip = 0.0005;
float chip_height = 0.016;
float chip_width = 0.016; 

void initializeArrays(float* pIn, float* tIn)
{
	int z,x,y;	

	for (z = 0; z < NZ; z++) {
		for (x = 0; x < NX; x++) {
			for (y = 0; y < NY; y++) {
				pIn[z * NX * NY + x * NY + y] = (float)((z * x * y + 1) % NZ) / NZ;
			}
		}
	}

	for (z = 0; z < NZ; z++) {
		for (x = 0; x < NX; x++) {
			for (y = 0; y < NY; y++) {
				tIn[z * NX * NY + x * NY + y] = (float)((z * x * y + 1) % NX) / NX;
			}
		}
	}
}

float accuracy(float *arr1, float *arr2, int len)
{
    float err = 0.0; 
    int i;
    for(i = 0; i < len; i++)
    {
        err += abs(arr1[i] - arr2[i]);
    }

    return (float)err / len;
}

int main(void)
{
    float dx = chip_height / NX;
    float dy = chip_width / NY;
    float dz = t_chip / NZ;
    float Cap = FACTOR_CHIP * SPEC_HEAT_SI * t_chip * dx * dy;
    float Rx = dy / (2.0 * K_SI * t_chip * dx);
    float Ry = dx / (2.0 * K_SI * t_chip * dy);
    float Rz = dz / (K_SI * dx * dy);
    float max_slope = MAX_PD / (FACTOR_CHIP * t_chip * SPEC_HEAT_SI);
    float dt = PRECISION / max_slope;

    int size = NX * NY * NZ;

    float* powerIn = (float*)calloc(size, sizeof(float));
    float* tempCopy = (float*)malloc(size * sizeof(float));
    float* tempIn = (float*)calloc(size,sizeof(float));
    float* tempOut = (float*)calloc(size, sizeof(float));
    float* answer = (float*)calloc(size, sizeof(float));

	initializeArrays(powerIn, tempIn);
    memcpy(tempCopy, tempIn, size * sizeof(float));

	timespec timer = tic();
    computeTempFPGA(powerIn, tempIn, tempOut, Cap, Rx, Ry, Rz, dt, NUM_ITERATIONS);
	toc(&timer, "Execution Time (s)");
   // computeTempFPGA(powerIn, tempCopy, answer, Cap, Rx,Ry,Rz,dt,NUM_ITERATIONS);
   computeTempCPU(powerIn, tempCopy, answer, NX, NY, NZ, Cap, Rx, Ry, Rz, dt, NUM_ITERATIONS);

    float acc = accuracy(tempOut, answer, size);
    printf("Accuracy: %e\n", acc);

	free(answer);
	free(tempOut);
	free(tempIn);
	free(tempCopy);
	free(powerIn);

	return 0;
}
