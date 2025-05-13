#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../common/my_timer.h"

#define MAX_PD	(3.0e6)
#define PRECISION	0.001
#define SPEC_HEAT_SI 1.75e6
#define K_SI 100
#define FACTOR_CHIP	0.5
#define AMB_TEMP 80.0

__global__ void computeTemp(float *p, float* tIn, float *tOut, float sdc,
        int nx, int ny, int nz,
        float ce, float cw, 
        float cn, float cs,
        float ct, float cb, 
        float cc) 
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;  
    int j = blockDim.y * blockIdx.y + threadIdx.y;

    int c = i + j * nx;
    int xy = nx * ny;

    int W = (i == 0)        ? c : c - 1;
    int E = (i == nx-1)     ? c : c + 1;
    int N = (j == 0)        ? c : c - nx;
    int S = (j == ny-1)     ? c : c + nx;

    float bottom, center, top;
    bottom = center = tIn[c];
    top = tIn[c+xy];
    tOut[c] = cc * center + cw * tIn[W] + ce * tIn[E] + cs * tIn[S]
        + cn * tIn[N] + cb * bottom + ct * top + sdc * p[c] + ct * AMB_TEMP;
    c += xy;
    W += xy;
    E += xy;
    N += xy;
    S += xy;

    for (int k = 1; k < nz-1; ++k) {
        bottom = center;
        center = top;
        top = tIn[c+xy];
        tOut[c] = cc * center + cw * tIn[W] + ce * tIn[E] + cs * tIn[S]
            + cn * tIn[N] + cb * bottom + ct * top + sdc * p[c] + ct * AMB_TEMP;
        c += xy;
        W += xy;
        E += xy;
        N += xy;
        S += xy;
    }
    bottom = center;
    center = top;
    tOut[c] = cc * center + cw * tIn[W] + ce * tIn[E] + cs * tIn[S]
        + cn * tIn[N] + cb * bottom + ct * top + sdc * p[c] + ct * AMB_TEMP;
    return;
}

void computeTempGPU(float *p, float *tIn, float *tOut,
        int nx, int ny, int nz,
        float Cap, 
        float Rx, float Ry, float Rz, 
        float dt, int numiter) 
{
    float ce, cw, cn, cs, ct, cb, cc;
    float stepDivCap = dt / Cap;
    ce = cw =stepDivCap/ Rx;
    cn = cs =stepDivCap/ Ry;
    ct = cb =stepDivCap/ Rz;

    cc = 1.0 - (2.0*ce + 2.0*cn + 3.0*ct);

    size_t s = sizeof(float) * nx * ny * nz;  
    float  *tIn_d, *tOut_d, *p_d;

	cudaStream_t stream;
	cudaStreamCreate(&stream);

	timespec timer = tic();

    cudaMalloc((void**)&p_d,s);
    cudaMalloc((void**)&tIn_d,s);
    cudaMalloc((void**)&tOut_d,s);
    cudaMemcpyAsync(tIn_d, tIn, s, cudaMemcpyHostToDevice, stream);
    cudaMemcpyAsync(p_d, p, s, cudaMemcpyHostToDevice, stream);

    cudaFuncSetCacheConfig(computeTemp, cudaFuncCachePreferL1);

    dim3 block_dim(64, 4, 1);
    dim3 grid_dim(nx / 64, ny / 4, 1);

    for (int i = 0; i < numiter; ++i) {
        computeTemp<<<grid_dim, block_dim>>>
            (p_d, tIn_d, tOut_d, stepDivCap, nx, ny, nz, ce, cw, cn, cs, ct, cb, cc);

        float *t = tIn_d;
        tIn_d = tOut_d;
        tOut_d = t;
    }

	cudaMemcpyAsync(tOut, tOut_d, s, cudaMemcpyDeviceToHost, stream);
	cudaStreamSynchronize(stream);

	toc(&timer, "GPU execution time (including data transfer) in seconds");

    cudaFree(p_d);
    cudaFree(tIn_d);
    cudaFree(tOut_d);
    return;
}

void initializeArrays(float* pIn, float* tIn, const int nx, const int ny, const int nz)
{
	int z,x,y;	

	for (z = 0; z < nz; z++) {
		for (x = 0; x < nx; x++) {
			for (y = 0; y < ny; y++) {
				pIn[z * nx * ny + x * ny + y] = (float)((z * x * y + 1) % nz) / nz;
			}
		}
	}

	for (z = 0; z < nz; z++) {
		for (x = 0; x < nx; x++) {
			for (y = 0; y < ny; y++) {
				tIn[z * nx * ny + x * ny + y] = (float)((z * x * y + 1) % nx) / nx;
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

void computeTempCPU(float *pIn, float* tIn, float *tOut, 
        int nx, int ny, int nz, float Cap, 
        float Rx, float Ry, float Rz, 
        float dt, int numiter) 
{   float ce, cw, cn, cs, ct, cb, cc;
    float stepDivCap = dt / Cap;
    ce = cw =stepDivCap/ Rx;
    cn = cs =stepDivCap/ Ry;
    ct = cb =stepDivCap/ Rz;

    cc = 1.0 - (2.0*ce + 2.0*cn + 3.0*ct);

    int c,w,e,n,s,b,t;
    int x,y,z;
    int i = 0;
    do{
        for(z = 0; z < nz; z++)
            for(y = 0; y < ny; y++)
                for(x = 0; x < nx; x++)
                {
                    c = x + y * nx + z * nx * ny;

                    w = (x == 0) ? c      : c - 1;
                    e = (x == nx - 1) ? c : c + 1;
                    n = (y == 0) ? c      : c - nx;
                    s = (y == ny - 1) ? c : c + nx;
                    b = (z == 0) ? c      : c - nx * ny;
                    t = (z == nz - 1) ? c : c + nx * ny;


                    tOut[c] = tIn[c]*cc + tIn[n]*cn + tIn[s]*cs + tIn[e]*ce + tIn[w]*cw + tIn[t]*ct + tIn[b]*cb + (dt/Cap) * pIn[c] + ct*AMB_TEMP;
                }
        float *temp = tIn;
        tIn = tOut;
        tOut = temp; 
        i++;
    }
    while(i < numiter);
}

int main(int argc, char* argv[]) {
	if (argc != 4) {
		printf("Usage: ./computeTempFPGA <NX/NY> <NZ> <NUM_ITERATIONS>\n");
	}

	const int NX = strtol(argv[1], NULL, 10);
	const int NY = NX;
	const int NZ = strtol(argv[2], NULL, 10);
	const int NUM_ITERATIONS = strtol(argv[3], NULL, 10);

	const float t_chip = 0.0005;
	const float chip_height = 0.016;
	const float chip_width = 0.016; 
    const float dx = chip_height / NX;
    const float dy = chip_width / NY;
    const float dz = t_chip / NZ;
    const float Cap = FACTOR_CHIP * SPEC_HEAT_SI * t_chip * dx * dy;
    const float Rx = dy / (2.0 * K_SI * t_chip * dx);
    const float Ry = dx / (2.0 * K_SI * t_chip * dy);
    const float Rz = dz / (K_SI * dx * dy);
    const float max_slope = MAX_PD / (FACTOR_CHIP * t_chip * SPEC_HEAT_SI);
    const float dt = PRECISION / max_slope;

    const int size = NX * NY * NZ;

    float* powerIn;
	cudaHostAlloc(&powerIn, size * sizeof(float), cudaHostAllocDefault);
    float* tempIn;
	cudaHostAlloc(&tempIn, size * sizeof(float), cudaHostAllocDefault);
    float* tempOut;
	cudaHostAlloc(&tempOut, size * sizeof(float), cudaHostAllocDefault);
    float* tempCopy = (float*)malloc(size * sizeof(float));
    float* answer = (float*)calloc(size, sizeof(float));

	initializeArrays(powerIn, tempIn, NX, NY, NZ);
    memcpy(tempCopy, tempIn, size * sizeof(float));

	computeTempGPU(powerIn, tempIn, tempOut, NX, NY, NZ, Cap, Rx, Ry, Rz, dt, NUM_ITERATIONS);
	cudaFreeHost(tempIn);
	computeTempCPU(powerIn, tempCopy, answer, NX, NY, NZ, Cap, Rx, Ry, Rz, dt, NUM_ITERATIONS);
	printf("Accuracy: %e\n", accuracy(tempOut, answer, size));

	cudaFreeHost(powerIn);
	cudaFreeHost(tempOut);
	free(tempCopy);
	free(answer);

	return 0;
}
