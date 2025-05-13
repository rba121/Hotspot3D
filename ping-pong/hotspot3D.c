#include <stdio.h>
#include <time.h>
#include <assert.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h> 
#include <sys/time.h>
#include <string.h>

#include "hotspot3D.h"

#define AMB_TEMP 80.0
#define TEMP_COLS (TILE_COLS + 2)

void initTempBuffer(float tempBuffer[3][3][TEMP_COLS], float* tIn, const int tileX)
{
	int z,y;
	for (z = 0; z < 3; z++) {
		for (y = 0; y < 2; y++) {
			memcpy(tempBuffer[z][y + 1], tIn + z * NY * NX + y * NX + tileX - 1, TEMP_COLS * sizeof(float));
		}
	}
}

void shiftRows(float tempBuffer[3][3][TEMP_COLS], float* tIn, const int tileX, const int y, const Direction direction)
{
	int z;
	for (z = 0; z < 3; z++) {
		memcpy(tempBuffer[z][0], tempBuffer[z][1], TEMP_COLS * sizeof(float));
		memcpy(tempBuffer[z][1], tempBuffer[z][2], TEMP_COLS * sizeof(float));

		const int layer = (direction == FRONT_TO_BACK) ? z : NZ - 1 - (2 - z);
		memcpy(tempBuffer[z][2], tIn + layer * NY * NX + (y + 1) * NX + tileX - 1, TEMP_COLS * sizeof(float));
	}
}

void loadPower(float* powerBuffer, float* pIn, const int y, const int tileX, const int z, const Direction direction)
{
	const int layer = (direction == FRONT_TO_BACK) ? z : NZ - 1 - z;
	memcpy(powerBuffer, pIn + layer * NY * NX + y * NX + tileX, TILE_COLS * sizeof(float)); 
}

void shiftLayers(float tempBuffer[3][3][TEMP_COLS], float* tIn, const int tileX, const int z, const int startRow, const Direction direction)
{
	int y;

	if (direction == FRONT_TO_BACK) {
		memcpy(tempBuffer[0], tempBuffer[1], 3 * TEMP_COLS * sizeof(float));
		memcpy(tempBuffer[1], tempBuffer[2], 3 * TEMP_COLS * sizeof(float));
		for (y = 0; y < 3; y++) {
			memcpy(tempBuffer[2][y], tIn + z * NY * NX + (y + startRow) * NX + tileX - 1, TEMP_COLS * sizeof(float));
		}
	} else {
		memcpy(tempBuffer[2], tempBuffer[1], 3 * TEMP_COLS * sizeof(float));
		memcpy(tempBuffer[1], tempBuffer[0], 3 * TEMP_COLS * sizeof(float));
		for (y = 0; y < 3; y++) {
			memcpy(tempBuffer[0][y], tIn + (NZ - 1 - z) * NY * NX + (y + startRow) * NX + tileX - 1, TEMP_COLS * sizeof(float));
		}
	}
}



void compute(float tempBuffer[3][3][TEMP_COLS], float* powerBuffer, float* resultBuffer, const float ce, const float cw, const float cn, const float cs, const float ct, const float cb, const float cc, const float dt, const float Cap)
{
	for (int x = 1; x < TEMP_COLS - 1; x++) {	
		resultBuffer[x - 1] = tempBuffer[1][1][x] * cc
							+ tempBuffer[1][1][x - 1] * cw
							+ tempBuffer[1][1][x + 1] * ce
							+ tempBuffer[1][0][x] * cn
							+ tempBuffer[1][2][x] * cs
							+ tempBuffer[0][1][x] * cb
							+ tempBuffer[2][1][x] * ct
							+ (dt/Cap) * powerBuffer[x - 1]
							+ ct * AMB_TEMP;
	}
}

void store(float* tOut, float* resultBuffer, const int tileX, const int y, const int z, const Direction direction)
{
	const int layer = (direction == FRONT_TO_BACK) ? z : NZ - 1 - z;
	memcpy(tOut + layer * NY * NX + y * NX + tileX, resultBuffer, TILE_COLS * sizeof(float));
}

void computeTempBorders(float *pIn, float* tIn, float *tOut, int nx, int ny, int nz, float Cap, 
        float Rx, float Ry, float Rz, 
        float dt) 
{   float ce, cw, cn, cs, ct, cb, cc;
    float stepDivCap = dt / Cap;
    ce = cw =stepDivCap/ Rx;
    cn = cs =stepDivCap/ Ry;
    ct = cb =stepDivCap/ Rz;

    cc = 1.0 - (2.0*ce + 2.0*cn + 3.0*ct);

    int c,w,e,n,s,b,t;
    int x,y,z;
        for(z = 0; z < nz; z++)
            for(y = 0; y < ny; y++)
                for(x = 0; x < nx; x++)
                {
					if (!(z == 0 || z == nz - 1 || x == 0 || x == nx - 1 || y == 0 || y == ny - 1)) {
						continue;
					}

                    c = x + y * nx + z * nx * ny;

                    w = (x == 0) ? c      : c - 1;
                    e = (x == nx - 1) ? c : c + 1;
                    n = (y == 0) ? c      : c - nx;
                    s = (y == ny - 1) ? c : c + nx;
                    b = (z == 0) ? c      : c - nx * ny;
                    t = (z == nz - 1) ? c : c + nx * ny;


                    tOut[c] = tIn[c]*cc + tIn[n]*cn + tIn[s]*cs + tIn[e]*ce + tIn[w]*cw + tIn[t]*ct + tIn[b]*cb + (dt/Cap) * pIn[c] + ct*AMB_TEMP;
                }
}

void loadLayers(float destBuffer[3][3][TEMP_COLS], float srcBuffer[3][3][TEMP_COLS], float* tIn, const int tileX, const int startRow, const int z, Direction direction){
	int y;

	if (direction == FRONT_TO_BACK) {
		memcpy(destBuffer[0], srcBuffer[1], 3 * TEMP_COLS * sizeof(float));
		memcpy(destBuffer[1], srcBuffer[2], 3 * TEMP_COLS * sizeof(float));
		for (y = 0; y < 3; y++) {
			memcpy(destBuffer[2][y], tIn + z * NY * NX + (y + startRow) * NX + tileX - 1, TEMP_COLS * sizeof(float));
		}
	} else {
		memcpy(destBuffer[2], srcBuffer[1], 3 * TEMP_COLS * sizeof(float));
		memcpy(destBuffer[1], srcBuffer[0], 3 * TEMP_COLS * sizeof(float));
		for (y = 0; y < 3; y++) {
			memcpy(destBuffer[0][y], tIn + (NZ - 1 - z) * NY * NX + (y + startRow) * NX + tileX - 1, TEMP_COLS * sizeof(float));
		}
	}
}


void computeTempFPGA(float *pIn, float* tIn, float *tOut, 
		float Cap, float Rx, float Ry,
		float Rz, float dt, int numiter) 
{
	float ce, cw, cn, cs, ct, cb, cc;
    float stepDivCap = dt / Cap;
    ce = cw =stepDivCap/ Rx;
    cn = cs =stepDivCap/ Ry;
    ct = cb =stepDivCap/ Rz;

    cc = 1.0 - (2.0*ce + 2.0*cn + 3.0*ct);

    int x,y,z;
    int i = 0;
	float tempBuffer[2][3][3][TEMP_COLS];
	#pragma HLS array_partition complete dim=1 variable=tempBuffer
	#pragma HLS array_partition complete dim=2 variable=tempBuffer
	float powerBuffer[2][TILE_COLS];
	float resultBuffer[2][TILE_COLS];

	Direction direction;
	
    do {
		// This function is not optimized (testing only)
		computeTempBorders(pIn, tIn, tOut, NX, NY, NZ, Cap, Rx, Ry, Rz, dt);
	
		for (x = 1; x < NX - 1; x += TILE_COLS) {
			initTempBuffer(tempBuffer[0], tIn, x);
			for (y = 1; y < NY - 1; y++) {
				direction = (y % 2) == 1 ? FRONT_TO_BACK : BACK_TO_FRONT;
				//Load initial power so that resultBuffer[0] can be computed
				//Currently I am setting NZ = 9 because I want to make sure that the copy in the last layer is the same copy I am applying shift rows to
				shiftRows(tempBuffer[0], tIn, x, y, direction);
				loadPower(powerBuffer[0], pIn, y, x, z -1, direction);
				for (z = 3; z < NZ; z++) {
					if(z == 3){
						//At this point there is nothing to store because there was no computation in the initial setup
						compute(tempBuffer[0], powerBuffer[0], resultBuffer[0], ce, cw, cn, cs, ct, cb, cc, dt, Cap);	
						loadPower(powerBuffer[1], pIn, y, x, z-1, direction);
						loadLayers(tempBuffer[1],tempBuffer[0], tIn, x, y-1, z, direction);			
					}else{
						switch((z) %2){
							//Use modulus for the indicies for the copies because the copies can fall out of sync with the z iterator
							case 0:
								store(tOut, resultBuffer[0], x, y, z - 2, direction);  // Store to z -2 because you are storing computed values of previous z. Store operations are always 1 iteration behind.
								loadPower(powerBuffer[0], pIn, y, x, z - 1, direction);
								loadLayers(tempBuffer[0], tempBuffer[1], tIn, x, y-1, z, direction);
								compute(tempBuffer[1], powerBuffer[1], resultBuffer[1], ce, cw, cn, cs, ct, cb, cc, dt, Cap);	
								break;

							case 1:

								store(tOut, resultBuffer[1], x, y, z-2, direction);
	
								loadLayers(tempBuffer[1], tempBuffer[0], tIn, x, y-1, z, direction); //Load current layers
								loadPower(powerBuffer[1], pIn, y, x, z - 1, direction);  //load current power
								compute(tempBuffer[0], powerBuffer[0], resultBuffer[0], ce, cw, cn, cs, ct, cb, cc, dt, Cap); //Compute from previously loaded buffers	
								break;
						}

					}
				}
				compute(tempBuffer[0], powerBuffer[0], resultBuffer[0], ce, cw, cn, cs, ct, cb, cc, dt, Cap); //Compute the last loaded  copy
				store(tOut, resultBuffer[0], x, y, z-2, direction); //Store the last loaded copy
				store(tOut, resultBuffer[1], x, y, z-3, direction); //Store the computed copy 2 iterations ago
			}
		}

        float *temp = tIn;
        tIn = tOut;
        tOut = temp; 
        i++;
    }
    while(i < numiter);

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
