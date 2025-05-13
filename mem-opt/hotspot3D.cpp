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

void computeVerticalStrip(float* pIn, float* tIn, float* tOut, const float ce, const float cw, const float cn, const float cs, const float ct, const float cb, const float cc, const float dt, const float Cap, const int x, const int z)
{
	int y;
	float south, center, north, west, east;
	int c = z * NY * NX + x; north = center = tIn[c];
	south = tIn[c + NX];
	west = (x == 0) ? center : tIn[c - 1];
	east = (x == NX - 1) ? center : tIn[c + 1];
	
	tOut[c] = center * cc + west * cw + east * ce + north * cn + south * cs + tIn[c - NY * NX] * cb + tIn[c + NY * NX] * ct + (dt/Cap) * pIn[c] + ct * AMB_TEMP;

	c += NX;

	for (y = 1; y < NY - 1; y++) {
		north = center;
		center = south;
		south = tIn[c + NX];

		west = (x == 0) ? center : tIn[c - 1];
		east = (x == NX - 1) ? center : tIn[c + 1];

		tOut[c] = center * cc + west * cw + east * ce + north * cn + south * cs + tIn[c - NY * NX] * cb + tIn[c + NY * NX] * ct + (dt/Cap) * pIn[c] + ct * AMB_TEMP;
		
		c += NX;
	}
	
	north = center;
	center = south;
	
	west = (x == 0) ? center : tIn[c - 1];
	east = (x == NX - 1) ? center : tIn[c + 1];

	tOut[c] = center * cc + west * cw + east * ce + north * cn + south * cs + tIn[c - NY * NX] * cb + tIn[c + NY * NX] * ct + (dt/Cap) * pIn[c] + ct * AMB_TEMP;
}

void computeTopFace(float* pIn, float* tIn, float* tOut, float tempBuffer[3][3][TEMP_COLS], float* powerBuffer, float* resultBuffer, float ce, float cw, float cn, float cs, float ct, float cb, float cc, float dt, float Cap)
{
	// TOP EDGE
	int x;
	for (x = 0; x < NX; x++) {
		const int c = (NZ - 1) * NY * NX + x;
		float center = tIn[c];
		float west = (x == 0)		? center : tIn[c - 1];
		float east = (x == NX - 1) 	? center : tIn[c + 1];
		tOut[c] = center * cc + west * cw + east * ce + center * cn + tIn[c + NX] * cs + tIn[c - NY * NX] * cb + center * ct
			+ (dt/Cap) * pIn[c] + ct * AMB_TEMP; 
	}

	// LEFT VERTICAL STRIP
	int y;
	float north, center, south; 
	int c = (NZ - 1) * NY * NX + NX;
	north = tIn[c - NX];
	center = tIn[c];
	south = tIn[c + NX];

	for (y = 1; y < NY - 1; y++) {
		tOut[c] = center * cc + center * cw + tIn[c + 1] * ce + north * cn + south * cs + tIn[c - NY * NX] * cb + center * ct 
			+ (dt/Cap) * pIn[c] + ct * AMB_TEMP;

		north = center;
		center = south;
		c += NX;
		south = tIn[c + NX];
	}

	// RIGHT VERTICAL STRIP
	c = (NZ - 1) * NY * NX + 2 * NX - 1;
	north = tIn[c - NX];
	center = tIn[c];
	south = tIn[c + NX];
	
	for (y = 1; y < NY - 1; y++) {
		tOut[c] = center * cc + tIn[c - 1] * cw + center * ce + north * cn + south * cs + tIn[c - NY * NX] * cb + center * ct
			+ (dt/Cap) * pIn[c] + ct * AMB_TEMP;

		north = center;
		center = south;
		c += NX;
		south = tIn[c + NX];
	}

	// INNER SQUARE
	const int z = (NZ - 1) * NY * NX;
	for (x = 1; x < NX - 1; x += TILE_COLS) {
		memcpy(tempBuffer[0][1], tIn + z + (x - 1), TEMP_COLS * sizeof(float));
		memcpy(tempBuffer[1][1], tIn + z - NY * NX + (x - 1), TEMP_COLS * sizeof(float));
		memcpy(tempBuffer[0][2], tIn + z + NX + (x - 1), TEMP_COLS * sizeof(float));
		memcpy(tempBuffer[1][2], tIn + z - NY * NX + NX + (x - 1), TEMP_COLS * sizeof(float));
		for (y = 1; y < NY - 1; y++) {
			// shift
			memcpy(tempBuffer[0][0], tempBuffer[0][1], TEMP_COLS * sizeof(float));
			memcpy(tempBuffer[1][0], tempBuffer[1][1], TEMP_COLS * sizeof(float));
			memcpy(tempBuffer[0][1], tempBuffer[0][2], TEMP_COLS * sizeof(float));
			memcpy(tempBuffer[1][1], tempBuffer[1][2], TEMP_COLS * sizeof(float));
			memcpy(tempBuffer[0][2], tIn + z + (y + 1) * NX + (x - 1), TEMP_COLS * sizeof(float));
			memcpy(tempBuffer[1][2], tIn + z - NY * NX + (y + 1) * NX + (x - 1), TEMP_COLS * sizeof(float));
			// load power
			memcpy(powerBuffer, pIn + z + y * NX + x, TILE_COLS * sizeof(float));
			// compute
			int i;
			for (i = 1; i < TEMP_COLS - 1; i++) {
				resultBuffer[i - 1] = tempBuffer[0][1][i] * cc + tempBuffer[0][1][i - 1] * cw + tempBuffer[0][1][i + 1] * ce
					+ tempBuffer[0][0][i] * cn + tempBuffer[0][2][i] * cs + tempBuffer[1][1][i] * cb + tempBuffer[0][1][i] * ct
					+ (dt/Cap) * powerBuffer[i - 1] + ct * AMB_TEMP;
			}
			// store
			memcpy(tOut + z + y * NX + x, resultBuffer, TILE_COLS * sizeof(float));
		}
	}

	// BOTTOM EDGE
	for (x = 0; x < NX; x++) {
		const int c = z + (NY - 1) * NX + x;
		float center = tIn[c];
		float west = (x == 0)		? center : tIn[c - 1];
		float east = (x == NX - 1) 	? center : tIn[c + 1];
		tOut[c] = center * cc + west * cw + east * ce + tIn[c - NX] * cn + center * cs + tIn[c - NY * NX] * cb + center * ct
			+ (dt/Cap) * pIn[c] + ct * AMB_TEMP; 
	}
}

void computeBottomFace(float* pIn, float* tIn, float* tOut, float tempBuffer[3][3][TEMP_COLS], float* powerBuffer, float* resultBuffer, float ce, float cw, float cn, float cs, float ct, float cb, float cc, float dt, float Cap)
{
	// TOP EDGE
	int x;
	for (x = 0; x < NX; x++) {	
		float center = tIn[x];
		float west = (x == 0)		? center : tIn[x - 1];
		float east = (x == NX - 1) 	? center : tIn[x + 1];
		tOut[x] = center * cc + west * cw + east * ce + center * cn + tIn[x + NX] * cs + center * cb + tIn[x + NY * NX] * ct
			+ (dt/Cap) * pIn[x] + ct * AMB_TEMP; 
	}

	// LEFT VERTICAL STRIP
	int y;
	float north, center, south; 
	int c = NX;
	north = tIn[c - NX];
	center = tIn[c];
	south = tIn[c + NX];

	for (y = 1; y < NY - 1; y++) {
		tOut[c] = center * cc + center * cw + tIn[c + 1] * ce + north * cn + south * cs + center * cb + tIn[c + NX * NY] * ct 
			+ (dt/Cap) * pIn[c] + ct * AMB_TEMP;

		north = center;
		center = south;
		south = tIn[(y + 1) * NX];
		c += NX;
	}

	// RIGHT VERTICAL STRIP
	c = 2 * NX - 1;
	north = tIn[c - NX];
	center = tIn[c];
	south = tIn[c + NX];
	
	for (y = 1; y < NY - 1; y++) {
		tOut[c] = center * cc + tIn[c - 1] * cw + center * ce + north * cn + south * cs + center * cb + tIn[c + NX * NY] * ct
			+ (dt/Cap) * pIn[c] + ct * AMB_TEMP;

		north = center;
		center = south;
		south = tIn[(y + 1) * NX + NX - 1];
	}

	// INNER SQUARE
	for (x = 1; x < NX - 1; x += TILE_COLS) {
		memcpy(tempBuffer[0][1], tIn + (x - 1), TEMP_COLS * sizeof(float));
		memcpy(tempBuffer[1][1], tIn + NY * NX + (x - 1), TEMP_COLS * sizeof(float));
		memcpy(tempBuffer[0][2], tIn + NX + (x - 1), TEMP_COLS * sizeof(float));
		memcpy(tempBuffer[1][2], tIn + NY * NX + NX + (x - 1), TEMP_COLS * sizeof(float));
		for (y = 1; y < NY - 1; y++) {
			// shift
			memcpy(tempBuffer[0][0], tempBuffer[0][1], TEMP_COLS * sizeof(float));
			memcpy(tempBuffer[1][0], tempBuffer[1][1], TEMP_COLS * sizeof(float));
			memcpy(tempBuffer[0][1], tempBuffer[0][2], TEMP_COLS * sizeof(float));
			memcpy(tempBuffer[1][1], tempBuffer[1][2], TEMP_COLS * sizeof(float));
			memcpy(tempBuffer[0][2], tIn + (y + 1) * NX + (x - 1), TEMP_COLS * sizeof(float));
			memcpy(tempBuffer[1][2], tIn + NY * NX + (y + 1) * NX + (x - 1), TEMP_COLS * sizeof(float));
			// load power
			memcpy(powerBuffer, pIn + y * NX + x, TILE_COLS * sizeof(float));
			// compute
			int i;
			for (i = 1; i < TEMP_COLS - 1; i++) {
				resultBuffer[i - 1] = tempBuffer[0][1][i] * cc + tempBuffer[0][1][i - 1] * cw + tempBuffer[0][1][i + 1] * ce
					+ tempBuffer[0][0][i] * cn + tempBuffer[0][2][i] * cs + tempBuffer[0][1][i] * cb + tempBuffer[1][1][i] * ct
					+ (dt/Cap) * powerBuffer[i - 1] + ct * AMB_TEMP;
			}
			// store
			memcpy(tOut + y * NX + x, resultBuffer, TILE_COLS * sizeof(float));
		}
	}

	// BOTTOM EDGE
	for (x = 0; x < NX; x++) {
		const int c = (NY - 1) * NX + x;
		float center = tIn[c];
		float west = (x == 0)		? center : tIn[c - 1];
		float east = (x == NX - 1) 	? center : tIn[c + 1];
		tOut[c] = center * cc + west * cw + east * ce + tIn[c - NX] * cn + center * cs + center * cb + tIn[c + NY * NX] * ct
			+ (dt/Cap) * pIn[c] + ct * AMB_TEMP; 
	}
}

void computeNorthFace(float* pIn, float* tIn, float* tOut, float tempBuffer[3][3][TEMP_COLS], float* powerBuffer, float* resultBuffer, float ce, float cw, float cn, float cs, float ct, float cb, float cc, float dt, float Cap)
{
	int x,y,z;
	for (x = 1; x < NX - 1; x += TILE_COLS) {
		// Initialize first two layers, accounting for a future shift
		for (z = 0; z < 2; z++) {
			for (y = 0; y < 2; y++) {
				memcpy(tempBuffer[z + 1][y], tIn + z * NY * NX + y * NX + (x - 1), TEMP_COLS * sizeof(float));
			}
		}
		
		for (z = 1; z < NZ - 1; z++) {
			// shift
			memcpy(tempBuffer[0], tempBuffer[1], 3 * TEMP_COLS * sizeof(float));
			memcpy(tempBuffer[1], tempBuffer[2], 3 * TEMP_COLS * sizeof(float));
			for (y = 0; y < 2; y++) {
				memcpy(tempBuffer[2][y], tIn + (z + 1) * NY * NX + y * NX + (x - 1), TEMP_COLS * sizeof(float));
			}
			// load power
			memcpy(powerBuffer, pIn + z * NY * NX + x, TILE_COLS * sizeof(float));
			// compute
			int i;
			for (i = 1; i < TEMP_COLS - 1; i++) {
				resultBuffer[i - 1] = tempBuffer[1][0][i] * cc + tempBuffer[1][0][i - 1] * cw + tempBuffer[1][0][i + 1] * ce
					+ tempBuffer[1][0][i] * cn + tempBuffer[1][1][i] * cs + tempBuffer[0][0][i] * cb + tempBuffer[2][0][i] * ct
					+ (dt/Cap) * powerBuffer[i - 1] + ct * AMB_TEMP;
			}
			// store
			memcpy(tOut + z * NY * NX + x, resultBuffer, TILE_COLS * sizeof(float));
		}
	}
}

void computeSouthFace(float* pIn, float* tIn, float* tOut, float tempBuffer[3][3][TEMP_COLS], float* powerBuffer, float* resultBuffer, float ce, float cw, float cn, float cs, float ct, float cb, float cc, float dt, float Cap)
{
	int x,y,z;
	for (x = 1; x < NX - 1; x += TILE_COLS) {
		// Initialize first two layers, accounting for a future shift
		for (z = 0; z < 2; z++) {
			for (y = 0; y < 2; y++) {
				const int row = (NY - 2) + y;
				memcpy(tempBuffer[z + 1][y], tIn + z * NY * NX + row * NX + (x - 1), TEMP_COLS * sizeof(float));
			}
		}
		
		for (z = 1; z < NZ - 1; z++) {
			// shift
			memcpy(tempBuffer[0], tempBuffer[1], 3 * TEMP_COLS * sizeof(float));
			memcpy(tempBuffer[1], tempBuffer[2], 3 * TEMP_COLS * sizeof(float));
			for (y = 0; y < 2; y++) {
				const int row = (NY - 2) + y;
				memcpy(tempBuffer[2][y], tIn + (z + 1) * NY * NX + row * NX + (x - 1), TEMP_COLS * sizeof(float));
			}
			// load power
			memcpy(powerBuffer, pIn + z * NY * NX + (NY - 1) * NX + x, TILE_COLS * sizeof(float));
			// compute
			int i;
			for (i = 1; i < TEMP_COLS - 1; i++) {
				resultBuffer[i - 1] = tempBuffer[1][1][i] * cc + tempBuffer[1][1][i - 1] * cw + tempBuffer[1][1][i + 1] * ce
					+ tempBuffer[1][0][i] * cn + tempBuffer[1][1][i] * cs + tempBuffer[0][1][i] * cb + tempBuffer[2][1][i] * ct
					+ (dt/Cap) * powerBuffer[i - 1] + ct * AMB_TEMP;
			}
			// store
			memcpy(tOut + z * NY * NX + (NY - 1) * NX + x, resultBuffer, TILE_COLS * sizeof(float));
		}
	}
}

void computeEastWestFaces(float *pIn, float* tIn, float *tOut, const float ce, const float cw, const float cn, const float cs, const float ct, const float cb, const float cc, float Cap, float dt) 
{   
	int z;
	for (z = 1; z < NZ - 1; z++) {
		// EAST FACE
		computeVerticalStrip(pIn, tIn, tOut, ce, cw, cn, cs, ct, cb, cc, dt, Cap, NX - 1, z);
		// WEST FACE
		computeVerticalStrip(pIn, tIn, tOut, ce, cw, cn, cs, ct, cb,cc, dt, Cap, 0, z);
	}
}

void initTempBuffer(float tempBuffer[3][3][TEMP_COLS], float* tIn, const int tileX)
{
	int z,y;
	for (z = 0; z < 3; z++) {
		for (y = 0; y < 2; y++) {
			memcpy_wide_bus_read_float(tempBuffer[z][y + 1], (MARS_WIDE_BUS_TYPE*)(tIn + z * NY * NX + y * NX + tileX - 1), 0, TEMP_COLS * sizeof(float));
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
		memcpy_wide_bus_read_float(tempBuffer[z][2], (MARS_WIDE_BUS_TYPE*)(tIn + layer * NY * NX + (y + 1) * NX + tileX - 1), 0, TEMP_COLS * sizeof(float));
	}
}

void loadPower(float* powerBuffer, float* pIn, const int y, const int tileX, const int z, const Direction direction)
{
	const int layer = (direction == FRONT_TO_BACK) ? z : NZ - 1 - z;
	memcpy_wide_bus_read_float(powerBuffer, (MARS_WIDE_BUS_TYPE*)(pIn + layer * NY * NX + y * NX + tileX), 0, TILE_COLS * sizeof(float));
}

void shiftLayers(float tempBuffer[3][3][TEMP_COLS], float* tIn, const int tileX, const int z, const int startRow, const Direction direction)
{
	int y;

	if (direction == FRONT_TO_BACK) {
		memcpy(tempBuffer[0], tempBuffer[1], 3 * TEMP_COLS * sizeof(float));
		memcpy(tempBuffer[1], tempBuffer[2], 3 * TEMP_COLS * sizeof(float));
		for (y = 0; y < 3; y++) {
			memcpy_wide_bus_read_float(tempBuffer[2][y], (MARS_WIDE_BUS_TYPE*)(tIn + z * NY * NX + (y + startRow) * NX + tileX - 1), 0, TEMP_COLS * sizeof(float));
		}
	} else {
		memcpy(tempBuffer[2], tempBuffer[1], 3 * TEMP_COLS * sizeof(float));
		memcpy(tempBuffer[1], tempBuffer[0], 3 * TEMP_COLS * sizeof(float));
		for (y = 0; y < 3; y++) {
			memcpy_wide_bus_read_float(tempBuffer[0][y], (MARS_WIDE_BUS_TYPE*)(tIn + (NZ - 1 - z) * NY * NX + (y + startRow) * NX + tileX - 1), 0, TEMP_COLS * sizeof(float));
		}
	}
}

void compute(float tempBuffer[3][3][TEMP_COLS], float* powerBuffer, float* resultBuffer, const float ce, const float cw, const float cn, const float cs, const float ct, const float cb, const float cc, const float dt, const float Cap)
{
COMPUTE: for (int x = 1; x < TEMP_COLS - 1; x++) {
		#pragma HLS unroll factor=5
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
	memcpy_wide_bus_write_float((MARS_WIDE_BUS_TYPE*)(tOut + layer * NY * NX + y * NX + tileX), resultBuffer, 0, TILE_COLS * sizeof(float));
}

void computeTempBorder(float* pIn, float* tIn, float* tOut, float ce, float cw, float cn, float cs, float ct, float cb, float cc, float dt, float Cap)
{
    int c,w,e,n,s,b,t;
    int x,y,z;
	for(z = 0; z < NZ; z++) 
    	for(y = 0; y < NY; y++)
        	for(x = 0; x < NX; x++)
            {
				if (!(z == 0 || z == NZ - 1 || x == 0 || x == NX - 1 || y == 0 || y == NY -1)) {
					continue;
				}

            	c = x + y * NX + z * NX * NY;

                w = (x == 0) ? c      : c - 1;
                e = (x == NX - 1) ? c : c + 1;
                n = (y == 0) ? c      : c - NX;
                s = (y == NY - 1) ? c : c + NX;
                b = (z == 0) ? c      : c - NX * NY;
                t = (z == NZ - 1) ? c : c + NX * NY;


                tOut[c] = tIn[c]*cc + tIn[n]*cn + tIn[s]*cs + tIn[e]*ce + tIn[w]*cw + tIn[t]*ct + tIn[b]*cb + (dt/Cap) * pIn[c] + ct*AMB_TEMP;
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
	float tempBuffer[3][3][TEMP_COLS];
	#pragma HLS ARRAY_PARTITION complete dim=1 variable=tempBuffer
	#pragma HLS ARRAY_PARTITION complete dim=2 variable=tempBuffer
	#pragma HLS ARRAY_PARTITION cyclic factor=5 dim=3 variable=tempBuffer
	float powerBuffer[TILE_COLS];
	#pragma HLS ARRAY_PARTITION cyclic factor=8 dim=1 variable=powerBuffer
	float resultBuffer[TILE_COLS];
	#pragma HLS ARRAY_PARTITION cyclic factor=8 dim=1 variable=resultBuffer

	Direction direction;
	
    do {
		//computeBottomFace(pIn, tIn, tOut, tempBuffer, powerBuffer, resultBuffer, ce, cw, cn, cs, ct, cb, cc, dt, Cap);
		//computeTopFace(pIn, tIn, tOut, tempBuffer, powerBuffer, resultBuffer, ce, cw, cn, cs, ct, cb, cc, dt, Cap);
		//computeNorthFace(pIn, tIn, tOut, tempBuffer, powerBuffer, resultBuffer, ce, cw, cn, cs, ct, cb, cc, dt, Cap);
		//computeSouthFace(pIn, tIn, tOut, tempBuffer, powerBuffer, resultBuffer, ce, cw, cn, cs, ct, cb, cc, dt, Cap);
		//computeEastWestFaces(pIn, tIn, tOut, ce, cw, cn, cs, ct, cb, cc, Cap, dt);
		computeTempBorder(pIn, tIn, tOut, ce, cw, cn, cs, ct, cb, cc, dt, Cap);

		for (x = 1; x < NX - 1; x += TILE_COLS) {
			initTempBuffer(tempBuffer, tIn, x);
			for (y = 1; y < NY - 1; y++) {
				direction = (y % 2) == 1 ? FRONT_TO_BACK : BACK_TO_FRONT;
				shiftRows(tempBuffer, tIn, x, y, direction);
				for (z = 2; z < NZ; z++) {
					if (z > 2) {
						shiftLayers(tempBuffer, tIn, x, z, y - 1, direction);
					}
					loadPower(powerBuffer, pIn, y, x, z - 1, direction);
					compute(tempBuffer, powerBuffer, resultBuffer, ce, cw, cn, cs, ct, cb, cc, dt, Cap);
					store(tOut, resultBuffer, x, y, z - 1, direction);
				}
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
