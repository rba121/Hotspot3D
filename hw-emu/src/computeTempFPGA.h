#ifndef COMPUTE_TEMP_FPGA_H
#define COMPUTE_TEMP_FPGA_H

#define LARGE_BUS 512
#define MARS_WIDE_BUS_TYPE ap_uint<LARGE_BUS>

#include "ap_int.h"
#include "../../common/mars_wide_bus.h"

#define NX 66
#define NY 66
#define NZ 8

#define NUM_ITERATIONS 100

#define TILE_COLS 16

#define AMB_TEMP 80.0

#endif
