#include "xcl2.hpp"

#include "computeTempFPGA.h"
#include "../../common/my_timer.h"

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

float accuracy(float *arr1, float *arr2)
{
	const int len = NZ * NY * NX;
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

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cout << "Usage: " << argv[0] << " <XCLBIN File>" << std::endl;
        return EXIT_FAILURE;
    }

    std::string binaryFile = argv[1];

	size_t grid_size_bytes = sizeof(float) * NX * NY * NZ;
    cl_int err;
    cl::CommandQueue q;
    cl::Context context;
    cl::Kernel kernel;

    float dx = chip_height / NX;
    float dy = chip_width / NY;
    float dz = t_chip / NZ;
    float Cap = FACTOR_CHIP * SPEC_HEAT_SI * t_chip * dx * dy;
    float Rx = dy / (2.0 * K_SI * t_chip * dx);
    float Ry = dx / (2.0 * K_SI * t_chip * dy);
    float Rz = dz / (K_SI * dx * dy);
    float max_slope = MAX_PD / (FACTOR_CHIP * t_chip * SPEC_HEAT_SI);
    float dt = PRECISION / max_slope;

    float* powerIn = (float*)malloc(grid_size_bytes);
    float* tempIn = (float*)malloc(grid_size_bytes);
    //float* tempCopy = (float*)malloc(grid_size_bytes);
    float* tempOut = (float*)malloc(grid_size_bytes);
    //float* answer = (float*)malloc(grid_size_bytes);

	initializeArrays(powerIn, tempIn);

    // OPENCL HOST CODE AREA START
    auto devices = xcl::get_xil_devices();
    // read_binary_file() is a utility API which will load the binaryFile
    // and will return the pointer to file buffer.
    auto fileBuf = xcl::read_binary_file(binaryFile);
    cl::Program::Binaries bins{{fileBuf.data(), fileBuf.size()}};
    bool valid_device = false;
    for (unsigned int i = 0; i < devices.size(); i++) {
        auto device = devices[i];
        // Creating Context and Command Queue for selected Device
        OCL_CHECK(err, context = cl::Context(device, nullptr, nullptr, nullptr, &err));
        OCL_CHECK(err, q = cl::CommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE, &err));

        std::cout << "Trying to program device[" << i << "]: " << device.getInfo<CL_DEVICE_NAME>() << std::endl;
        cl::Program program(context, {device}, bins, nullptr, &err);
        if (err != CL_SUCCESS) {
            std::cout << "Failed to program device[" << i << "] with xclbin file!\n";
        } else {
            std::cout << "Device[" << i << "]: program successful!\n";
            OCL_CHECK(err, kernel = cl::Kernel(program, "computeTempFPGA", &err));
            valid_device = true;
            break; // we break because we found a valid device
        }
    }
    if (!valid_device) {
        std::cout << "Failed to program any device found, exit!\n";
        exit(EXIT_FAILURE);
    }

    // Allocate Buffer in Global Memory
    OCL_CHECK(err, cl::Buffer buffer_pIn(context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, grid_size_bytes, powerIn, &err));
    OCL_CHECK(err, cl::Buffer buffer_tIn(context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, grid_size_bytes, tempIn, &err));
    OCL_CHECK(err, cl::Buffer buffer_tOut(context, CL_MEM_USE_HOST_PTR | CL_MEM_WRITE_ONLY, grid_size_bytes, tempOut, &err));

    OCL_CHECK(err, err = kernel.setArg(0, buffer_pIn));
    OCL_CHECK(err, err = kernel.setArg(1, buffer_tIn));
    OCL_CHECK(err, err = kernel.setArg(2, buffer_tOut));
    OCL_CHECK(err, err = kernel.setArg(3, Cap));
    OCL_CHECK(err, err = kernel.setArg(4, Rx));
    OCL_CHECK(err, err = kernel.setArg(5, Ry));
    OCL_CHECK(err, err = kernel.setArg(6, Rz));
    OCL_CHECK(err, err = kernel.setArg(7, dt));
    OCL_CHECK(err, err = kernel.setArg(8, NUM_ITERATIONS));

    // Copy input data to device global memory
    OCL_CHECK(err, err = q.enqueueMigrateMemObjects({buffer_pIn, buffer_tIn, buffer_tOut}, 0 /* 0 means from host*/));
    q.finish();

    // Start timer
    timespec timer = tic();

    // Launch the Kernel
    OCL_CHECK(err, err = q.enqueueTask(kernel));
    q.finish();

    // Stop timer
    toc(&timer, "Kernel Execution Time (s)");

    // Copy Result from Device Global Memory to Host Local Memory
    OCL_CHECK(err, err = q.enqueueMigrateMemObjects({buffer_tOut}, CL_MIGRATE_MEM_OBJECT_HOST));
    q.finish();

    // OPENCL HOST CODE AREA END

	// Verify result using known-good CPU
	// computeTempCPU(powerIn, tempCopy, answer, NX, NY, NZ, Cap, Rx, Ry, Rz, dt, NUM_ITERATIONS);

    //float acc = accuracy(tempOut, answer);
    //printf("Accuracy: %e\n", acc);

	//free(answer);
	free(tempOut);
	free(tempIn);
	//free(tempCopy);
	free(powerIn);

    return 0;
}
