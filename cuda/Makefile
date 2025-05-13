all: computeTempGPU

# PASCAL GT1080 and Volts
CUDA_ARCH=-gencode arch=compute_61,code=sm_61 -gencode arch=compute_70,code=sm_70 -gencode arch=compute_80,code=sm_80

computeTempGPU: computeTempGPU.o
	/usr/local/cuda-11.0/bin/nvcc -ccbin g++ -m64 ${CUDA_ARCH} -o computeTempGPU computeTempGPU.o -std=c++11

computeTempGPU.o:
	/usr/local/cuda-11.0/bin/nvcc -ccbin g++ -m64 ${CUDA_ARCH} -o computeTempGPU.o -c computeTempGPU.cu -std=c++11

run: computeTempGPU
	./computeTempGPU 512 8 100

clean:
	rm -rf *.o computeTempGPU *~
