CUDA_ROOT=/opt/cuda-7.0
INC=-I./inc -I$(CUDA_ROOT)/include
LIB=-L$(CUDA_ROOT)/lib64
NVCC=$(CUDA_ROOT)/bin/nvcc
NVCC_FLAGS=-O2 -arch=sm_20 --compiler-options "-O2 -Wall -Wextra"

.PHONY: build
build: ./bin/MemorySpeedD2D

.PHONY: clean
clean:
	rm ./bin/*

.PHONY: rebuild
rebuild: clean build


./bin/MemorySpeedD2D: ./src/MemorySpeedD2D.cu
	$(NVCC) $(NVCC_FLAGS) -o $@ $^ $(INC) $(LIB)
