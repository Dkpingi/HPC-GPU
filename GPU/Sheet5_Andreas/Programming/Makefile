CUDA_ROOT=/opt/cuda-7.0
INC=-I./inc -I. -I$(CUDA_ROOT)/include
LIB=-L$(CUDA_ROOT)/lib64
NVCC=$(CUDA_ROOT)/bin/nvcc
NVCC_FLAGS=-O2 -arch=sm_20 --compiler-options "-O2 -Wall -Wextra"

.PHONY: build
build: ./bin/matMul

.PHONY: clean
clean:
	rm ./bin/*
	
.PHONY: rebuild
rebuild: clean build

./bin/matMul: ./src/main.cu
	$(NVCC) $(NVCC_FLAGS) -o $@ ./src/main.cu ./src/mmult_cpu.c $(INC) $(LIB)