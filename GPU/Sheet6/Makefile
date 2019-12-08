NVCC=nvcc
INC=-I./inc -I/usr/local/cuda/include
LIB=-L/usr/local/cuda/lib64

.PHONY: build
build: ./bin/reduction

.PHONY: clean
clean:
	rm -f ./bin/*
	
.PHONY: rebuild
rebuild: clean build

./bin/reduction: ./src/main.cpp ./src/kernel.cu
	$(NVCC) -O2 -arch=sm_20 --compiler-options "-O2 -Wall -Wextra" -o $@ $^ $(INC) $(LIB)
