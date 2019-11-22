.PHONY: build
build: ./bin/matMul

.PHONY: clean
clean:
	rm ./bin/*
	
.PHONY: rebuild
rebuild: clean build

./bin/matMul: ./src/main.cu
	nvcc -arch=sm_20 -O2 --compiler-options "-O2 -Wall -Wextra" -o ./bin/matMul ./src/main.cu ./src/mmult_cpu.c -I./inc/
