.PHONY: build
build: ./bin/nbody

.PHONY: clean
clean:
	rm ./bin/*
	
.PHONY: rebuild
rebuild: clean build

./bin/nbody: ./src/main.cu
	nvcc -O3 -arch=sm_20 -ftz=true --compiler-options "-O3 -Wextra" -o ./bin/nbody ./src/main.cu -I./inc/
