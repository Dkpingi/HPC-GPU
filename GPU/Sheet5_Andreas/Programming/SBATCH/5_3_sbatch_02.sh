#!/usr/bin/env bash

# 5.3, Task 2: With optimal thread count vary problem size (matrix size)
#SBATCH -o 5_3_02_benchmark.txt
#SBATCH --open-mode=append

# Args -p pinned mem -s matrix width -t threads per block -c do not check results
for i in 16 32 64 128 256 512 1024 2048 4096 8192
do
	./bin/matMul -p -s $i -t 32 -c --shared
done

