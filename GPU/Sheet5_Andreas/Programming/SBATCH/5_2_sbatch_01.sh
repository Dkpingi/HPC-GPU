#!/usr/bin/env bash

# 5.2, Task 1: Determine optimal thread per block count
#SBATCH -o 5_2_01_benchmark.txt
#SBATCH --open-mode=append



# Vary the thread size to determine an optimal thread count
# Args -p pinned mem -s matrix width -t threads per block -c do not check results
for i in {2..32..2}
do
	./bin/matMul  -p -s 1024 -t $i -c
done

