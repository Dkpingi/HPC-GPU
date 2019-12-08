#!/usr/bin/env bash

# 5.3, Task 1: Determine optimal thread per block count
#SBATCH -o 5_3_01_benchmark.txt
#SBATCH --open-mode=append


# We use two different problem sizes
# Vary the thread size to determine an optimal thread count
# Args -p pinned mem -s matrix width -t threads per block -c do not check results
for i in {2..32..2}
do
	./bin/matMul  -p -s 1024 -t $i -c --shared
done

for i in {2..32..2}
do
	./bin/matMul  -p -s 2048 -t $i -c --shared
done

