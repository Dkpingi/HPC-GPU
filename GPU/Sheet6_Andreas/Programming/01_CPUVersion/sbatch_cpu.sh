#!/usr/bin/env bash

#SBATCH -o 6_2_cpu.txt
#SBATCH --open-mode=append

# Iterate the problem size
for i in 128 256 512 1024 2048 4096 8192 16384 32768 65536 131072
do
   ./ReductionCPU $i
done

