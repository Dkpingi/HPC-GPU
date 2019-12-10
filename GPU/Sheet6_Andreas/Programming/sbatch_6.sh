#!/usr/bin/env bash

#SBATCH -o 6_3_4.txt
#SBATCH --open-mode=append

# Iterate threads and problemsize
for t in 32 64 128 256 512 1024
do
	for i in 128 256 512 1024 2048 4096 8192 16384 32768 65536 131072
	do
	   ./bin/reduction -t $t -s $i
	done
done



