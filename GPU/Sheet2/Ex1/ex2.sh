#!/usr/bin/env bash
#SBATCH --gres=gpu
#SBATCH -o ex2_out.txt

rm SyncTime.txt AsyncTime.txt
for i in {1..10}; do
k=$((2**$i))
bin/nullKernelsync $k 1
done
for j in {1..10}; do
l=$((2**$j))
bin/nullKernelsync 1 $l
done
