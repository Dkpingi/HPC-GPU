#!/usr/bin/env bash
#SBATCH --gres=gpu
#SBATCH -o ex2_out.txt

rm Memory*.txt
for i in {0..20}
do
j=$((2**$i))
echo $j
bin/MemorySpeedD2D $j
done
