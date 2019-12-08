#!/usr/bin/env bash
#SBATCH --gres=gpu
#SBATCH -o ex2_out.txt

rm BreakEven.txt
for j in {0..10000..10}; do
bin/BreakEven $j
done
