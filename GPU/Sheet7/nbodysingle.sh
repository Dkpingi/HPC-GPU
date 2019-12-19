#!/usr/bin/env bash
#SBATCH --gres=gpu
#SBATCH -o out.txt
#SBATCH -w creek04



s=$((2**16))

bin/nbody -s $s -i 10 -t 192
bin/nbody -s $s -i 10 -opt -t 192
