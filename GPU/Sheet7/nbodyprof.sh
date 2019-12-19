#!/usr/bin/env bash
#SBATCH --gres=gpu
#SBATCH -o out.txt
#SBATCH -w creek04

nvprof bin/nbody -s $s -i 100
nvprof bin/nbody -s $s -i 100 -opt
