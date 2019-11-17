#!/usr/bin/env bash
#SBATCH --gres=gpu
#SBATCH -o ex2_out.txt

bin/memCpy
