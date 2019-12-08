#!/usr/bin/env bash
#SBATCH --gres=gpu
#SBATCH -o out.txt

bin/matMul -shared -s 8192 -t 27
bin/matMul -s 8192 -t 9





