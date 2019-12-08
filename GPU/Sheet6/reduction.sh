#!/usr/bin/env bash
#SBATCH --gres=gpu
#SBATCH -o out.txt
#SBATCH -w creek04

rm greduction.txt
bin/reduction -s 1048576 -t 1024


