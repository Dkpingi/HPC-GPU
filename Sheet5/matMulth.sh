#!/usr/bin/env bash
#SBATCH --gres=gpu
#SBATCH -o out.txt
for i in {4..32..1}
do
bin/matMul -shared -no-check -s 8192 -t $i
done

for i in {4..32..1}
do
bin/matMul -no-check -s 8192 -t $i
done




