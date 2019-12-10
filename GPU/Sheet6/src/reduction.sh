#!/usr/bin/env bash
#SBATCH --gres=gpu
#SBATCH -o out.txt
#SBATCH -w creek04

rm greduction.txt
for i in {8..17..1}
do
k=$((2**$i))
bin/reduction -s $k -t 64
done


for j in {3..10..1}
do
l=$((2**$j))
bin/reduction -s 8192 -t $l
done
