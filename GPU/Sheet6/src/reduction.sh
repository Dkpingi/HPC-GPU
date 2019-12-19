#!/usr/bin/env bash
#SBATCH --gres=gpu
#SBATCH -o out.txt
#SBATCH -w creek04

rm greduction.txt
for j in {3..10..1}
do
for i in {8..17..1}
do
k=$((2**$i))
l=$((2**$j))
bin/reduction -s $k -t $l
done
done

