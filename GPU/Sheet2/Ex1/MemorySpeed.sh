#!/usr/bin/env bash
#SBATCH --gres=gpu
#SBATCH -o ex2_out.txt

rm Memory*.txt
for i in {0..10}
do
j=$((4000*(3**$i)))
echo $j
bin/MemorySpeed $j
done
