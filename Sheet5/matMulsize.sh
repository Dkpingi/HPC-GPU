#!/usr/bin/env bash
#SBATCH --gres=gpu
#SBATCH -o out.txt
rm matMul.txt
for i in {1..13..1}
do
j=$((2**$i))
bin/matMul -shared -no-check -s $j -t 27
done

for i in {1..13..1}
do
j=$((2**$i))
bin/matMul -no-check -s $j -t 9
done

