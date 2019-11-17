#!/usr/bin/env bash
#SBATCH --gres=gpu
#SBATCH -o out.txt
rm shared2register.txt
for j in {1024..49152..2048}
do
for k in {0..10..1}
do
for l in {0..10..1}
do
a=$((2**$k))
b=$((2**$l))
bin/memCpy -shared2register -i 1000 -s $j -t $a -g $b
done
done
done
