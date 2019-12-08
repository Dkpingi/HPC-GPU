#!/usr/bin/env bash
#SBATCH --gres=gpu
#SBATCH -o out.txt
rm global2shared.txt
for j in {1024..49152..2048}
do
for k in {0..10..1}
do
l=$((2**$k))
bin/memCpy -global2shared -i 10000 -s $j -t $l -g 1
done
done
