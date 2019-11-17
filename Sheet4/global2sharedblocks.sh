#!/usr/bin/env bash
#SBATCH --gres=gpu
#SBATCH -o out.txt
rm global2sharedblocks.txt
for k in {0..10..1}
do
l=$((2**$k))
bin/memCpy -global2shared -i 10000 -s 10240 -t 1024 -g $l
done
