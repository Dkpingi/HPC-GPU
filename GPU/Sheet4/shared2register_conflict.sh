#!/usr/bin/env bash
#SBATCH --gres=gpu
#SBATCH -o out.txt
rm shared2register_conflict.txt
for s in {0..64..1}
do
bin/memCpy -shared2register_conflict -i 10000 -s 49152 -t 192 -g 100 -stride $s
done
