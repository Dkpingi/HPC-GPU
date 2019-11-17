#!/usr/bin/env bash
#SBATCH --gres=gpu
#SBATCH -o out.txt
rm register2shared.txt
for j in {1024..49152..2048}
do
for k in {0..10..1}
do
for l in {0..10..1}
do
m=$((2**$k))
n=$((2**$l))
bin/memCpy -register2shared -i 1000 -s $j -t $m -g $n
done
done
done
