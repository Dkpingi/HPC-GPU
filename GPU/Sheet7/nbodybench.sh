#!/usr/bin/env bash
#SBATCH --gres=gpu
#SBATCH -o out.txt
#SBATCH -w creek04
rm nBody_opt.txt
rm nBody.txt

for i in {0..30..1}
do
s=$((2**$i))

bin/nbody -s $s -i 100 -t 192
bin/nbody -s $s -i 100 -opt -t 192
done
