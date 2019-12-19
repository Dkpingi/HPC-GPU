#!/usr/bin/env bash
#SBATCH --gres=gpu
#SBATCH -o out.txt
#SBATCH -w creek04
rm nBody_opt.txt
rm nBody.txt

for t in {32..1024..32}
do
s=$((2**16))
bin/nbody -s $s -i 100 -t $t
bin/nbody -s $s -i 100 -opt -t $t
done
