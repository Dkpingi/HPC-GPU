#!/bin/bash

# name of this programm
#SBATCH --job-name=test-mpi

# Total amount of nodes
#SBATCH --nodes=2

# Runtime of this jobs is less then 15 minutes.
#SBATCH --time=0:30:00

# Name of output file
#SBATCH --output=out.txt

# Clear the environment from any previously loaded modules
module purge > /dev/null 2>&1

# Load the module environment suitable for the job
module load cmake
module load gcc
module load mpi

for i in {2..16..1}
do
for j in {1..100..1}
do
mpirun -np $i ./bin/relaxation $j 100
done
done

for i in {1,10,100,1000,10000,100000,1000000}
do
mpirun -np 8 ./bin/relaxation 100 $i
done


