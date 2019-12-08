#!/bin/bash

# name of this programm
#SBATCH --job-name=test-mpi

# Total amount of nodes
#SBATCH --nodes=2

# Total amount of processes, spread over nodes
#SBATCH --ntasks=2

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

rm Ex1res.txt
for i in {0..20}
do
j=$((2**$i))
echo $j
mpirun ./Ex1.bin $j
done
