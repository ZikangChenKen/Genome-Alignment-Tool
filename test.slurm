#!/bin/bash

#SBATCH --job-name=comp416-regular
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --exclusive
#SBATCH --threads-per-core=1
#SBATCH --time=01:00:00
#SBATCH --mail-type=ALL
#SBATCH --export=ALL
#SBATCH --partition=scavenge
#SBATCH --output=/projects/comp416/BJYZ/test-result.out

cd /projects/comp416/BJYZ/tests
python3 test_index_parallel.py