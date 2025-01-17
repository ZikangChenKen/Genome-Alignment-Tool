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
#SBATCH --output=/projects/comp416/BJYZ/generate-sam-result.out

cd /projects/comp416/BJYZ
python src/read_mapper.py --fasta short_reads_reference_genome.fasta --fastq challenging_dataset_1.fastq  -o output.sam