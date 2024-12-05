#!/bin/bash
#SBATCH --job-name=parallel
#SBATCH --output=parallel.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=1024
#SBATCH --time=00:05:00
#SBATCH --partition=normal

mpirun ./parallel grafo_150_vertices.txt