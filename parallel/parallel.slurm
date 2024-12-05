#!/bin/bash
#SBATCH --job-name=parallel
#SBATCH --output=outputs/parallel_250_vertices.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=1024
#SBATCH --time=00:05:00
#SBATCH --partition=normal

mpirun ./parallel ../grafos/grafo_250_vertices.txt