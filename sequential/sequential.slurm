#!/bin/bash
#SBATCH --job-name=sequencial
#SBATCH --output=outputs/sequencial_150_vertices.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1024
#SBATCH --time=00:15:00
#SBATCH --partition=normal

mpirun ./sequential ../grafos/grafo_150_vertices.txt