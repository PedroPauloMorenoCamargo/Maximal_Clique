#!/bin/bash
#SBATCH --job-name=sequencial
#SBATCH --output=sequencial.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1024
#SBATCH --time=00:05:00
#SBATCH --partition=normal

mpirun ./sequential grafo_200_vertices.txt