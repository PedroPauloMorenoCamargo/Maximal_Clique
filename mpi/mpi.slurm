#!/bin/bash
#SBATCH --job-name=mpi
#SBATCH --output=outputs/mpi_250_vertices.txt
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem=1024
#SBATCH --time=00:15:00
#SBATCH --partition=normal

mpirun ./mpi ../grafos/grafo_250_vertices.txt