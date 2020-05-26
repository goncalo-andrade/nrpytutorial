#!/bin/bash
#SBATCH --job-name=fields
#SBATCH --output=/home/goncaloa/StaticScalarField/output/out/fields_%j.out
#SBATCH --error=/home/goncaloa/StaticScalarField/output/err/fields_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --mail-user=goncalo.j.c.andrade@tecnico.ulisboa.pt
#SBATCH --mail-type=ALL
#SBATCH --time=72:00:00
#SBATCH --mem=32G

RUNPATH=/home/goncaloa/StaticScalarField/output ; cd $RUNPATH
mkdir symbolic_pi_1024_32_64 ; cd symbolic_pi_1024_32_64
export OMP_NUM_THREADS=24
../exec 1024 32 64 0.5 1 64 8 0
