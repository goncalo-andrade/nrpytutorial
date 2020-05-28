#!/bin/bash
#SBATCH --job-name=fields
#SBATCH --output=/home/goncaloa/StaticScalarField_Sommerfeld/output/out/fields_%j.out
#SBATCH --error=/home/goncaloa/StaticScalarField_Sommerfeld/output/err/fields_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --mail-user=goncalo.j.c.andrade@tecnico.ulisboa.pt
#SBATCH --mail-type=ALL
#SBATCH --time=72:00:00
#SBATCH --mem=32G

RUNPATH=/home/goncaloa/StaticScalarField_Sommerfeld/output ; cd $RUNPATH
mkdir sommerfeld_512_16_32 ; cd sommerfeld_512_16_32
export OMP_NUM_THREADS=24
../exec 512 32 64 0.5 1 64 8 0
