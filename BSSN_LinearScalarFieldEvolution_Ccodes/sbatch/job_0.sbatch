#!/bin/bash
#SBATCH --job-name=nonlinear_980_8_16_0.0
#SBATCH --output=/home/goncaloa/NonLinear/sinhw0.2/output/out/nonlinear_980_8_16_0.0%j.out
#SBATCH --error=/home/goncaloa/NonLinear/sinhw0.2/output/err/nonlinear_980_8_16_0.0%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --mail-user=goncalo.j.c.andrade@tecnico.ulisboa.pt
#SBATCH --mail-type=ALL
#SBATCH --time=72:00:00
#SBATCH --mem=32G

RUNPATH=/home/goncaloa/NonLinear/sinhw0.2 ; cd $RUNPATH
mkdir output/fields_980_8_16_0.0 ; cd output/fields_980_8_16_0.0 ; mkdir checkpoint_data
export OMP_NUM_THREADS=48
../exec 980 8 16 0.5 0 0.0 0.075 5 0.5 0.0
