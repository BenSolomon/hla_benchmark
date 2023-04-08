#!/bin/sh
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=solomonb@stanford.edu
#SBATCH --time=13-23:05 # Runtime in D-HH:MM
#SBATCH --job-name=metaHLA
#SBATCH --nodes=1 # Ensure that all cores are reserved on one machine
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=25G
#SBATCH --partition=khatrilab # Partition allocated for the lab
#SBATCH --error=%x.err
#SBATCH --output=%x.out

module load R/4.0.4
Rscript /labs/khatrilab/solomonb/hla_project/hla_benchmark/7_HLA_ASE/sc_meta.R