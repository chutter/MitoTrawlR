#!/bin/bash
#SBATCH --job-name=data-pipe-frogs
#SBATCH --partition=bi
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=c111h652@ku.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=80gb
#SBATCH --time=168:00:00
#SBATCH -D /home/c111h652/scratch/MitoGenomes

module load anaconda/4.7

source activate /panfs/pfs.local/work/bi/c111h652/conda/frogcap

Rscript Run-Pipeline_Frogs_Dec14-2020.R
