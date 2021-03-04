#PBS -N proc-mito
#PBS -l nodes=1:ppn=16,mem=250g,walltime=48:00:00
#PBS -M hutter@lsu.edu
#PBS -m e
#PBS -d /scratch/chutter/Murinae
#PBS -j oe
#PBS -A hpc_jake_2020
#PBS -q bigmem

eval "$(conda shell.bash hook)"

conda activate /work/chutter/Programs/conda/mitocap

Rscript Process_Mitogenomes.R
