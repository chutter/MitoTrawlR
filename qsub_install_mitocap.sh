#!/bin/bash
#PBS -N install-mitocap
#PBS -l nodes=1:ppn=1,mem=10g,walltime=72:00:00
#PBS -M hutter@lsu.edu
#PBS -m e
#PBS -d /work/chutter/Programs
#PBS -j oe
#PBS -A hpc_jake_2020
#PBS -q single

mkdir conda

eval "$(conda shell.bash hook)"

conda create --prefix /work/chutter/Programs/conda/mitocap

conda activate /work/chutter/Programs/conda/mitocap

conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --set channel_priority flexible

conda install -c conda-forge r-base=3.5 r-devtools r-ape r-stringr r-data.table r-seqinr r-foreach r-doparallel r-rdrop2

conda install -c bioconda bioconductor-rsamtools bioconductor-genomicranges bioconductor-biostrings bioconductor-genbankr fastp spades=3.14.1 bbmap mafft trnascan-se bwa samtools=1.10 gatk4 trimal cap3 iqtree
