#!/bin/bash

mkdir conda

eval "$(conda shell.bash hook)"

conda create --name mitocap

conda activate mitocap

conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --set channel_priority flexible

conda install -c conda-forge r-base=4.0.2 r-devtools r-ape r-stringr r-data.table r-seqinr r-foreach r-doparallel r-rdrop2

conda install -c bioconda bioconductor-rsamtools bioconductor-genomicranges bioconductor-biostrings bioconductor-genbankr blast fastp spades=3.14.1 bbmap mafft trnascan-se bwa samtools=1.10 gatk4 trimal cap3 iqtree
