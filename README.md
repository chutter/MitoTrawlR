# MitoCap
An R package for processing high-throughput sequence data from raw reads to alignments from mitochondrial genomes sequenced as by-catch

The R package can be run as a single pipeline with a configuration file or through separate functions to create a customized pipeline. The R package has functions for: 

1) Iteratively assemble mitochondrial genomes from processed reads using a reference (mitochondrial genome or gene)
2) Annotate mitochondrial genomes based on reference
3) Align different features from the mitochondrial genome (rRNA, tRNA, CDS)
4) Trim alignments and concatenate for tree building 
5) Create finished mitochondrial genomes ready for submission to GenBank or other repositories


# Installation of prerequisites 

MitoCap uses several R packages and other outside programs for certain functions. 

1. R version 4.0.2 (and above likely work)

Requires:
- devtools
- ape
- stringr
- data.table
- seqinr
- biostrings (Bioconductor)
- genomicranges (Bioconductor)

Imports:
- foreach
- doparallel
- rdrop2



2. Outside programs

- fastp: adaptor trimming and paired-end read merging
- bwa: read mapping
- spades: assembly
- BLAST: matching assembled contigs to targets, other utilities
- mafft: creating alignments
- trimal: trimming alignments
- IqTree: gene tree and concatenation trees
- GATK4: variant calling functions
- SamTools: variant calling and read mapping tools
- tRNAScan: annotate tRNAs

First, you will want to clone this repository to your computer to obtain the setup files. Or alternatively go to the green "Code" button in top right of this repository and select "download ZIP".

```bash
git clone https://github.com/chutter/MitoCap.git
```

Second, change your working directory in the terminal to the downloaded repository. The key file here is the "MitoCap.yml" anaconda environment file, which must be present in the working directory being used. 

```bash
cd /MitoCap/setup-configuration_files/
```

The R packages and outside programs can be installed manually or more easily through the anaconda environment file provided (version numbers are provided in environment file if manual installation is desired). To install with the environment file, the easiest and quickest way is to first install the Anaconda package manager. Anaconda can be downloaded and installed for different operating systems from https://anaconda.org. Miniconda is recommended. Once installed, you can create a new environment for PhyloCap by: 

```bash
conda env create -f MitoCap.yml -n MitoCap
```

OR if a specific location for the environment directory is needed:

```bash
conda env create -f MitoCap.yml -p /PLACE/YOUR/DIRECTORY/HERE/MitoCap
```

And finally, you may delete the cloned GitHub directory after installing the prerequisites through the conda env file that manually installs the anaconda environment. There are some useful examples (also in the tutorial here), which could be saved.   

To use the environment, it must first be activated in your current terminal session or placed in your cluster job script. 

```bash
conda activate MitoCap
```

OR if a specific location for the environment directory is needed:

```bash
conda activate /PLACE/YOUR/DIRECTORY/HERE/MitoCap
```

# Installation of R package

The main functions of MitoCap are contained in an R package that has been tested on R version 4.0.2 and use the listed programs above along with custom scripts. To install PhyloCap from GitHub, you can use the R package devtools included in the environment above. When running in a cluster environment, the code for installation here should be included at the top of your R script with your selected PhyloCap functions. Here are step-by-step instructions for installation:

1) Install MitoCap by typing in your R console: 

```R
devtools::install_github("chutter/MitoCap", update = "never", dependencies = FALSE)
```

The update = "never" flag ensures that packages already installed via the anaconda environment are not changed, which will often break things. Additionally, dependencies = FALSE is set for the same reason. 


2) Devtools should finish and say the package loaded properly with no errors. Load the package in your R script with:

```R
library(MitoCap)
```

And installation should be done! All the functions for MitoCap should be ready to go! It is recommended to keep the install line above in your R script as the package is frequently updated for bugs and other features. In the future when there is a stable release, the R package will be available through Anaconda. 


3) You can run the following function to see if PhyloCap can find the dependencies: 

< coming soon a function to test if they can found >


# MitoCap pipeline tutorials 

[Installation: detailed installation instructions and trouble-shooting ](https://github.com/chutter/PhyloCap/wiki/Installation:-detailed-installation-instructions-and-trouble-shooting)

[Tutorial 1: PhyloCap configuration](https://github.com/chutter/PhyloCap/wiki/Tutorial-1:-PhyloCap-configuration)

[Tutorial 2: PhyloCap quick start pipeline](https://github.com/chutter/PhyloCap/wiki/Tutorial-2:-PhyloCap-quick-start-pipeline)

[Tutorial 3: Advanced function use](https://github.com/chutter/PhyloCap/wiki/Tutorial-3:-Advanced-function-use)

