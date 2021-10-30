# MitoCap
An R package for iterative assembly and annotating mitochondrial genomes from short-read sequence data. 

The R package can be run as a single pipeline with a configuration file or through separate functions to create a customized pipeline. In one command, the pipeline will perform: 

1a) Iteratively assemble mitochondrial genomes from processed reads using a reference (mitochondrial genome or gene) [Slower but more complete data]

OR

1b) Extract mitochondrial genome from assembled contigs if its present [Faster but more missing data]

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
- tRNAScan: annotate tRNAs

First, you will want to clone this repository to your computer to obtain the setup files. Or alternatively go to the green "Code" button in top right of this repository and select "download ZIP".Second, change your working directory in the terminal to the downloaded repository. The key file here is the "MitoCap.yml" anaconda environment file, which must be present in the working directory being used. 

```bash
git clone https://github.com/chutter/MitoCap.git
cd /MitoCap/setup-configuration_files/
```

The R packages and outside programs can be installed manually or more easily through the anaconda environment file provided (version numbers are provided in environment file if manual installation is desired). To install with the environment file, the easiest and quickest way is to first install the Anaconda package manager. Anaconda can be downloaded and installed for different operating systems from https://anaconda.org. Miniconda is recommended. Once installed, you can create a new environment for MitoCap by: 

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

The main functions of MitoCap are contained in an R package that has been tested on R version 4.0.2 and use the listed programs above along with custom scripts. To install MitoCap from GitHub, you can use the R package devtools included in the environment above. When running in a cluster environment, the code for installation here should be included at the top of your R script with your selected PhyloCap functions. Here are step-by-step instructions for installation:

1) Install MitoCap by typing in your R console or including at the top of your R script (its already included in the automatic pipeline script): 

```R
devtools::install_github("chutter/MitoCap", upgrade = "never", dependencies = FALSE)
```

The update = "never" flag ensures that packages already installed via the anaconda environment are not changed, which will often break things. Additionally, dependencies = FALSE is set for the same reason. 


2) Devtools should finish and say the package loaded properly with no errors. Load the package in your R script with:

```R
library(MitoCap)
```

And installation should be done! All the functions for MitoCap should be ready to go! It is recommended to keep the install line above in your R script as the package is frequently updated for bugs and other features. In the future when there is a stable release, the R package will be available through Anaconda. 


3) You can run the following function to see if MitoCap can find the dependencies: 

MitoCap::setupCheck(anaconda.environment =  "path/to/anaconda/environment/MitoCap)

And it should return the list of programs and whether they were found or not. You can use this information to either install the program manually, and you'll have to add in the new non-anaconda path into the configuration file. 
                                

# MitoCap pipeline tutorials 

[Tutorial 1: MitoCap configuration](https://github.com/chutter/MitoCap/wiki/Tutorial-1:-MitoCap-configuration)

[Tutorial 2: MitoCap quick start pipeline](https://github.com/chutter/MitoCap/wiki/Tutorial-2:-MitoCap-quick-start-pipeline)

[Tutorial 3: Advanced function use](https://github.com/chutter/MitoCap/wiki/Tutorial-3:-Advanced-function-use)

