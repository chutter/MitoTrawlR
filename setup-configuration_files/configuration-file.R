#################################################
## Configuration file for PhyloCap
#################################################

#Directories and input files
#########################
# *** Full paths should be used whenever possible
#The main working directory
working.directory = "/directory/for/the/work"
#The input raw read directory
read.directory = "/directory/with/the/cleaned/reads"
#The name for the dataset
dataset.name = "Name"

#Reference configuration
#########################
#Provide a fasta file to a single sequence for the mitochondrial genome
reference.fasta = "reference.fa"
#provide a table or GFF file annotating the features (genes, tRNAs, etc) for reference.fasta
#GFF file can be downloaded directly from GenBank
annotation.file = "reference.gff"
#Only GFF available for now
annotation.type = "gff"

#Global settings
#########################
#number of threads
threads = 4
#Amount of memory to allocate in GB
memory = 8
#Whether to overwrite previous runs
overwrite = TRUE
#Print verbose output for each function
quiet = TRUE

#Mitogenome assembly configuration
#########################
#TRUE to iterative assemble; FALSE to assemble a single time
iterative.assemble = TRUE
#minimum number of iterations before terminating
min.iterations = 5
#maximum number of iterations before terminating
max.iterations = 30
#minimum allowable length in basepairs of a single contig before considering "complete"
min.length = 16000
#maximum allowable length in basepairs of a single contig before quitting
max.length = 40000
#minimum match proportion of the reads to the reference
min.read.match = 0.8

#Alignment and trimming settings
#############################
#TRUE = to run alignments for the matching targets from above
align.mitochondrial.markers = TRUE
#TRUE = to run alignment trimming function batchTrimAlignments
trim.alignments = TRUE
#The minimum number of taxa to keep an alignment
min.taxa.alignment = 4
#The minimum alignment basepairs to keep an alignment
min.alignment.length = 10
#The maximum gaps from throughout the entire alignment to keep an alignment
max.alignment.gap.percent = 50
#run the trimming program TAPER to trim out unalignment sample segments
run.TAPER = TRUE
#run the trimming program TrimAl to remove high variable or misaligned columns
run.TrimAl = TRUE
#Whether to trim out columns below a certain threshold
trim.column = TRUE
#The percent of bases that must be present to keep a column
min.column.gap.percent = 50
#Resolves ambiguous sites to the same arbitrary base
convert.ambiguous.sites = TRUE
#TRUE = to output an alignment assessment spreadsheet
alignment.assess = TRUE
#TRUE = to externally trim alignment edges
trim.external = TRUE
#The minimum percent of bases that must be present to keep a column on the edges
min.external.percent = 50
#TRUE = to trim samples below a certain coverage (percent bases present out of total alignment) threshold
trim.coverage = TRUE
#The minimum percent of bases that must be present to keep a sample
min.coverage.percent = 30
#The minimum number of bases that must be present to keep a sample
min.coverage.bp = 10


#Program paths
#########################
### *** Modify any of these from NULL to the path that the program is found if R is not detecting system paths
### e.g. fastp.path = "/conda/PhyloCap/bin
### e.g. fastp.path = "/conda/PhyloCap/bin
conda.env = "/panfs/pfs.local/work/bi/c111h652/conda/envs/MitoTrawlR/bin"
samtools.path = conda.env
bwa.path = conda.env
spades.path = conda.env
bbmap.path = conda.env
blast.path = conda.env
mafft.path = conda.env
iqtree.path = conda.env
trimAl.path = conda.env
tRNAscan.path = conda.env
cap3.path = conda.env

#### End configuration
