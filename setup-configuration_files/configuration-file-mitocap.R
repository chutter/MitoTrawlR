#################################################
## Configuration file for PhyloCap
#################################################

#Directories and input files
#########################
# *** Full paths should be used whenever possible
#The main working directory
work.dir = "/Volumes/Rodents/Mitogenomes/Crocidura"
#The input raw read directory
read.dir = "/Volumes/Rodents/Mitogenomes/Crocidura/processed-reads/pe-merged-reads"
#The name for the dataset
dataset.name = "Crocidura"

#Reference configuration
#########################
#Provide a fasta file to a single sequence for the mitochondrial genome
reference.fasta = "Mus_mus.fa"
#provide a table or GFF file annotating the features (genes, tRNAs, etc) for reference.fasta
#GFF file can be downloaded directly from GenBank
annotation.file = "Mus_mus.gff"
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
#Resume from previous runs (does not overwrite)
resume = FALSE
#Print verbose output for each function
quiet = TRUE



#Mitogenome assembly configuration
#########################
#
# mitochondrialCapture(input.reads = input.reads,
#                      reference.name = "reference",
#                      output.dir = "draftContigs",
#                      min.iterations = 5,
#                      max.iterations = 30,
#                      min.length = 16000,
#                      max.length = 40000,
#                      min.ref.id = 0.8,
#                      memory = memory,
#                      threads = threads,
#                      spades.path = spades.path,
#                      bbmap.path = bbmap.path,
#                      resume = TRUE,
#                      overwrite = FALSE)


#Alignment settings
#########################
#TRUE = to run alignments for the matching targets from above
align.mitochondrial.markers = FALSE
#The minimum number of taxa to keep an alignment
min.taxa.alignment = 4

#Alignment trimming settings
#########################
#TRUE = to run alignment trimming function batchTrimAlignments
trim.alignments = TRUE
#The minimum number of taxa to keep an alignment
min.taxa.alignment.trim = 5
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


#Tree settings
#########################
#TRUE = to estimate gene trees for each alignment
estimate.gene.trees = FALSE
#The minimum number of taxa to keep a gene tree
min.taxa.tree = 4
#TRUE = Clean up iqtree files except the ML tree
cleanup.genetrees = TRUE #Only saves the ML tree; FALSE saves all IQTree files for each gene tree

#Program paths
#########################
### *** Modify any of these from NULL to the path that the program is found if R is not detecting system paths
### e.g. fastp.path = "/conda/PhyloCap/bin
samtools.path = "/Users/chutter/conda/PhyloCap/bin"
bwa.path = "/Users/chutter/conda/PhyloCap/bin"
spades.path = "/Users/chutter/conda/PhyloCap/bin"
bbmap.path = "/Users/chutter/conda/PhyloCap/bin"
blast.path = "/Users/chutter/conda/PhyloCap/bin"
mafft.path = "/Users/chutter/conda/PhyloCap/bin"
iqtree.path = "/Users/chutter/conda/PhyloCap/bin"
trimAl.path = "/Users/chutter/conda/PhyloCap/bin"
tRNAscan.path = "/Users/chutter/conda/PhyloCap/bin"
cap3.path = "/Users/chutter/conda/PhyloCap/bin"

#### End configuration
