##############################################################################
## MitoTrawlR — Configuration File
## Edit the values below, then run mitoGenome-workflow.R
## Full paths should be used wherever possible
##############################################################################


##############################################################################
## 1. Directories and dataset
##############################################################################

# Working directory — top-level project folder
working.directory = "/path/to/your/project"

# All pipeline outputs go inside this sub-directory of working.directory
# Change to suit your project layout; it will be created automatically
analysis.dir = "data-analysis/mitogenome"

# Folder of processed reads (one sub-directory per sample)
read.directory = "/path/to/processed-reads"

# Short label used in output file names (no spaces)
dataset.name = "myDataset"


##############################################################################
## 2. Reference genome
##############################################################################

# FASTA file of the reference mitochondrial genome
# Download from NCBI: Send to > File > Format: FASTA
reference.fasta = "reference_mitogenome.fasta"

# Annotation file for the reference (GFF3 recommended)
# Download from NCBI: Send to > File > Format: GFF3
annotation.file = "reference_mitogenome.gff3"

# Format of annotation.file: "gff", "genbank", or "table"
annotation.type = "gff"

# tRNAscan-SE organism model: "vertebrate", "mammal", or "eukaryotic"
organism.type = "vertebrate"

# Include the replication origin as a reference marker?
rep.origin = FALSE


##############################################################################
## 3. Global settings
##############################################################################

# Number of CPU threads
threads = 4

# RAM to allocate in GB
memory = 8

# Overwrite previous outputs? (FALSE = resume / skip completed samples)
overwrite = FALSE

# Suppress external tool screen output
quiet = TRUE


##############################################################################
## 4. Step toggles — set TRUE to run each step, FALSE to skip
##############################################################################

# Barcode scan: identify sample taxonomy via 16S before assembly
# Requires a separate barcode.fasta (see Section 5 below)
run.barcode.scan = FALSE

# Alignment trimming (trimMitoAlignments)
run.trim.alignments = TRUE

# Build polished final mitogenome FASTA files per sample (buildMitogenomes)
run.build.mitogenomes = TRUE

# Build maximum-likelihood phylogeny (buildPhylogeny / IQ-TREE)
run.phylogeny = TRUE

# Plot phylogeny + completeness heatmap (plotMitoGenomes)
run.plot.tree = TRUE

# Plot genome comparison / synteny viewer (plotGenomeComparison)
run.plot.genome.comparison = TRUE


##############################################################################
## 5. Barcode scan settings (used only when run.barcode.scan = TRUE)
##############################################################################

# FASTA file of the barcode locus (e.g. 16S) used to seed barcode assembly
barcode.fasta = "16S_barcode.fa"

# Source for taxonomic ID: "GenBank" (remote NCBI BLAST) or "File" (local db)
barcode.database = "GenBank"

# Path to a local BLAST database FASTA (only used when barcode.database = "File")
barcode.database.file = NULL

# Maximum GenBank hits to retrieve per sample
barcode.hits.per.sample = 5

# Maximum contig length as a fraction above the reference barcode length
barcode.per.max.length = 0.25


##############################################################################
## 6. Assembly settings (mitochondrialCapture)
##############################################################################

# Minimum assembly iterations before checking convergence
min.iterations = 5

# Maximum assembly iterations
max.iterations = 30

# Minimum assembled length (bp) to consider a contig plausible
min.length = 16000

# Maximum assembled length (bp) — contigs longer than this trigger stricter filtering
max.length = 40000

# Minimum read-to-reference identity for BBMap read recruitment (0–1)
min.read.match = 0.80


##############################################################################
## 7. Annotation settings (annotateMitoContigs)
##############################################################################

# Write a GFF3 annotation file per sample alongside the CSV?
save.gff = TRUE

# Minimum BLAST threads for annotation
# (uses global threads above; override here if needed)
# annot.threads = threads


##############################################################################
## 8. Alignment settings (markerAlignment / alignMitogenomes)
##############################################################################

# Maximum pairwise distance from reference to keep a sample sequence (0–1)
max.alignment.distance = 0.40


##############################################################################
## 9. Trimming settings (trimMitoAlignments; used when run.trim.alignments = TRUE)
##############################################################################

# Run TrimAl automated column trimming
run.TrimAl = TRUE

# Trim leading/trailing gap columns
trim.external = TRUE

# Minimum percent of samples that must have sequence to keep an edge column
min.external.percent = 50

# Remove samples with too little coverage
trim.coverage = TRUE

# Minimum percent of alignment columns a sample must cover to be kept
min.coverage.percent = 50

# Minimum number of unambiguous base pairs required per sample
min.coverage.bp = 0

# Remove alignment columns with too many gaps
trim.column = TRUE

# Maximum gap percent allowed in a column before it is removed
min.column.gap.percent = 100

# Resolve ambiguous IUPAC bases to the consensus character
convert.ambiguous.sites = FALSE

# Apply final quality assessment filters (discard failing alignments)
alignment.assess = TRUE

# Minimum alignment length (bp) to keep
min.alignment.length = 100

# Minimum number of taxa to keep an alignment
min.taxa.alignment = 4

# Maximum overall gap percent allowed for an alignment to pass assessment
max.alignment.gap.percent = 50


##############################################################################
## 10. Phylogeny settings (buildPhylogeny; used when run.phylogeny = TRUE)
##############################################################################

# Use trimmed or untrimmed alignment for the phylogeny?
# "trimmed" requires run.trim.alignments = TRUE; falls back to "untrimmed"
phylo.dataset = "trimmed"

# IQ-TREE ultrafast bootstrap replicates
bootstrap = 1000

# Partition scheme: "codon" (CDS pos1+2 / pos3), "byMarker", or "auto"
partition.scheme = "codon"


##############################################################################
## 11. Figure settings
##############################################################################

# Outgroup sample name(s) for rooting the tree in plotMitoGenomes
# Set to NULL to leave the tree unrooted
outgroup = NULL

# Output format for figures: "pdf" or "png"
figure.format = "pdf"

# Colour gene arrows by feature "type" (CDS/tRNA/rRNA/D-loop) or by "name"
genome.color.by = "type"

# Draw synteny ribbon links between adjacent genomes?
genome.show.links = TRUE

# Auto-flip / shift genomes to minimise link crossings?
genome.sync = TRUE

# Print gene name labels on arrows?
genome.label.genes = FALSE

# Show marker-type strip above the completeness heatmap?
show.marker.types = TRUE


##############################################################################
## 12. Program paths
##############################################################################
# Set conda.env to the root of your MitoTrawlR conda environment and all
# tool paths will be set automatically. Override individual paths below if
# a tool is installed outside the conda environment.

conda.env    = "/path/to/conda/envs/MitoTrawlR"

fastp.path   = paste0(conda.env, "/bin")
samtools.path = paste0(conda.env, "/bin")
bwa.path     = paste0(conda.env, "/bin")
spades.path  = paste0(conda.env, "/bin")
bbmap.path   = paste0(conda.env, "/bin")
blast.path   = paste0(conda.env, "/bin")
cap3.path    = paste0(conda.env, "/bin")
mafft.path   = paste0(conda.env, "/bin")
iqtree.path  = paste0(conda.env, "/bin")
trimAl.path  = paste0(conda.env, "/bin")
tRNAscan.path = paste0(conda.env, "/bin")

#### End configuration
