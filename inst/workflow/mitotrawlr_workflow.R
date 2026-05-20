################################################################################
# MitoTrawlR — complete end-to-end workflow
#
# Adjust the path variables in Step 0 before running. All steps after that
# are self-contained. Steps marked OPTIONAL can be skipped.
################################################################################

library(MitoTrawlR)

################################################################################
# Step 0: Paths — set once, used everywhere
################################################################################

# Root directory for this project (set working directory here)
# setwd("/path/to/project")

# Path to the conda environment created from setup-configuration_files/environment.yml
conda.env     = "/opt/miniconda3/envs/MitoTrawlR"

# Input reads folder: one sub-directory per sample containing fastq/fq files
reads.dir     = "processed-reads/adaptor-removed-reads"

# Reference name (folder created by buildReference)
ref.name      = "reference"

# Draft contigs output folder
contigs.dir   = "draftContigs"

# Annotation output folder
annot.dir     = "Annotations"

# Alignment output folder root
align.dir     = "Alignments"

# MitoGenomes output folder
mito.dir      = "MitoGenomes"

# Phylogeny output folder
phylo.dir     = "Phylogeny"

# Dataset label used for file names
dataset       = "untrimmed"

# Number of CPU threads
n.threads     = 4

# RAM in GB
n.memory      = 8


################################################################################
# Step 1: Check that all required tools are installed
################################################################################

all.found = MitoTrawlR::setupCheck(
  anaconda.environment = conda.env
)
if (!all.found) stop("One or more tools not found — check paths before continuing.")


################################################################################
# Step 2 (OPTIONAL): Barcode scan — identify samples by 16S before assembly
################################################################################

# MitoTrawlR::barcodeSampleScan(
#   input.reads      = reads.dir,
#   barcode.fasta    = "16S_barcode.fa",
#   output.directory = "barcodeScan",
#   barcode.database = "GenBank",
#   hits.per.sample  = 5,
#   memory           = n.memory,
#   threads          = n.threads,
#   overwrite        = FALSE
# )


################################################################################
# Step 3: Build reference from a GenBank accession or local FASTA
################################################################################

MitoTrawlR::buildReference(
  genbank.accession = NULL,          # e.g. "NC_012345" — or supply fasta.file
  fasta.file        = "reference_mitogenome.fa",
  reference.name    = ref.name,
  overwrite         = FALSE
)


################################################################################
# Step 4: Iterative mitochondrial genome assembly for all samples
################################################################################

MitoTrawlR::mitochondrialCapture(
  input.reads    = reads.dir,
  reference.name = ref.name,
  output.dir     = contigs.dir,
  min.iterations = 5,
  max.iterations = 20,
  min.length     = 15000,
  max.length     = 30000,
  min.ref.id     = 0.75,
  memory         = n.memory,
  threads        = n.threads,
  overwrite      = FALSE,
  quiet          = TRUE
)


################################################################################
# Step 5: Annotate mitochondrial contigs
################################################################################

MitoTrawlR::annotateMitoContigs(
  input.dir      = contigs.dir,
  reference.name = ref.name,
  output.dir     = annot.dir,
  threads        = n.threads,
  overwrite      = FALSE
)


################################################################################
# Step 6: Align each marker across all samples
################################################################################

MitoTrawlR::markerAlignment(
  input.folder   = paste0(annot.dir, "/sample-markers"),
  reference.name = ref.name,
  threads        = n.threads,
  max.distance   = 0.40,
  overwrite      = FALSE
)


################################################################################
# Step 7 (OPTIONAL): Trim alignments with trimAl / TAPER
################################################################################

# MitoTrawlR::trimMitoAlignments(
#   input.folder  = paste0(align.dir, "/untrimmed-alignments"),
#   output.folder = paste0(align.dir, "/trimmed-alignments"),
#   codon.trim    = FALSE,
#   overwrite     = FALSE
# )


################################################################################
# Step 8: Concatenate per-marker alignments into a whole-mitogenome alignment
################################################################################

MitoTrawlR::alignMitogenomes(
  alignment.folder = paste0(align.dir, "/untrimmed-alignments"),
  draft.contigs    = contigs.dir,
  reference.name   = ref.name,
  output.dir       = mito.dir,
  dataset.name     = dataset,
  threads          = n.threads,
  overwrite        = FALSE
)


################################################################################
# Step 9 (OPTIONAL): Build polished final mitogenome FASTA files per sample
################################################################################

# MitoTrawlR::buildMitogenomes(
#   annotation.dir   = annot.dir,
#   alignment.folder = paste0(align.dir, "/untrimmed-alignments"),
#   genome.alignment = paste0(mito.dir, "/alignments/", dataset, "_mitogenome_alignment.phy"),
#   genome.dir       = mito.dir,
#   reference.name   = ref.name,
#   output.dir       = "untrimmed-finished",
#   overwrite        = FALSE
# )


################################################################################
# Step 10: Build maximum-likelihood phylogeny with IQ-TREE
################################################################################

MitoTrawlR::buildPhylogeny(
  alignment.file   = paste0(mito.dir, "/alignments/", dataset, "_mitogenome_alignment.phy"),
  feature.table    = paste0(mito.dir, "/alignments/", dataset, "_alignment_feature_table.txt"),
  output.dir       = phylo.dir,
  dataset.name     = dataset,
  threads          = n.threads,
  bootstrap        = 1000,
  partition.scheme = "codon",
  overwrite        = FALSE
)


################################################################################
# Step 11: Plot phylogeny + completeness heatmap
################################################################################

MitoTrawlR::plotMitoGenomes(
  tree.file         = paste0(phylo.dir, "/", dataset, ".treefile"),
  bp.table          = paste0(mito.dir, "/logs/", dataset, "_mito-alignment_bp-count.csv"),
  feature.table     = paste0(mito.dir, "/alignments/", dataset, "_alignment_feature_table.txt"),
  genome.dir        = mito.dir,
  output.format     = "pdf",
  width             = 16,
  height            = 10,
  outgroup          = NULL,
  show.marker.types = TRUE,
  overwrite         = TRUE
)
# Saves to: MitoGenomes/figures/mitogenome_completeness.pdf


################################################################################
# Step 12: Plot genome comparison (synteny / rearrangement viewer)
################################################################################

MitoTrawlR::plotGenomeComparison(
  annotation.dir = annot.dir,
  genome.dir     = mito.dir,
  output.format  = "pdf",
  color.by       = "type",
  show.links     = TRUE,
  sync.genomes   = TRUE,
  width          = 14,
  overwrite      = TRUE
)
# Saves to: MitoGenomes/figures/genome_comparison.pdf
