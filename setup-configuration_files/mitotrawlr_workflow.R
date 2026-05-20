##############################################################################
## MitoTrawlR — Main Workflow Script
## Edit configuration-file.R first, then run this script end-to-end
## (or source individual sections as needed).
##############################################################################

# Install / update packages
devtools::install_github("chutter/PhyloProcessR", upgrade = "never", dependencies = FALSE)
devtools::install_github("chutter/MitoTrawlR",    upgrade = "never", dependencies = FALSE)

library(PhyloProcessR)
library(MitoTrawlR)

# Load configuration and set working directory
source("configuration-file.R")
setwd(working.directory)


##############################################################################
## Step 1: Check that all required tools are installed
##############################################################################

pass.fail = MitoTrawlR::setupCheck(
  fastp.path    = fastp.path,
  samtools.path = samtools.path,
  bwa.path      = bwa.path,
  spades.path   = spades.path,
  bbmap.path    = bbmap.path,
  blast.path    = blast.path,
  cap3.path     = cap3.path,
  mafft.path    = mafft.path,
  iqtree.path   = iqtree.path,
  trimAl.path   = trimAl.path,
  tRNAscan.path = tRNAscan.path
)
if (!pass.fail) stop("One or more required programs are missing. Check paths in configuration-file.R.")


##############################################################################
## Step 2 (OPTIONAL): Barcode scan — identify sample taxonomy before assembly
##############################################################################

if (run.barcode.scan) {
  MitoTrawlR::barcodeSampleScan(
    input.reads      = read.directory,
    barcode.fasta    = barcode.fasta,
    output.directory = "barcodeScan",
    barcode.database = barcode.database,
    database.file    = barcode.database.file,
    hits.per.sample  = barcode.hits.per.sample,
    per.max.length   = barcode.per.max.length,
    memory           = memory,
    threads          = threads,
    spades.path      = spades.path,
    bbmap.path       = bbmap.path,
    blast.path       = blast.path,
    cap3.path        = cap3.path,
    overwrite        = overwrite,
    quiet            = quiet
  )
}


##############################################################################
## Step 3: Build reference from annotated mitogenome
##############################################################################

MitoTrawlR::buildReference(
  reference.fasta = reference.fasta,
  annotation.file = annotation.file,
  annotation.type = annotation.type,
  reference.name  = "reference",
  overwrite       = overwrite,
  rep.origin      = rep.origin
)


##############################################################################
## Step 4: Iterative mitochondrial genome assembly
##############################################################################

MitoTrawlR::mitochondrialCapture(
  input.reads    = read.directory,
  reference.name = "reference",
  output.dir     = "draftContigs",
  min.iterations = min.iterations,
  max.iterations = max.iterations,
  min.length     = min.length,
  max.length     = max.length,
  min.ref.id     = min.read.match,
  memory         = memory,
  threads        = threads,
  spades.path    = spades.path,
  bbmap.path     = bbmap.path,
  cap3.path      = cap3.path,
  blast.path     = blast.path,
  overwrite      = overwrite,
  quiet          = quiet
)


##############################################################################
## Step 5: Annotate mitochondrial contigs
##############################################################################

MitoTrawlR::annotateMitoContigs(
  contig.folder  = "draftContigs",
  reference.name = "reference",
  blast.path     = blast.path,
  tRNAscan.path  = tRNAscan.path,
  cap3.path      = cap3.path,
  organism.type  = organism.type,
  threads        = threads,
  overwrite      = overwrite,
  save.gff       = save.gff,
  quiet          = quiet
)


##############################################################################
## Step 6: Align each marker across all samples
##############################################################################

MitoTrawlR::markerAlignment(
  input.folder   = "Annotations/sample-markers",
  reference.name = "reference",
  threads        = threads,
  mafft.path     = mafft.path,
  max.distance   = max.alignment.distance,
  overwrite      = overwrite
)


##############################################################################
## Step 7 (OPTIONAL): Trim alignments
##############################################################################

if (run.trim.alignments) {
  MitoTrawlR::trimMitoAlignments(
    alignment.dir          = "Alignments/untrimmed-alignments",
    alignment.format       = "phylip",
    output.dir             = "Alignments/trimmed-alignments",
    output.format          = "phylip",
    TrimAl                 = run.TrimAl,
    TrimAl.path            = trimAl.path,
    TAPER                  = run.TAPER,
    TAPER.path             = taper.path,
    julia.path             = julia.path,
    trim.external          = trim.external,
    min.external.percent   = min.external.percent,
    trim.coverage          = trim.coverage,
    min.coverage.percent   = min.coverage.percent,
    min.coverage.bp        = min.coverage.bp,
    trim.column            = trim.column,
    min.column.gap.percent = min.column.gap.percent,
    convert.ambiguous.sites = convert.ambiguous.sites,
    alignment.assess       = alignment.assess,
    min.alignment.length   = min.alignment.length,
    min.taxa.alignment     = min.taxa.alignment,
    max.alignment.gap.percent = max.alignment.gap.percent,
    threads                = threads,
    memory                 = memory,
    overwrite              = overwrite
  )
}


##############################################################################
## Step 8: Concatenate per-marker alignments into whole-mitogenome alignment
##   Run once on untrimmed, and again on trimmed if trimming was done.
##############################################################################

MitoTrawlR::alignMitogenomes(
  alignment.folder = "Alignments/untrimmed-alignments",
  draft.contigs    = "draftContigs",
  reference.name   = "reference",
  output.dir       = "MitoGenomes",
  dataset.name     = "untrimmed",
  mafft.path       = mafft.path,
  threads          = threads,
  overwrite        = overwrite
)

if (run.trim.alignments) {
  MitoTrawlR::alignMitogenomes(
    alignment.folder = "Alignments/trimmed-alignments",
    draft.contigs    = "draftContigs",
    reference.name   = "reference",
    output.dir       = "MitoGenomes",
    dataset.name     = "trimmed",
    mafft.path       = mafft.path,
    threads          = threads,
    overwrite        = FALSE
  )
}


##############################################################################
## Step 9 (OPTIONAL): Build polished final mitogenome FASTA files
##############################################################################

if (run.build.mitogenomes) {
  MitoTrawlR::buildMitogenomes(
    annotation.dir   = "Annotations",
    alignment.folder = "Alignments/untrimmed-alignments",
    genome.alignment = paste0("MitoGenomes/alignments/untrimmed_mitogenome_alignment.phy"),
    genome.dir       = "MitoGenomes",
    reference.name   = "reference",
    output.dir       = "untrimmed-finished",
    blast.path       = blast.path,
    threads          = threads,
    overwrite        = overwrite
  )
}


##############################################################################
## Step 10 (OPTIONAL): Build maximum-likelihood phylogeny with IQ-TREE
##############################################################################

if (run.phylogeny) {

  # Choose alignment: prefer trimmed if available, fall back to untrimmed
  phylo.align.dir = if (run.trim.alignments && phylo.dataset == "trimmed") "trimmed" else "untrimmed"

  MitoTrawlR::buildPhylogeny(
    alignment.file   = paste0("MitoGenomes/alignments/", phylo.align.dir, "_mitogenome_alignment.phy"),
    feature.table    = paste0("MitoGenomes/alignments/", phylo.align.dir, "_alignment_feature_table.txt"),
    output.dir       = "Phylogeny",
    dataset.name     = dataset.name,
    iqtree.path      = iqtree.path,
    threads          = threads,
    bootstrap        = bootstrap,
    partition.scheme = partition.scheme,
    overwrite        = overwrite
  )
}


##############################################################################
## Step 11 (OPTIONAL): Plot phylogeny + per-marker completeness heatmap
##############################################################################

if (run.plot.tree) {

  phylo.align.dir = if (run.trim.alignments && phylo.dataset == "trimmed") "trimmed" else "untrimmed"

  MitoTrawlR::plotMitoGenomes(
    tree.file         = paste0("Phylogeny/", dataset.name, ".treefile"),
    bp.table          = paste0("MitoGenomes/logs/untrimmed_mito-alignment_bp-count.csv"),
    feature.table     = paste0("MitoGenomes/alignments/", phylo.align.dir, "_alignment_feature_table.txt"),
    genome.dir        = "MitoGenomes",
    output.format     = figure.format,
    outgroup          = outgroup,
    show.marker.types = show.marker.types,
    overwrite         = TRUE
  )
  # Plot saved to: MitoGenomes/figures/mitogenome_completeness.<format>
}


##############################################################################
## Step 12 (OPTIONAL): Plot genome comparison (synteny / rearrangement viewer)
##############################################################################

if (run.plot.genome.comparison) {
  MitoTrawlR::plotGenomeComparison(
    annotation.dir = "Annotations",
    genome.dir     = "MitoGenomes",
    output.format  = figure.format,
    color.by       = genome.color.by,
    show.links     = genome.show.links,
    sync.genomes   = genome.sync,
    label.genes    = genome.label.genes,
    overwrite      = TRUE
  )
  # Plot saved to: MitoGenomes/figures/genome_comparison.<format>
}
