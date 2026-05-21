##############################################################################
## MitoTrawlR — Main Workflow Script
## Edit mitotrawlr_configuration-file.R first, then run this script
## end-to-end (or source individual sections as needed).
##############################################################################

# Install / update packages
devtools::install_github("chutter/PhyloProcessR", upgrade = "never", dependencies = FALSE)
devtools::install_github("chutter/MitoTrawlR",    upgrade = "never", dependencies = FALSE)

library(PhyloProcessR)
library(MitoTrawlR)

# Load configuration and set working directory
source("mitotrawlr_configuration-file.R")
setwd(working.directory)

# Create the analysis output directory and define path shortcuts
.adir = file.path(working.directory, analysis.dir)
dir.create(.adir, recursive = TRUE, showWarnings = FALSE)

.ref.dir    = file.path(.adir, "reference")
.contigs    = file.path(.adir, "draftContigs")
.annot      = file.path(.adir, "Annotations")
.align.un   = file.path(.adir, "Alignments", "untrimmed-alignments")
.align.tr   = file.path(.adir, "Alignments", "trimmed-alignments")
.genomes    = file.path(.adir, "MitoGenomes")
.phylo      = file.path(.adir, "Phylogeny")
.barcode    = file.path(.adir, "barcodeScan")


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
if (!pass.fail) stop("One or more required programs are missing. Check paths in the configuration file.")


##############################################################################
## Step 2 (OPTIONAL): Barcode scan — identify sample taxonomy before assembly
##############################################################################

if (run.barcode.scan) {
  MitoTrawlR::barcodeSampleScan(
    input.reads      = read.directory,
    barcode.fasta    = barcode.fasta,
    output.directory = .barcode,
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
  reference.name  = .ref.dir,
  overwrite       = overwrite,
  rep.origin      = rep.origin
)


##############################################################################
## Step 4: Iterative mitochondrial genome assembly
##############################################################################

MitoTrawlR::mitochondrialCapture(
  input.reads    = read.directory,
  reference.name = .ref.dir,
  output.dir     = .contigs,
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
  contig.folder  = .contigs,
  reference.name = .ref.dir,
  output.dir     = .annot,
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
  input.folder   = file.path(.annot, "sample-markers"),
  reference.name = .ref.dir,
  output.dir     = file.path(.adir, "Alignments"),
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
    alignment.dir            = .align.un,
    alignment.format         = "phylip",
    output.dir               = .align.tr,
    output.format            = "phylip",
    TrimAl                   = run.TrimAl,
    TrimAl.path              = trimAl.path,
    trim.external            = trim.external,
    min.external.percent     = min.external.percent,
    trim.coverage            = trim.coverage,
    min.coverage.percent     = min.coverage.percent,
    min.coverage.bp          = min.coverage.bp,
    trim.column              = trim.column,
    min.column.gap.percent   = min.column.gap.percent,
    convert.ambiguous.sites  = convert.ambiguous.sites,
    alignment.assess         = alignment.assess,
    min.alignment.length     = min.alignment.length,
    min.taxa.alignment       = min.taxa.alignment,
    max.alignment.gap.percent = max.alignment.gap.percent,
    threads                  = threads,
    memory                   = memory,
    overwrite                = overwrite
  )
}


##############################################################################
## Step 8: Concatenate per-marker alignments into whole-mitogenome alignment
##   Run on untrimmed, and again on trimmed if trimming was done.
##############################################################################

MitoTrawlR::alignMitogenomes(
  alignment.folder = .align.un,
  draft.contigs    = .contigs,
  reference.name   = .ref.dir,
  output.dir       = .genomes,
  dataset.name     = "untrimmed",
  mafft.path       = mafft.path,
  threads          = threads,
  overwrite        = overwrite
)

if (run.trim.alignments) {
  MitoTrawlR::alignMitogenomes(
    alignment.folder = .align.tr,
    draft.contigs    = .contigs,
    reference.name   = .ref.dir,
    output.dir       = .genomes,
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
    annotation.dir   = .annot,
    alignment.folder = .align.un,
    genome.alignment = file.path(.genomes, "alignments", "untrimmed_mitogenome_alignment.phy"),
    genome.dir       = .genomes,
    reference.name   = .ref.dir,
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

  .phylo.set = if (run.trim.alignments && phylo.dataset == "trimmed") "trimmed" else "untrimmed"
  if (.phylo.set == "trimmed" &&
      !file.exists(file.path(.genomes, "alignments", "trimmed_mitogenome_alignment.phy"))) {
    message("Trimmed alignment not found, falling back to untrimmed.")
    .phylo.set = "untrimmed"
  }

  MitoTrawlR::buildPhylogeny(
    alignment.file   = file.path(.genomes, "alignments", paste0(.phylo.set, "_mitogenome_alignment.phy")),
    feature.table    = file.path(.genomes, "alignments", paste0(.phylo.set, "_alignment_feature_table.txt")),
    output.dir       = .phylo,
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

  .phylo.set = if (run.trim.alignments && phylo.dataset == "trimmed") "trimmed" else "untrimmed"
  if (.phylo.set == "trimmed" &&
      !file.exists(file.path(.genomes, "alignments", "trimmed_mitogenome_alignment.phy"))) {
    message("Trimmed alignment not found, falling back to untrimmed.")
    .phylo.set = "untrimmed"
  }

  MitoTrawlR::plotMitoGenomes(
    tree.file         = file.path(.phylo, paste0(dataset.name, ".treefile")),
    bp.table          = file.path(.genomes, "logs", "untrimmed_mito-alignment_bp-count.csv"),
    feature.table     = file.path(.genomes, "alignments", paste0(.phylo.set, "_alignment_feature_table.txt")),
    genome.dir        = .genomes,
    output.format     = figure.format,
    outgroup          = outgroup,
    show.marker.types = show.marker.types,
    overwrite         = TRUE
  )
  # Plot saved to: <analysis.dir>/MitoGenomes/figures/mitogenome_completeness.<format>
}


##############################################################################
## Step 12 (OPTIONAL): Plot genome comparison (synteny / rearrangement viewer)
##############################################################################

if (run.plot.genome.comparison) {
  MitoTrawlR::plotGenomeComparison(
    annotation.dir = .annot,
    genome.dir     = .genomes,
    output.format  = figure.format,
    color.by       = genome.color.by,
    show.links     = genome.show.links,
    sync.genomes   = genome.sync,
    label.genes    = genome.label.genes,
    overwrite      = TRUE
  )
  # Plot saved to: <analysis.dir>/MitoGenomes/figures/genome_comparison.<format>
}
