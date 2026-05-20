#' @title plotGenomeComparison
#'
#' @description Plots multiple mitochondrial genome annotations stacked
#'   vertically using \pkg{gggenomes}, with gene arrows coloured by feature
#'   type or gene name. When \code{show.links = TRUE}, ribbon links connect
#'   identically-named genes between adjacent samples, making gene-order
#'   rearrangements immediately visible as crossing ribbons. When
#'   \code{sync.genomes = TRUE}, sequences are automatically flipped and
#'   shifted by \code{gggenomes::sync()} to minimise link crossings before
#'   plotting.
#'
#'   Multi-contig samples are handled by stacking contigs into a single linear
#'   coordinate space (each contig's coordinates are offset by the cumulative
#'   length of preceding contigs) so they appear as one continuous track.
#'
#' @param annotation.dir path to the annotations directory produced by
#'   \code{annotateMitoContigs}, which must contain a
#'   \code{sample-summary/} subdirectory of per-sample CSV files.
#'
#' @param contig.dir path to the folder of per-sample contig FASTA files used
#'   to obtain exact sequence lengths. Defaults to
#'   \code{annotation.dir/sample-contigs/}. If a sample's FASTA is absent,
#'   the length is estimated as \code{max(end)} across all its annotations.
#'
#' @param samples character vector of sample names to include. NULL includes
#'   all samples found in \code{annotation.dir/sample-summary/}.
#'
#' @param output.file path (without extension) for the saved plot.
#'
#' @param output.format output file format: \code{"pdf"} or \code{"png"}.
#'
#' @param width plot width in inches.
#'
#' @param height plot height in inches. NULL auto-scales with sample count.
#'
#' @param color.by how to colour gene arrows. \code{"type"} uses four colours
#'   for CDS, tRNA, rRNA, and D-loop; \code{"name"} colours each gene
#'   individually and adds short labels.
#'
#' @param show.links logical; if TRUE, draw ribbon links connecting the same
#'   gene between adjacent samples to reveal rearrangements.
#'
#' @param sync.genomes logical; if TRUE, call \code{gggenomes::sync()} to
#'   flip and shift sequences so that links cross as little as possible.
#'
#' @param label.genes logical; if TRUE, print gene names inside or above
#'   arrows (only when \code{color.by = "name"} or when \code{color.by =
#'   "type"} and there are few enough genes to label without overlap).
#'
#' @param overwrite logical; if TRUE, an existing output file is overwritten.
#'
#' @return Invisibly returns the gggenomes plot object. Writes a PDF or PNG
#'   to \code{output.file}.
#'
#' @examples
#' \dontrun{
#' plotGenomeComparison(
#'   annotation.dir = "Annotations",
#'   samples        = c("Sample_A", "Sample_B", "Sample_C"),
#'   output.file    = "mito_rearrangements",
#'   output.format  = "pdf",
#'   color.by       = "type",
#'   show.links     = TRUE,
#'   sync.genomes   = TRUE,
#'   overwrite      = TRUE
#' )
#' }
#'
#' @export

plotGenomeComparison = function(annotation.dir = "Annotations",
                                contig.dir = NULL,
                                genome.dir = "MitoGenomes",
                                samples = NULL,
                                output.file = NULL,
                                output.format = c("pdf", "png"),
                                width = 14,
                                height = NULL,
                                color.by = c("type", "name"),
                                show.links = TRUE,
                                sync.genomes = TRUE,
                                label.genes = FALSE,
                                overwrite = TRUE) {

  for (pkg in c("ggplot2", "gggenomes")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("Package '", pkg, "' is required for plotGenomeComparison(). Install it with: install.packages('", pkg, "')")
  }

  output.format = match.arg(output.format)
  color.by      = match.arg(color.by)

  if (is.null(output.file)) {
    figures.dir = paste0(genome.dir, "/figures")
    dir.create(figures.dir, recursive = TRUE, showWarnings = FALSE)
    output.file = paste0(figures.dir, "/genome_comparison")
  }

  if (is.null(contig.dir)) {
    contig.dir = paste0(annotation.dir, "/sample-contigs")
  }

  summary.dir = paste0(annotation.dir, "/sample-summary")
  if (!dir.exists(summary.dir)) {
    stop("Annotation summary directory not found: ", summary.dir)
  }

  out.path = paste0(output.file, ".", output.format)
  if (file.exists(out.path) && !overwrite) {
    stop("Output file exists and overwrite = FALSE: ", out.path)
  }

  ##############################################################################
  # Step 1: Discover samples
  ##############################################################################
  csv.files   = list.files(summary.dir, pattern = "_sample-summary\\.csv$",
                            full.names = TRUE)
  if (length(csv.files) == 0) {
    stop("No sample-summary CSV files found in: ", summary.dir)
  }
  all.samples = gsub("_sample-summary\\.csv$", "", basename(csv.files))

  if (!is.null(samples)) {
    missing = setdiff(samples, all.samples)
    if (length(missing) > 0) {
      warning("Samples not found in annotation directory and will be skipped: ",
              paste(missing, collapse = ", "))
    }
    all.samples = samples[samples %in% all.samples]
  }

  if (length(all.samples) == 0) { stop("No valid samples to plot.") }

  ##############################################################################
  # Step 2: Build seqs and genes data frames (one row in seqs per sample)
  ##############################################################################

  # Helper: strip the "NNN_TYPE_" prefix from reference marker names
  strip.prefix = function(x) gsub("^[0-9]+_[^_]+_", "", x)

  # Helper: classify stripped gene name into a feature type
  classify.type = function(x) {
    ifelse(grepl("^tRNA",           x, ignore.case = FALSE), "tRNA",
    ifelse(grepl("rRNA|12S|16S",   x, ignore.case = TRUE),  "rRNA",
    ifelse(grepl("D.loop|D_loop",  x, ignore.case = TRUE),  "D-loop",
           "CDS")))
  }

  seqs.list  = vector("list", length(all.samples))
  genes.list = vector("list", length(all.samples))

  for (i in seq_along(all.samples)) {
    spp = all.samples[i]

    ann = utils::read.csv(paste0(summary.dir, "/", spp, "_sample-summary.csv"),
                          stringsAsFactors = FALSE)

    if (nrow(ann) == 0) {
      warning(spp, ": annotation table is empty. Skipping.")
      next
    }

    ##########################################################################
    # Get per-contig lengths and compute cumulative offsets so all contigs
    # appear as a single linear sequence per sample.
    ##########################################################################
    fa.file = paste0(contig.dir, "/", spp, "_sampleContigs.fa")
    if (file.exists(fa.file)) {
      fa      = Biostrings::readDNAStringSet(fa.file)
      ctg.len = setNames(as.integer(Biostrings::width(fa)), names(fa))
    } else {
      # Estimate length from annotation endpoints
      ctg.len = tapply(ann$end, ann$contig, max, na.rm = TRUE)
      ctg.len = as.integer(ctg.len)
      names(ctg.len) = names(tapply(ann$end, ann$contig, max, na.rm = TRUE))
    }

    # Sort contigs in natural order (seq1, seq2, …)
    ctg.order  = names(ctg.len)[order(as.integer(gsub("[^0-9]", "", names(ctg.len))))]
    ctg.len    = ctg.len[ctg.order]
    ctg.offset = c(0L, cumsum(as.integer(ctg.len[-length(ctg.len)])))
    names(ctg.offset) = ctg.order
    total.len  = sum(ctg.len)

    seqs.list[[i]] = data.frame(seq_id = spp,
                                length = total.len,
                                stringsAsFactors = FALSE)

    # Apply offsets to annotation coordinates
    offset.vec = ctg.offset[ann$contig]
    offset.vec[is.na(offset.vec)] = 0L

    gene.name = strip.prefix(ann$name)
    gene.type = classify.type(gene.name)

    genes.list[[i]] = data.frame(
      seq_id    = spp,
      start     = ann$start + offset.vec,
      end       = ann$end   + offset.vec,
      strand    = ann$direction,
      name      = gene.name,
      type      = gene.type,
      stringsAsFactors = FALSE
    )
  }

  # Drop any NULLs from skipped samples
  seqs.list  = Filter(Negate(is.null), seqs.list)
  genes.list = Filter(Negate(is.null), genes.list)

  if (length(seqs.list) == 0) { stop("No samples with valid annotations to plot.") }

  seqs.df  = do.call(rbind, seqs.list)
  genes.df = do.call(rbind, genes.list)

  ##############################################################################
  # Step 3: Build links from shared gene names between adjacent samples
  ##############################################################################
  links.df = data.frame()

  if (show.links && length(genes.list) >= 2) {
    links.list = list()

    for (i in seq_len(length(genes.list) - 1)) {
      g1 = genes.list[[i]]
      g2 = genes.list[[i + 1]]

      shared = intersect(unique(g1$name), unique(g2$name))
      if (length(shared) == 0) { next }

      for (nm in shared) {
        r1 = g1[g1$name == nm, ][1, ]
        r2 = g2[g2$name == nm, ][1, ]
        links.list[[length(links.list) + 1]] = data.frame(
          seq_id  = seqs.df$seq_id[i],
          start   = r1$start,
          end     = r1$end,
          seq_id2 = seqs.df$seq_id[i + 1],
          start2  = r2$start,
          end2    = r2$end,
          name    = nm,
          type    = r1$type,
          stringsAsFactors = FALSE
        )
      }
    }

    if (length(links.list) > 0) {
      links.df = do.call(rbind, links.list)
    }
  }

  ##############################################################################
  # Step 4: Colour palette
  ##############################################################################
  type.colors = c(CDS      = "#4472C4",
                  tRNA     = "#ED7D31",
                  rRNA     = "#70AD47",
                  `D-loop` = "#C00000")

  ##############################################################################
  # Step 5: Build gggenomes plot
  ##############################################################################
  if (is.null(height)) {
    height = max(5, length(seqs.list) * 1.5 + 2)
  }

  have.links = show.links && nrow(links.df) > 0

  if (have.links) {
    p = gggenomes::gggenomes(seqs = seqs.df, genes = genes.df, links = links.df)
  } else {
    p = gggenomes::gggenomes(seqs = seqs.df, genes = genes.df)
  }

  # Sync sequences to minimise link crossings before adding layers
  if (sync.genomes && have.links) {
    p = gggenomes::sync(p)
  }

  # Sequence backbone lines and sample-name labels on the left
  p = p +
    gggenomes::geom_seq() +
    gggenomes::geom_seq_label()

  # Ribbon links (drawn behind gene arrows)
  if (have.links) {
    if (color.by == "type") {
      p = p + gggenomes::geom_link(ggplot2::aes(fill = type, colour = type),
                                   alpha = 0.35)
    } else {
      p = p + gggenomes::geom_link(ggplot2::aes(fill = name, colour = name),
                                   alpha = 0.35)
    }
  }

  # Gene arrows
  if (color.by == "type") {
    p = p + gggenomes::geom_gene(ggplot2::aes(fill = type), size = 4) +
      ggplot2::scale_fill_manual(name = "Feature type",
                                 values = type.colors, guide = "legend") +
      ggplot2::scale_colour_manual(name = "Feature type",
                                   values = type.colors, guide = "none")
  } else {
    p = p + gggenomes::geom_gene(ggplot2::aes(fill = name), size = 4) +
      ggplot2::guides(fill  = ggplot2::guide_legend(ncol = 6, title = "Gene"),
                      colour = "none")
  }

  # Optional gene name labels
  if (label.genes) {
    p = p + gggenomes::geom_gene_label(ggplot2::aes(label = name),
                                       size = 2, check_overlap = TRUE)
  }

  p = p +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position  = "bottom",
      legend.key.size  = ggplot2::unit(0.4, "cm"),
      legend.text      = ggplot2::element_text(size = 8),
      axis.text.y      = ggplot2::element_blank(),
      axis.ticks.y     = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank()
    )

  ##############################################################################
  # Step 6: Save
  ##############################################################################
  ggplot2::ggsave(out.path, plot = p,
                  width = width, height = height, units = "in", dpi = 300)
  message("Plot saved: ", out.path)

  return(invisible(p))

}#end function
