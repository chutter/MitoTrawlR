#' @title plotMitoGenomes
#'
#' @description Plots a maximum-likelihood phylogeny alongside a per-marker
#'   mitogenome completeness heatmap. Each cell in the heatmap shows the
#'   proportion of base pairs recovered for that sample and marker relative to
#'   the best-covered sample. An optional coloured strip above the heatmap
#'   indicates the functional type of each marker (CDS, tRNA, rRNA, D-loop).
#'   The plot is saved as a PDF or PNG file and returned invisibly.
#'
#' @param tree.file path to a tree file in Newick format (e.g., the
#'   \code{.treefile} or \code{.contree} produced by \code{buildPhylogeny}).
#'
#' @param bp.table path to the per-locus base-pair count CSV produced by
#'   \code{alignMitogenomes} (e.g.,
#'   \code{MitoGenomes/logs/untrimmed_mito-alignment_bp-count.csv}).
#'
#' @param feature.table path to the alignment feature coordinate table
#'   (\code{*_alignment_feature_table.txt}) produced by
#'   \code{alignMitogenomes}. Required when \code{show.marker.types = TRUE}.
#'
#' @param output.file path (without extension) for the saved plot.
#'
#' @param output.format output file format: \code{"pdf"} or \code{"png"}.
#'
#' @param width plot width in inches.
#'
#' @param height plot height in inches.
#'
#' @param tip.label.size font size for tip labels on the tree.
#'
#' @param col.label.size font size for marker name column labels on the
#'   heatmap.
#'
#' @param tree.width relative width of the tree panel versus the heatmap
#'   panel (0--1; e.g. 0.35 gives the tree 35\% of the total width).
#'
#' @param outgroup character vector of outgroup tip label(s) used to root the
#'   tree. NULL leaves the tree unrooted.
#'
#' @param show.marker.types logical; if TRUE and \code{feature.table} is
#'   provided, a narrow coloured strip above the heatmap shows the functional
#'   type of each marker (CDS = blue, tRNA = orange, rRNA = green,
#'   D-loop = red).
#'
#' @param color.complete fill colour for 100\% completeness cells.
#'
#' @param color.missing fill colour for 0\% completeness cells.
#'
#' @param overwrite logical; if TRUE, an existing output file is overwritten.
#'
#' @return Invisibly returns the combined plot object. Writes a PDF or PNG to
#'   \code{output.file}.
#'
#' @examples
#' \dontrun{
#' plotMitoGenomes(
#'   tree.file        = "Phylogeny/mitogenome.treefile",
#'   bp.table         = "MitoGenomes/logs/untrimmed_mito-alignment_bp-count.csv",
#'   feature.table    = "MitoGenomes/alignments/untrimmed_alignment_feature_table.txt",
#'   output.file      = "mitogenome_completeness",
#'   output.format    = "pdf",
#'   width            = 16,
#'   height           = 10,
#'   outgroup         = "OutgroupSpecies",
#'   show.marker.types = TRUE,
#'   overwrite        = TRUE
#' )
#' }
#'
#' @export

plotMitoGenomes = function(tree.file = NULL,
                           bp.table = NULL,
                           feature.table = NULL,
                           genome.dir = "MitoGenomes",
                           output.file = NULL,
                           output.format = c("pdf", "png"),
                           width = 16,
                           height = 10,
                           tip.label.size = 2,
                           col.label.size = 6,
                           tree.width = 0.50,
                           outgroup = NULL,
                           show.marker.types = TRUE,
                           color.complete = "steelblue4",
                           color.missing = "white",
                           overwrite = TRUE) {

  # #Debug
  # tree.file     = "Phylogeny/mitogenome.treefile"
  # bp.table      = "MitoGenomes/logs/untrimmed_mito-alignment_bp-count.csv"
  # feature.table = "MitoGenomes/alignments/untrimmed_alignment_feature_table.txt"
  # output.file   = "mitogenome_completeness"
  # outgroup      = NULL
  # show.marker.types = TRUE

  for (pkg in c("ggplot2", "ggtree", "patchwork", "scales")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("Package '", pkg, "' is required for plotMitoGenomes(). Install it with: install.packages('", pkg, "')")
  }

  output.format = match.arg(output.format)

  if (is.null(output.file)) {
    figures.dir = paste0(genome.dir, "/figures")
    dir.create(figures.dir, recursive = TRUE, showWarnings = FALSE)
    output.file = paste0(figures.dir, "/mitogenome_completeness")
  }

  if (is.null(tree.file) || !file.exists(tree.file)){
    stop("Tree file not found: ", tree.file)
  }
  if (is.null(bp.table) || !file.exists(bp.table)){
    stop("bp.table file not found: ", bp.table)
  }

  out.path = paste0(output.file, ".", output.format)
  if (file.exists(out.path) && overwrite == FALSE){
    stop("Output file already exists and overwrite = FALSE: ", out.path)
  }

  ##############################################################################
  # Step 1: Read and prepare tree
  ##############################################################################
  tree = ape::read.tree(tree.file)

  if (!is.null(outgroup)){
    valid.out = outgroup[outgroup %in% tree$tip.label]
    if (length(valid.out) == 0){
      warning("None of the supplied outgroup labels found in tree. Skipping rooting.")
    } else {
      tree = ape::root(tree, outgroup = valid.out, resolve.root = TRUE)
    }
  }

  tree = ape::ladderize(tree)

  ##############################################################################
  # Step 2: Build completeness matrix from bp-count table
  ##############################################################################
  bp.data = utils::read.csv(bp.table, stringsAsFactors = FALSE, check.names = FALSE)
  rownames(bp.data) = bp.data$Sample
  bp.data$Sample = NULL
  bp.mat = as.matrix(bp.data)

  # Per-locus normalisation: proportion of maximum coverage seen across samples
  col.max = apply(bp.mat, 2, max)
  col.max[col.max == 0] = 1
  comp.mat = sweep(bp.mat, 2, col.max, "/")

  # Fill any tree tips absent from the bp table with 0
  missing.tips = setdiff(tree$tip.label, rownames(comp.mat))
  if (length(missing.tips) > 0){
    fill.mat = matrix(0, nrow = length(missing.tips), ncol = ncol(comp.mat),
                      dimnames = list(missing.tips, colnames(comp.mat)))
    comp.mat = rbind(comp.mat, fill.mat)
    message(length(missing.tips), " tip(s) not found in bp.table — filled with 0.")
  }

  # Strip numeric prefix for display: "003_CDS_ND1" -> "ND1"
  # make.unique handles duplicate short names (e.g. two tRNA-Leu)
  col.short = make.unique(gsub("^[0-9]+_[^_]+_", "", colnames(comp.mat)))
  colnames(comp.mat) = col.short

  ##############################################################################
  # Step 3: Build tree panel (pure ggplot2 – ggtree used for layout data only)
  # Rendering with ggtree layers triggers an incompatibility with ggplot2 >= 3.5
  # where ggtree's old S3 guide objects crash the new R6 guide-build system.
  ##############################################################################
  p.tree.gt = ggtree::ggtree(tree, ladderize = FALSE)   # object only, never rendered
  td        = p.tree.gt$data
  tip.data  = td[td$isTip, ]
  n.tips    = nrow(tip.data)

  # Numeric y per tip (y=1 at bottom, y=n at top) – used to align heatmap rows
  tip.y.map = stats::setNames(tip.data$y, tip.data$label)
  # Tips in bottom-to-top order for subsetting comp.mat
  plot.tips = tip.data$label[order(tip.data$y)]
  plot.tips = plot.tips[plot.tips %in% rownames(comp.mat)]
  comp.mat  = comp.mat[plot.tips, , drop = FALSE]

  # Horizontal branch segments
  edge.df = td[!is.na(td$parent) & td$node != td$parent & !is.na(td$branch), ]
  h.segs  = data.frame(x    = edge.df$branch, xend = edge.df$x,
                        y    = edge.df$y,      yend = edge.df$y)

  # Vertical segments at each internal node connecting its children's y range
  int.nodes = unique(td$parent[!is.na(td$parent) & td$node != td$parent])
  v.list    = lapply(int.nodes, function(p) {
    nx = td$x[td$node == p]
    ch = td[!is.na(td$parent) & td$parent == p & td$node != p, ]
    if (length(nx) == 0 || nrow(ch) < 2) return(NULL)
    data.frame(x = nx[1], xend = nx[1], y = min(ch$y), yend = max(ch$y))
  })
  v.segs   = do.call(rbind, Filter(Negate(is.null), v.list))
  all.segs = rbind(h.segs, v.segs)

  p.tree = ggplot2::ggplot() +
    ggplot2::geom_segment(data = all.segs,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      colour = "black", linewidth = 0.4) +
    ggplot2::geom_text(data = tip.data,
      ggplot2::aes(x = x, y = y, label = label),
      hjust = -0.1, size = tip.label.size) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.5))) +
    ggplot2::scale_y_continuous(limits = c(0.5, n.tips + 0.5), expand = c(0, 0)) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(5, 0, 5, 5))

  ##############################################################################
  # Step 4: Build completeness heatmap panel
  ##############################################################################
  comp.long = data.frame(
    Sample       = rep(rownames(comp.mat), times = ncol(comp.mat)),
    Marker       = rep(colnames(comp.mat), each = nrow(comp.mat)),
    Completeness = as.vector(comp.mat),
    stringsAsFactors = FALSE
  )
  # Use same numeric y as the tree for pixel-perfect row alignment
  comp.long$y_num = tip.y.map[comp.long$Sample]
  comp.long$Marker = factor(comp.long$Marker, levels = colnames(comp.mat))

  p.heat = ggplot2::ggplot(comp.long,
                           ggplot2::aes(x = Marker, y = y_num, fill = Completeness)) +
    ggplot2::geom_tile(color = "grey85", linewidth = 0.1) +
    ggplot2::scale_fill_gradient(name = "Completeness",
                                 low  = color.missing,
                                 high = color.complete,
                                 limits = c(0, 1),
                                 labels = scales::percent_format(accuracy = 1)) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::scale_y_continuous(limits = c(0.5, n.tips + 0.5), expand = c(0, 0),
                                 breaks = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_text(angle = -55, hjust = 1, vjust = 0.5,
                                            size = col.label.size),
      axis.text.y  = ggplot2::element_blank(),
      axis.ticks   = ggplot2::element_blank(),
      axis.title   = ggplot2::element_blank(),
      panel.grid   = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.key.width = ggplot2::unit(1.5, "cm"),
      plot.margin  = ggplot2::margin(5, 5, 5, 0)
    )

  ##############################################################################
  # Step 5: Optional marker-type strip above the heatmap
  ##############################################################################
  if (show.marker.types && !is.null(feature.table) && file.exists(feature.table)){

    feat = utils::read.table(feature.table, header = TRUE,
                             stringsAsFactors = FALSE)
    feat$label = gsub("^[0-9]+_[^_]+_", "", feat$Marker)
    feat$type  = gsub("^[0-9]+_([^_]+)_.*", "\\1", feat$Marker)
    feat = feat[feat$label %in% col.short, ]
    feat$label = factor(feat$label, levels = col.short)

    type.colors = c(CDS    = "#4472C4",
                    tRNA   = "#ED7D31",
                    rRNA   = "#70AD47",
                    D_loop = "#C00000")
    unknown.types = setdiff(unique(feat$type), names(type.colors))
    if (length(unknown.types) > 0){
      type.colors = c(type.colors, setNames(rep("grey60", length(unknown.types)),
                                            unknown.types))
    }

    p.type = ggplot2::ggplot(feat,
                             ggplot2::aes(x = label, y = 1, fill = type)) +
      ggplot2::geom_tile(color = "grey85", linewidth = 0.1) +
      ggplot2::scale_fill_manual(name = "Marker Type", values = type.colors) +
      ggplot2::scale_x_discrete(limits = col.short) +
      ggplot2::theme_void() +
      ggplot2::theme(
        legend.position  = "bottom",
        legend.key.size  = ggplot2::unit(0.4, "cm"),
        legend.text      = ggplot2::element_text(size = 8),
        plot.margin      = ggplot2::margin(0, 5, 0, 0)
      )

    # Stack heatmap above type strip, place tree to the left
    heat.col   = p.heat / p.type +
      patchwork::plot_layout(heights = c(0.95, 0.05))
    final.plot = (p.tree | heat.col) +
      patchwork::plot_layout(widths = c(tree.width, 1 - tree.width),
                             guides = "collect") &
      ggplot2::theme(legend.position = "bottom")

  } else {

    final.plot = (p.tree | p.heat) +
      patchwork::plot_layout(widths = c(tree.width, 1 - tree.width),
                             guides = "collect") &
      ggplot2::theme(legend.position = "bottom")

  }

  ##############################################################################
  # Step 6: Save
  ##############################################################################
  ggplot2::ggsave(out.path, plot = final.plot,
                  width = width, height = height, units = "in", dpi = 300)
  message("Plot saved: ", out.path)

  return(invisible(final.plot))

}#end function
