#' @title buildPhylogeny
#'
#' @description Infers a maximum-likelihood phylogeny from a concatenated
#'   mitogenome alignment using IQ-TREE 3. A NEXUS partition file is
#'   automatically built from the feature coordinate table produced by
#'   \code{alignMitogenomes}: protein-coding genes (CDS) are split into codon
#'   positions 1+2 (grouped) and 3 (\code{partition.scheme = "codon"}), or each
#'   marker receives a single partition (\code{"byMarker"}), or no partition
#'   file is used and model selection runs on the whole alignment
#'   (\code{"auto"}). All partitioned schemes run \code{MFP+MERGE} so IQ-TREE
#'   finds the best model per partition and merges statistically
#'   indistinguishable ones.
#'
#' @param alignment.file path to the concatenated mitogenome alignment (phylip
#'   format) produced by \code{alignMitogenomes}.
#'
#' @param feature.table path to the alignment feature coordinate table
#'   (\code{*_alignment_feature_table.txt}) produced by
#'   \code{alignMitogenomes}. Required when \code{partition.scheme} is
#'   \code{"codon"} or \code{"byMarker"}.
#'
#' @param output.dir path to the directory where IQ-TREE output files will be
#'   written.
#'
#' @param dataset.name prefix used for all output file names.
#'
#' @param iqtree.path path to the directory containing \code{iqtree3}. NULL
#'   uses the system PATH.
#'
#' @param threads number of CPU threads to pass to IQ-TREE (\code{-T}).
#'
#' @param bootstrap number of ultrafast bootstrap replicates (\code{-B});
#'   set to 0 or NULL to skip bootstrapping.
#'
#' @param partition.scheme partitioning strategy. \code{"codon"} (default)
#'   splits each CDS marker into codon positions 1+2 and 3, with a single
#'   partition for each tRNA, rRNA, and D-loop region, then runs
#'   \code{MFP+MERGE}. \code{"byMarker"} assigns one partition per marker and
#'   runs \code{MFP+MERGE}. \code{"auto"} uses no partition file and runs
#'   \code{MFP} on the whole alignment.
#'
#' @param overwrite logical; if TRUE, existing output files are overwritten
#'   (passes \code{--redo} to IQ-TREE and recreates \code{output.dir}).
#'
#' @return Invisibly returns NULL. Writes IQ-TREE output files (tree, log,
#'   model report, partition file) to \code{output.dir}.
#'
#' @examples
#' \dontrun{
#' buildPhylogeny(
#'   alignment.file   = "MitoGenomes/alignments/untrimmed_mitogenome_alignment.phy",
#'   feature.table    = "MitoGenomes/alignments/untrimmed_alignment_feature_table.txt",
#'   output.dir       = "Phylogeny",
#'   dataset.name     = "mitogenome",
#'   threads          = 4,
#'   bootstrap        = 1000,
#'   partition.scheme = "codon",
#'   overwrite        = FALSE
#' )
#' }
#'
#' @export

buildPhylogeny = function(alignment.file = NULL,
                          feature.table = NULL,
                          output.dir = NULL,
                          dataset.name = "mitogenome",
                          iqtree.path = NULL,
                          threads = 1,
                          bootstrap = 1000,
                          partition.scheme = c("codon", "byMarker", "auto"),
                          overwrite = FALSE) {

  # #Debug
  # alignment.file   = "MitoGenomes/alignments/untrimmed_mitogenome_alignment.phy"
  # feature.table    = "MitoGenomes/alignments/untrimmed_alignment_feature_table.txt"
  # output.dir       = "Phylogeny"
  # dataset.name     = "mitogenome"
  # threads          = 4
  # bootstrap        = 1000
  # partition.scheme = "codon"
  # overwrite        = FALSE

  partition.scheme = match.arg(partition.scheme)

  if (is.null(alignment.file)){ stop("Please provide an alignment file.") }
  if (!file.exists(alignment.file)){ stop("Alignment file not found: ", alignment.file) }

  if (partition.scheme != "auto" && is.null(feature.table)){
    stop("feature.table is required when partition.scheme is '", partition.scheme, "'.")
  }
  if (!is.null(feature.table) && !file.exists(feature.table)){
    stop("Feature table not found: ", feature.table)
  }

  if (is.null(iqtree.path) == FALSE){
    b.string = unlist(strsplit(iqtree.path, ""))
    if (b.string[length(b.string)] != "/") {
      iqtree.path = paste0(append(b.string, "/"), collapse = "")
    }
  } else { iqtree.path = "" }

  if (dir.exists(output.dir) == FALSE) {
    dir.create(output.dir, recursive = TRUE)
  } else if (overwrite == TRUE) {
    unlink(output.dir, recursive = TRUE)
    dir.create(output.dir, recursive = TRUE)
  } else {
    stop("overwrite is FALSE and output directory '", output.dir, "' already exists.")
  }

  ##############################################################################
  # Build partition file
  ##############################################################################
  partition.file = NULL

  if (partition.scheme != "auto"){

    feat = read.table(feature.table, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

    charset.lines  = character(0)
    partition.names = character(0)

    for (i in seq_len(nrow(feat))){
      marker = feat$Marker[i]
      start  = feat$Start[i]
      end    = feat$End[i]

      # Parse functional type from NNN_type_name format
      type = gsub("^[0-9]+_([^_]+)_.*", "\\1", marker)

      # Clean label: strip numeric prefix and type, replace non-alphanumeric
      label = gsub("^[0-9]+_[^_]+_", "", marker)
      label = gsub("[^A-Za-z0-9]", "_", label)

      if (partition.scheme == "codon" && type == "CDS"){
        # Codon positions 1+2 grouped; position 3 separate
        pos12.label = paste0(label, "_pos12")
        pos3.label  = paste0(label, "_pos3")
        charset.lines = c(charset.lines,
          paste0("  charset ", pos12.label, " = ", start, "-", end, "\\3 ", (start + 1), "-", end, "\\3;"),
          paste0("  charset ", pos3.label,  " = ", (start + 2), "-", end, "\\3;"))
        partition.names = c(partition.names, pos12.label, pos3.label)
      } else {
        # Single partition for tRNA, rRNA, D_loop, or byMarker scheme
        charset.lines = c(charset.lines,
          paste0("  charset ", label, " = ", start, "-", end, ";"))
        partition.names = c(partition.names, label)
      }
    }#end for

    nexus.lines = c(
      "#nexus",
      "begin sets;",
      charset.lines,
      paste0("  partition MyPartition = ", length(partition.names), ": ",
             paste(partition.names, collapse = ", "), ";"),
      "  set partition = MyPartition;",
      "end;"
    )

    partition.file = paste0(output.dir, "/", dataset.name, "_partitions.nex")
    writeLines(nexus.lines, con = partition.file)
    message("Partition file written: ", partition.file,
            " (", length(partition.names), " partitions)")

  }#end partition build

  ##############################################################################
  # Build IQ-TREE command
  ##############################################################################
  prefix = paste0(output.dir, "/", dataset.name)

  if (partition.scheme == "auto"){
    model.flag = "-m MFP"
    input.flag = paste0("-s ", alignment.file)
  } else {
    model.flag = "-m MFP+MERGE"
    input.flag = paste0("-s ", alignment.file, " -p ", partition.file)
  }

  boot.flag = if (!is.null(bootstrap) && bootstrap > 0) paste0("-B ", bootstrap) else ""
  redo.flag = if (overwrite) "--redo" else ""

  iqtree.cmd = paste(
    paste0(iqtree.path, "iqtree3"),
    input.flag,
    model.flag,
    boot.flag,
    paste0("-T ", threads),
    paste0("--prefix ", prefix),
    redo.flag
  )

  # Strip extra spaces from optional empty flags
  iqtree.cmd = gsub(" +", " ", trimws(iqtree.cmd))

  message("Running IQ-TREE 3...")
  message(iqtree.cmd)
  system(iqtree.cmd)

  # Check that a tree was produced
  tree.file = paste0(prefix, ".treefile")
  if (!file.exists(tree.file)){
    warning("IQ-TREE did not produce a tree file. Check ", prefix, ".log for errors.")
  } else {
    message("Phylogeny complete. Tree: ", tree.file)
  }

  return(invisible(NULL))

}#end function
