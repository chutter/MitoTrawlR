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

  # Guard on the treefile, not the directory, so the function can be re-run
  # with a different dataset.name into the same output directory.
  out.tree = file.path(output.dir, paste0(dataset.name, ".treefile"))
  if (file.exists(out.tree) && overwrite == FALSE) {
    message("Tree file '", out.tree, "' already exists and overwrite = FALSE. Skipping.")
    return(invisible(NULL))
  }
  dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

  ##############################################################################
  # Build partition file
  ##############################################################################
  partition.file = NULL

  if (partition.scheme != "auto"){

    feat = read.table(feature.table, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

    part.lines = character(0)

    for (i in seq_len(nrow(feat))){
      marker = feat$Marker[i]
      start  = as.integer(feat$Start[i])
      end    = as.integer(feat$End[i])

      # Parse functional type from NNN_type_name format
      type = gsub("^[0-9]+_([^_]+)_.*", "\\1", marker)

      # Clean label: strip numeric prefix and type, replace non-alphanumeric
      label = gsub("^[0-9]+_[^_]+_", "", marker)
      label = gsub("[^A-Za-z0-9]", "_", label)

      if (partition.scheme == "codon" && type == "CDS"){
        # Three partitions per CDS gene (one per codon position)
        # Listed as individual sites to avoid stride-notation compatibility issues
        len = end - start + 1
        s1 = seq(start,     end, by = 3)
        s2 = seq(start + 1, end, by = 3)
        s3 = seq(start + 2, end, by = 3)
        part.lines = c(part.lines,
          paste0("DNA, ", label, "_pos1 = ", paste(s1, collapse = ", ")),
          paste0("DNA, ", label, "_pos2 = ", paste(s2, collapse = ", ")),
          paste0("DNA, ", label, "_pos3 = ", paste(s3, collapse = ", ")))
      } else {
        # Single partition per marker (tRNA, rRNA, D-loop, or byMarker scheme)
        part.lines = c(part.lines,
          paste0("DNA, ", label, " = ", start, "-", end))
      }
    }#end for

    partition.file = paste0(output.dir, "/", dataset.name, "_partitions.txt")
    writeLines(part.lines, con = partition.file)
    message("Partition file written: ", partition.file,
            " (", length(part.lines), " partitions)")

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
    "-st DNA",
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
