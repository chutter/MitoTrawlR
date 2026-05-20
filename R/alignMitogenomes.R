#' @title alignMitogenomes
#'
#' @description Concatenates per-marker phylip alignments from a folder into a
#'   single whole-mitogenome alignment. Each locus alignment is re-aligned
#'   against the reference marker using MAFFT (add mode), gaps relative to the
#'   reference are converted to Ns, missing samples are filled with Ns, and the
#'   resulting matrices are concatenated in reference order. Summary statistics
#'   (base-pair counts and proportions per locus per sample) are written to
#'   logs/.
#'
#' @param alignment.folder path to a folder of per-locus phylip alignments
#'   (untrimmed).
#'
#' @param draft.contigs path to the folder of draft contig FASTA files (one
#'   .fa file per sample); used to determine the set of sample names.
#'
#' @param reference.name name of the reference directory created by
#'   \code{buildReference}, which must contain \code{refMarkers.fa}.
#'
#' @param output.dir path to the output directory where the concatenated
#'   mitogenome alignment and log files will be written.
#'
#' @param dataset.name prefix string used for all output file names.
#'
#' @param mafft.path path to the directory containing \code{mafft}. NULL uses
#'   the system PATH.
#'
#' @param threads number of CPU threads to pass to MAFFT.
#'
#' @param overwrite logical; if TRUE, the output directory is deleted and
#'   recreated before writing results.
#'
#' @return Invisibly returns NULL. Writes a concatenated mitogenome phylip
#'   alignment, a feature coordinate table, and per-locus bp-count and
#'   bp-proportion CSV files to \code{output.dir}.
#'
#' @examples
#' \dontrun{
#' alignMitogenomes(
#'   alignment.folder = "Alignments/untrimmed-alignments",
#'   draft.contigs    = "draftContigs",
#'   reference.name   = "reference",
#'   output.dir       = "MitoGenomes",
#'   dataset.name     = "untrimmed",
#'   overwrite        = FALSE
#' )
#' }
#'
#' @export

alignMitogenomes = function(alignment.folder = NULL,
                            draft.contigs = "draftContigs",
                            reference.name = "reference",
                            output.dir = "MitoGenomes",
                            dataset.name = "untrimmed",
                            mafft.path = NULL,
                            threads = 1,
                            overwrite = FALSE) {

  # #Debug
  #  alignment.folder = "Alignments/untrimmed-alignments"
  #  reference.name = "reference"
  #  draft.contigs = "draftContigs"
  #  output.dir = "MitoGenomes"
  #  dataset.name = "untrimmed"
  #  overwrite = TRUE

  if (is.null(mafft.path) == FALSE){
    b.string = unlist(strsplit(mafft.path, ""))
    if (b.string[length(b.string)] != "/") {
      mafft.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { mafft.path = "" }

  if (dir.exists(output.dir) == FALSE) {
    dir.create(output.dir)
  } else if (overwrite == TRUE) {
    unlink(output.dir, recursive = TRUE)
    dir.create(output.dir)
  } else {
    stop("overwrite is FALSE and output directory '", output.dir, "' already exists.")
  }

  if (dir.exists(paste0(output.dir, "/logs")) == FALSE) {
    dir.create(paste0(output.dir, "/logs"))
  } else if (overwrite == TRUE) {
    unlink(paste0(output.dir, "/logs"), recursive = TRUE)
    dir.create(paste0(output.dir, "/logs"))
  }

  #Gets the samples
  locus.names = list.files(alignment.folder)
  locus.names = gsub(".phy$", "", locus.names)
  sample.names = list.files(draft.contigs)
  sample.names = gsub(".fa$", "", sample.names)

  #Create new directories
  if (dir.exists(paste0(output.dir, "/alignments")) == FALSE) { dir.create(paste0(output.dir, "/alignments")) }
  if (dir.exists(paste0(output.dir, "/sample-genomes")) == FALSE) { dir.create(paste0(output.dir, "/sample-genomes")) }

  #Sets up the loci to align
  ref.data = Biostrings::readDNAStringSet(paste0(reference.name, "/refMarkers.fa"))

  #Header data for features and whatnot
  header.data = c("Sample", locus.names)

  collect.data.bp = data.table::data.table(matrix(as.numeric(0), nrow = length(sample.names), ncol = length(header.data)))
  data.table::setnames(collect.data.bp, header.data)
  collect.data.bp[, Sample := as.character(sample.names)]

  collect.data.pr = data.table::data.table(matrix(as.numeric(0), nrow = length(sample.names), ncol = length(header.data)))
  data.table::setnames(collect.data.pr, header.data)
  collect.data.pr[, Sample := as.character(sample.names)]

  feat.headers = c("Marker", "Start", "End")
  feature.data = data.table::data.table(matrix(as.numeric(0), nrow = length(locus.names), ncol = length(feat.headers)))
  data.table::setnames(feature.data, feat.headers)
  feature.data[, Marker := as.character(Marker)]

  draft.genome = c()
  for (i in seq_along(locus.names)){

    align = Biostrings::DNAStringSet(Biostrings::readAAMultipleAlignment(
      file = paste0(alignment.folder, "/", locus.names[i], ".phy"), format = "phylip"))

    #Gathers reference
    ref.seq = ref.data[names(ref.data) %in% locus.names[i]]
    names(ref.seq) = "Reference"

    new.align = PhyloProcessR::runMafft(sequence.data = align,
                                        add.contigs = ref.seq,
                                        algorithm = "add",
                                        adjust.direction = TRUE,
                                        threads = threads,
                                        cleanup.files = TRUE,
                                        mafft.path = mafft.path,
                                        quiet = TRUE)

    mat.align = strsplit(as.character(new.align), "")
    mat.align = lapply(mat.align, tolower)

    # Fill gaps aligned to reference positions with Ns
    save.matrix = c()
    for (j in seq_along(mat.align)){

      if (names(mat.align)[j] == "Reference"){ next }

      tar.seq = mat.align[[j]]
      ref.row = mat.align[[grep("Reference", names(mat.align))]]
      for (k in seq_along(ref.row)){
        if (ref.row[k] != "-" && tar.seq[k] == "-"){ tar.seq[k] = "n" }
      }#end k loop

      save.seq = matrix(tar.seq, nrow = 1)
      rownames(save.seq) = names(mat.align)[j]
      save.matrix = rbind(save.matrix, save.seq)

    }#end j loop

    #Fill in missing samples
    miss.samples = sample.names[!sample.names %in% rownames(save.matrix)]
    if (length(miss.samples) != 0){
      for (j in seq_along(miss.samples)){
        save.seq = matrix(as.character("n"), nrow = 1, ncol = ncol(save.matrix))
        rownames(save.seq) = miss.samples[j]
        save.matrix = rbind(save.matrix, save.seq)
      }#end j loop
    }#end if

    #Sort and concatenate alignments
    order.matrix = save.matrix[order(rownames(save.matrix)),]

    if (is.null(ncol(draft.genome))){
      start = 1
      end = ncol(order.matrix)
    } else {
      start = ncol(draft.genome) + 1
      end = start + ncol(order.matrix) - 1
    }#end else

    draft.genome = cbind(draft.genome, order.matrix)

    #Saves location data
    data.table::set(feature.data, i = as.integer(i), j = match("Marker", feat.headers), value = locus.names[i])
    data.table::set(feature.data, i = as.integer(i), j = match("Start", feat.headers), value = start)
    data.table::set(feature.data, i = as.integer(i), j = match("End", feat.headers), value = end)

    #Collects stats data
    bp.data = apply(order.matrix, MARGIN = 1, FUN = function(x) length(x[x != "n"]))

    data.table::set(collect.data.bp, i = match(names(bp.data), collect.data.bp$Sample),
                    j = match(locus.names[i], header.data), value = bp.data)
    data.table::set(collect.data.pr, i = match(names(bp.data), collect.data.pr$Sample),
                    j = match(locus.names[i], header.data), value = round(bp.data / max(bp.data), 3))

  }#end i loop

  #Save files
  write.csv(collect.data.bp, file = paste0(output.dir, "/logs/", dataset.name, "_mito-alignment_bp-count.csv"), row.names = FALSE)
  write.csv(collect.data.pr, file = paste0(output.dir, "/logs/", dataset.name, "_mito-alignment_bp-prop.csv"), row.names = FALSE)

  write.genome = as.matrix(ape::as.DNAbin(draft.genome))
  PhyloProcessR::writePhylip(write.genome, file = paste0(output.dir, "/alignments/", dataset.name, "_mitogenome_alignment.phy"), interleave = FALSE)
  write.table(feature.data, file = paste0(output.dir, "/alignments/", dataset.name, "_alignment_feature_table.txt"), row.names = FALSE, quote = FALSE)

  return(invisible(NULL))

}#end function
