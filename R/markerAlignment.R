#' @title markerAlignment
#'
#' @description Aligns per-sample marker sequences across all samples for each
#'   reference marker using MAFFT (localpair, high accuracy). Sequences that
#'   are too divergent from the reference (pairwise distance > max.distance)
#'   are removed and the alignment is rerun. Leading and trailing gaps relative
#'   to the reference are converted to Ns before writing the final trimmed
#'   alignment as a phylip file. Alignments with fewer than three taxa are
#'   skipped.
#'
#' @param input.folder path to the folder of per-sample marker FASTA files as
#'   produced by \code{annotateMitoContigs} (typically
#'   \code{Annotations/sample-markers}).
#'
#' @param reference.name name of the reference directory created by
#'   \code{buildReference}, which must contain \code{refMarkers.fa}.
#'
#' @param threads number of CPU threads to pass to MAFFT.
#'
#' @param mafft.path path to the directory containing \code{mafft}. NULL uses
#'   the system PATH.
#'
#' @param max.distance maximum pairwise distance from the reference sequence
#'   above which a sample sequence is excluded before realignment (0--1).
#'
#' @param overwrite logical; if TRUE, existing \code{Alignments/} directories
#'   are deleted and recreated.
#'
#' @return Invisibly returns NULL. Writes per-locus phylip alignments to
#'   \code{Alignments/untrimmed-alignments/} and unaligned FASTA files to
#'   \code{Alignments/unaligned-markers/}.
#'
#' @examples
#' \dontrun{
#' markerAlignment(
#'   input.folder   = "Annotations/sample-markers",
#'   reference.name = "reference",
#'   threads        = 4,
#'   mafft.path     = "/path/to/mafft/bin",
#'   max.distance   = 0.40,
#'   overwrite      = FALSE
#' )
#' }
#'
#' @export

#Aligns all the different markers
markerAlignment = function(input.folder = NULL,
                           reference.name = NULL,
                           output.dir = "Alignments",
                           threads = 1,
                           mafft.path = NULL,
                           max.distance = 0.40,
                           overwrite = FALSE){

  #Debug
  #input.folder = "Annotations/sample-markers"
  #reference.name = "Reference"
  #threads = 4
  #overwrite = TRUE

  if (is.null(input.folder)) { stop("Please provide an input folder path.") }
  if (is.null(reference.name)) { stop("Please provide a reference name.") }

  if (is.null(mafft.path) == FALSE){
    b.string = unlist(strsplit(mafft.path, ""))
    if (b.string[length(b.string)] != "/") {
      mafft.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { mafft.path = "" }

  #Checks for overwrite
  if (dir.exists(output.dir) == FALSE) {
    dir.create(output.dir)
  } else if (overwrite == TRUE) {
    unlink(output.dir, recursive = TRUE)
    dir.create(output.dir)
  } else { message("Alignments directory already exists and overwrite = FALSE. Skipping."); return(invisible(NULL)) }

  if (dir.exists(paste0(output.dir, "/untrimmed-alignments")) == FALSE) {
    dir.create(paste0(output.dir, "/untrimmed-alignments"))
  } else if (overwrite == TRUE) {
    unlink(paste0(output.dir, "/untrimmed-alignments"), recursive = TRUE)
    dir.create(paste0(output.dir, "/untrimmed-alignments"))
  }#end dir exists

  if (dir.exists(paste0(output.dir, "/unaligned-markers")) == FALSE) {
    dir.create(paste0(output.dir, "/unaligned-markers"))
  } else if (overwrite == TRUE) {
    unlink(paste0(output.dir, "/unaligned-markers"), recursive = TRUE)
    dir.create(paste0(output.dir, "/unaligned-markers"))
  }#end dir exists

  # Helper: run MAFFT on a named marker, read result, fix orientation, rename
  run.mafft = function(marker.name) {
    in.fa  = paste0(output.dir, "/unaligned-markers/", marker.name, ".fa")
    out.fa = paste0(output.dir, "/unaligned-markers/", marker.name, "_align.fa")
    system(paste0(mafft.path, "mafft --localpair --maxiterate 1000 --adjustdirection --quiet --op 3 --ep 0.123",
                  " --thread ", threads, " ", in.fa, " > ", out.fa))
    aln = Biostrings::readDNAStringSet(out.fa)
    reversed = names(aln)[grep("_R_", names(aln))]
    if (length(reversed[grep("Reference", reversed)]) == 1) {
      aln = Biostrings::reverseComplement(aln)
    }
    names(aln) = gsub("_R_", "", names(aln))
    names(aln) = gsub("_\\|_.*", "", names(aln))
    aln
  }

  #Sets up the loci to align
  ref.data = Biostrings::readDNAStringSet(paste0(reference.name, "/refMarkers.fa"))

  #Gets the samples
  spp.samples = list.files(input.folder)
  spp.samples = spp.samples[spp.samples != ""]

  #Gets all species data
  final.contigs = Biostrings::DNAStringSet()
  for (i in seq_along(spp.samples)){
    spp.data = Biostrings::readDNAStringSet(paste0(input.folder, "/", spp.samples[i]))
    final.contigs = append(final.contigs, spp.data)
  }#end i loop

  #Aligns each potential locus
  for (i in seq_along(ref.data)){

    ##############
    # STEP 1: Set up sequences for this marker
    ##############
    temp.name = gsub("^[0-9]+_[^_]+_", "", names(ref.data)[i])
    sample.markers = final.contigs[grep(paste0(temp.name, "$"), names(final.contigs))]

    if (length(sample.markers) <= 2){
      message(names(ref.data)[i], " had too few taxa.")
      next
    }

    #Adds reference locus
    align.data = append(sample.markers, ref.data[i])
    names(align.data)[length(align.data)] = "Reference"
    final.save = as.list(as.character(align.data))

    PhyloProcessR::writeFasta(sequences = final.save, names = names(final.save),
               paste0(output.dir, "/unaligned-markers/", names(ref.data)[i], ".fa"), nbchar = 1000000, as.string = TRUE)

    ##############
    # STEP 2: Run MAFFT and filter divergent sequences
    ##############
    alignment = run.mafft(names(ref.data)[i])

    dist.data = PhyloProcessR::pairwiseDistanceTarget(alignment = alignment, target = "Reference")
    good.seqs = which(dist.data <= max.distance)
    rem.align = alignment[as.numeric(good.seqs),]
    rem.align = rem.align[!duplicated(names(rem.align))]

    if (length(rem.align) <= 2){
      message(names(ref.data)[i], " had too few taxa after distance filtering.")
      next
    }

    ##############
    # STEP 3: Realign if divergent sequences were removed
    ##############
    if (length(good.seqs) != length(alignment)){
      final.save = as.list(as.character(rem.align))
      PhyloProcessR::writeFasta(sequences = final.save, names = names(final.save),
                           paste0(output.dir, "/unaligned-markers/", names(ref.data)[i], ".fa"), nbchar = 1000000, as.string = TRUE)
      alignment = run.mafft(names(ref.data)[i])
      rem.align = alignment[!duplicated(names(alignment))]
    }

    ##############
    # STEP 4: Trim to reference bounds and replace edge gaps with Ns
    ##############
    ref.aligned = as.character(rem.align['Reference'])
    not.gaps = stringr::str_locate_all(ref.aligned, pattern = "[^-]")[[1]][,1]
    ref.start = min(not.gaps)
    ref.finish = max(not.gaps)
    trim.align = Biostrings::subseq(rem.align, ref.start, ref.finish)
    trim.align = trim.align[names(trim.align) != "Reference"]

    sample.n = names(trim.align)
    align.list = vector(mode = "list", length = length(sample.n))
    for (x in seq_along(sample.n)){
      ref.aligned = as.character(trim.align[sample.n[x]])
      split.align = unlist(strsplit(as.character(ref.aligned), ""))

      #Replace leading gaps with N
      for (y in seq_along(split.align)){
        if (split.align[y] != "-"){ break }
        split.align[y] = "N"
      }#end forward loop

      #Replace trailing gaps with N
      for (y in length(split.align):1){
        if (split.align[y] != "-"){ break }
        split.align[y] = "N"
      }#end reverse loop

      align.list[[x]] = split.align
      names(align.list)[x] = sample.n[x]
    }#end x loop

    aligned.set = as.matrix(ape::as.DNAbin(align.list))
    PhyloProcessR::writePhylip(aligned.set, file = paste0(output.dir, "/untrimmed-alignments/", names(ref.data)[i], ".phy"), interleave = FALSE)
    file.remove(paste0(output.dir, "/unaligned-markers/", names(ref.data)[i], "_align.fa"))

  }#end i loop

  return(invisible(NULL))

} #end function
