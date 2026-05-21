#' @title buildMitogenomes
#'
#' @description Constructs final mitochondrial genome sequences for each sample
#'   using per-locus alignments and annotated sample contigs. For samples where
#'   all markers fall on a single contig the contig is reoriented to match the
#'   reference gene order; for fragmented samples, individual contigs are
#'   oriented and written as separate sequences. A per-sample summary CSV
#'   (completeness, N-count, fragment count, etc.) is written to
#'   \code{genome.dir}.
#'
#' @param annotation.dir path to the annotations directory produced by
#'   \code{annotateMitoContigs}, which must contain
#'   \code{sample-contigs/} and \code{sample-markers/} subdirectories.
#'
#' @param alignment.folder path to a folder of per-locus phylip alignments
#'   (untrimmed), as produced by \code{markerAlignment}.
#'
#' @param genome.alignment path to the concatenated whole-mitogenome phylip
#'   alignment produced by \code{alignMitogenomes}.
#'
#' @param genome.dir path to the top-level output directory for final genomes
#'   (e.g., \code{"MitoGenomes"}).
#'
#' @param reference.name name of the reference directory created by
#'   \code{buildReference}.
#'
#' @param output.dir subdirectory name (inside \code{genome.dir}) where
#'   finished genome FASTA files are written.
#'
#' @param blast.path path to the directory containing the BLAST executables.
#'   NULL uses the system PATH.
#'
#' @param threads number of CPU threads to pass to BLAST.
#'
#' @param overwrite logical; if TRUE, the output directory is deleted and
#'   recreated.
#'
#' @return Invisibly returns NULL. Writes per-sample genome FASTA files to
#'   \code{genome.dir/final-genomes/} and a summary CSV to \code{genome.dir}.
#'
#' @examples
#' \dontrun{
#' buildMitogenomes(
#'   annotation.dir   = "Annotations",
#'   alignment.folder = "Alignments/untrimmed-alignments",
#'   genome.alignment = "MitoGenomes/alignments/untrimmed_mitogenome_alignment.phy",
#'   genome.dir       = "MitoGenomes",
#'   reference.name   = "reference",
#'   output.dir       = "untrimmed-finished",
#'   overwrite        = FALSE
#' )
#' }
#'
#' @export

buildMitogenomes = function(annotation.dir = "Annotations",
                            alignment.folder = "Alignments/untrimmed-alignments",
                            genome.alignment = "MitoGenomes/alignments/untrimmed_mitogenome_alignment.phy",
                            genome.dir = "MitoGenomes",
                            reference.name = "reference",
                            output.dir = NULL,
                            blast.path = NULL,
                            threads = 1,
                            overwrite = FALSE) {

  #Debug
  # alignment.folder = "Alignments/untrimmed-alignments"
  # genome.dir = "MitoGenomes"
  # annotation.dir = "Annotations"
  # genome.alignment = "MitoGenomes/alignments/untrimmed_mitogenome_alignment.phy"
  # output.dir = "untrimmed-finished"
  # overwrite = FALSE
  # reference.name = "reference"

  if (is.null(blast.path) == FALSE){
    b.string = unlist(strsplit(blast.path, ""))
    if (b.string[length(b.string)] != "/") {
      blast.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { blast.path = "" }

  #Prepares output directories
  output.dir = paste0(genome.dir, "/", output.dir)
  if (dir.exists(output.dir) == FALSE) {
    dir.create(output.dir)
  } else if (overwrite == TRUE) {
    unlink(output.dir, recursive = TRUE)
    dir.create(output.dir)
  } else {
    message("Output directory '", output.dir, "' already exists and overwrite = FALSE. Skipping."); return(invisible(NULL))
  }

  #Create new directories
  dir.create(paste0(genome.dir, "/reference-order"), showWarnings = FALSE)
  dir.create(paste0(genome.dir, "/sample-markers"), showWarnings = FALSE)

  #Read in reference and alignment
  ref.data = Biostrings::readDNAStringSet(paste0(reference.name, "/refMarkers.fa"))
  gen.align = Biostrings::DNAStringSet(Biostrings::readAAMultipleAlignment(file = genome.alignment, format = "phylip"))
  sample.names = names(gen.align)
  locus.names = names(ref.data)
  first.marker = names(ref.data)[1]

  #Header data for features and whatnot
  header.data = c("Sample", "Circular", "Fragments", "bp_Complete", "Total_Length",
                  "bp_Count", "N_Count", "Marker_Complete", "Marker_Count",
                  "Missing_Markers", "Missing_BP")

  collect.data = data.table::data.table(matrix(as.numeric(0), nrow = length(sample.names), ncol = length(header.data)))
  data.table::setnames(collect.data, header.data)
  collect.data[, Sample := as.character(sample.names)]

  #Loops through each sample to assemble genomes
  for (i in seq_along(sample.names)){

    contigs = Biostrings::readDNAStringSet(paste0(annotation.dir, "/sample-contigs/", sample.names[i], "_sampleContigs.fa"))
    n.fragments = length(contigs)

    ##################################
    # Saves reference-order sequence
    ##################################
    sample.align = gen.align[names(gen.align) %in% sample.names[i]]
    stand.order = as.character(sample.align)
    PhyloProcessR::writeFasta(sequences = as.list(stand.order), names = names(stand.order),
               paste0("MitoGenomes/reference-order/", sample.names[i], "_referenceOrder.fa"),
               nbchar = 1000000, as.string = TRUE)

    ##################################
    # Saves per-marker sequences in locus order
    ##################################
    sample.data = Biostrings::DNAStringSet()
    for (j in seq_along(locus.names)){
      align = Biostrings::DNAStringSet(Biostrings::readAAMultipleAlignment(
        file = paste0(alignment.folder, "/", locus.names[j], ".phy"), format = "phylip"))
      new.align = align[names(align) == sample.names[i]]
      if (length(new.align) == 0){ next }
      names(new.align) = locus.names[j]
      sample.data = append(sample.data, new.align)
    }#end j loop

    final.save = as.list(as.character(sample.data))
    final.save = lapply(final.save, function(x) gsub("-", "", x))
    delete.sample = names(final.save[final.save == ""])
    final.write = final.save[!names(final.save) %in% delete.sample]

    PhyloProcessR::writeFasta(sequences = final.write, names = names(final.write),
               paste0(genome.dir, "/sample-markers/", sample.names[i], "_Markers.fa"),
               nbchar = 1000000, as.string = TRUE)

    #Make blast database for the contig sequences
    headers = c("qName", "tName", "pident", "matches", "misMatches", "gapopen",
                "qStart", "qEnd", "tStart", "tEnd", "evalue", "bitscore", "qLen", "tLen", "gaps")

    system(paste0(blast.path, "makeblastdb -in ", genome.dir, "/sample-markers/", sample.names[i],
                  "_Markers.fa -parse_seqids -dbtype nucl -out ", genome.dir, "/sample_blast_db"),
           ignore.stdout = TRUE, ignore.stderr = TRUE)

    match.file = paste0(sample.names[i], "_match.txt")
    system(paste0(blast.path, "blastn -task dc-megablast -db ", genome.dir, "/sample_blast_db",
                  " -query ", annotation.dir, "/sample-contigs/", sample.names[i], "_sampleContigs.fa",
                  " -out ", match.file,
                  " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen gaps\"",
                  " -num_threads ", threads),
           ignore.stdout = TRUE, ignore.stderr = TRUE)

    if (!file.exists(match.file) || file.info(match.file)$size == 0){
      message(sample.names[i], ": no BLAST output. Skipping.")
      next
    }

    match.data = data.table::fread(match.file, sep = "\t", header = FALSE)
    data.table::setnames(match.data, headers)
    file.remove(match.file)

    if (nrow(match.data) == 0){
      message(sample.names[i], ": no matching mitochondrial genes found. Skipping.")
      next
    }#end if

    ##############################################################
    # Part A: Fix up the match table
    ##############################################################
    filt.data = match.data[match.data$matches > 10,]
    filt.data = filt.data[filt.data$evalue <= 0.05,]
    filt.data = filt.data[filt.data$pident >= 95,]
    filt.data = filt.data[order(filt.data$qStart),]

    #Vectorized direction fix
    neg.strand = filt.data$tStart > filt.data$tEnd
    filt.data$qDir = ifelse(neg.strand, "-", "+")
    filt.data[neg.strand, c("tStart", "tEnd")] = filt.data[neg.strand, c("tEnd", "tStart")]

    #If there is only one fragment, reorganize to match reference gene order
    if (length(unique(filt.data$qName)) == 1){

      first.written = names(final.write)[1]
      first.hits = filt.data[filt.data$tName == first.written,]

      if (nrow(first.hits) == 0){
        warning(sample.names[i], ": first marker absent from BLAST hits; cannot orient single-fragment contig. Skipping.")
        next
      }

      # Majority-vote direction: robust when a marker has hits on both strands
      fwd.votes = sum(first.hits$qDir == "+")
      rev.votes = sum(first.hits$qDir == "-")
      temp.contig = contigs[names(contigs) == filt.data$qName[1]]

      if (fwd.votes >= rev.votes){
        start.pos = first.hits$qStart[which(first.hits$qDir == "+")[1]]
        first.part  = Biostrings::subseq(temp.contig, start = start.pos, end = Biostrings::width(temp.contig))
        second.part = Biostrings::subseq(temp.contig, start = 1, end = start.pos - 1)
        combined.contig = Biostrings::DNAStringSet(paste0(as.character(first.part), as.character(second.part)))
      } else {
        start.pos = first.hits$qEnd[which(first.hits$qDir == "-")[1]] + 1
        first.part  = Biostrings::subseq(temp.contig, start = start.pos, end = Biostrings::width(temp.contig))
        second.part = Biostrings::subseq(temp.contig, start = 1, end = start.pos - 1)
        combined.contig = Biostrings::reverseComplement(
          Biostrings::DNAStringSet(paste0(as.character(first.part), as.character(second.part))))
      }

      save.contig = as.list(as.character(combined.contig))
      names(save.contig) = paste0(sample.names[i], "_sequence-", seq_along(save.contig))
      PhyloProcessR::writeFasta(sequences = save.contig, names = names(save.contig),
                 paste0(output.dir, "/", sample.names[i], "_Complete.fa"),
                 nbchar = 1000000, as.string = TRUE)

    } else {

      frag.contig = Biostrings::DNAStringSet()
      f.names = unique(filt.data$qName)

      # Sort fragments by their earliest reference gene index so the output
      # follows the reference gene order, not the arbitrary contig assembly order
      ref.pos = sapply(f.names, function(fn) {
        fn.markers = filt.data$tName[filt.data$qName == fn]
        min(match(fn.markers, locus.names), na.rm = TRUE)
      })
      f.names = f.names[order(ref.pos)]

      for (k in seq_along(f.names)){
        temp.filt = filt.data[filt.data$qName == f.names[k],]
        temp.contig = contigs[names(contigs) == f.names[k]]
        # Majority-vote direction per fragment
        fwd.k = sum(temp.filt$qDir == "+")
        rev.k = sum(temp.filt$qDir == "-")
        if (fwd.k >= rev.k){
          frag.contig = append(frag.contig, temp.contig)
        } else {
          frag.contig = append(frag.contig, Biostrings::reverseComplement(temp.contig))
        }
      }#end k loop

      save.contig = as.list(as.character(frag.contig))
      names(save.contig) = paste0(sample.names[i], "_sequence-", seq_along(save.contig))
      PhyloProcessR::writeFasta(sequences = save.contig, names = names(save.contig),
                 paste0(output.dir, "/", sample.names[i], ".fa"),
                 nbchar = 1000000, as.string = TRUE)

    }#end fragment if

    ### Summary data collection
    data.table::set(collect.data, i = as.integer(i), j = match("Sample", header.data), value = sample.names[i])

    no.gap = gsub("-", "", stand.order)
    total.nogap = nchar(no.gap)
    total.bp = nchar(gsub("N", "", no.gap))
    total.n = total.nogap - total.bp

    data.table::set(collect.data, i = as.integer(i), j = match("Fragments", header.data), value = n.fragments)
    data.table::set(collect.data, i = as.integer(i), j = match("bp_Complete", header.data), value = round(total.bp / total.nogap, 3))
    data.table::set(collect.data, i = as.integer(i), j = match("Total_Length", header.data), value = total.nogap)
    data.table::set(collect.data, i = as.integer(i), j = match("bp_Count", header.data), value = total.bp)
    data.table::set(collect.data, i = as.integer(i), j = match("N_Count", header.data), value = total.n)

    total.markers = length(final.save)
    total.count = length(final.write)

    data.table::set(collect.data, i = as.integer(i), j = match("Marker_Complete", header.data), value = round(total.count / total.markers, 3))
    data.table::set(collect.data, i = as.integer(i), j = match("Marker_Count", header.data), value = total.count)
    data.table::set(collect.data, i = as.integer(i), j = match("Missing_Markers", header.data), value = total.markers - total.count)
    data.table::set(collect.data, i = as.integer(i), j = match("Missing_BP", header.data), value = round(total.n / total.nogap, 3))

  }#end i loop

  #Clean up BLAST database files
  blast.db.files = list.files(genome.dir, pattern = "sample_blast_db", full.names = TRUE)
  if (length(blast.db.files) > 0){ file.remove(blast.db.files) }

  write.csv(collect.data, file = paste0(genome.dir, "/sample_mito_genome_summary.csv"), row.names = FALSE)

  return(invisible(NULL))

}#end function
