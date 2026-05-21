#' @title prepareAlignments
#'
#' @description Prepares per-sample marker FASTA files for alignment by
#'   BLASTing draft contigs against the reference marker database, extracting
#'   matched locus regions, and handling fragmented matches by merging with
#'   CAP3 or joining with N spacers. Output per-sample FASTA files are written
#'   to \code{Alignments/sample-fastas/} for subsequent use by
#'   \code{markerAlignment}.
#'
#' @param contig.folder path to a folder of draft contig FASTA files (one .fa
#'   per sample).
#'
#' @param genbank.file path to a GenBank annotation file used by
#'   \code{makeReference} to build the reference marker database.
#'
#' @param blast.path path to the directory containing the BLAST executables.
#'   NULL uses the system PATH.
#'
#' @param overwrite logical; if TRUE, existing output directories are deleted
#'   and recreated.
#'
#' @param quiet logical; if TRUE, BLAST screen output is suppressed.
#'
#' @return Invisibly returns NULL. Writes per-sample marker FASTA files to
#'   \code{Alignments/sample-fastas/}.
#'
#' @examples
#' \dontrun{
#' prepareAlignments(
#'   contig.folder = "draftContigs",
#'   genbank.file  = "Crocidura.gb",
#'   blast.path    = "/path/to/blast/bin",
#'   overwrite     = FALSE,
#'   quiet         = TRUE
#' )
#' }
#'
#' @export

prepareAlignments = function(contig.folder = NULL,
                             genbank.file = NULL,
                             blast.path = NULL,
                             overwrite = FALSE,
                             quiet = TRUE) {

  if (!is.null(blast.path)){
    b.string = unlist(strsplit(blast.path, ""))
    if (b.string[length(b.string)] != "/") {
      blast.path = paste0(append(b.string, "/"), collapse = "")
    }
  } else { blast.path = "" }

  if (dir.exists("Alignments") == FALSE) {
    dir.create("Alignments")
  } else if (overwrite) {
    unlink("Alignments", recursive = TRUE)
    dir.create("Alignments")
  } else { message("Alignments directory already exists and overwrite = FALSE. Skipping."); return(invisible(NULL)) }

  if (dir.exists("Alignments/sample-fastas") == FALSE) {
    dir.create("Alignments/sample-fastas")
  } else if (overwrite) {
    unlink("Alignments/sample-fastas", recursive = TRUE)
    dir.create("Alignments/sample-fastas")
  }

  spp.samples = list.files(contig.folder)
  spp.samples = gsub(".fa$", "", spp.samples)

  makeReference(genbank.file = genbank.file, overwrite = overwrite)

  headers = c("qName", "tName", "pident", "matches", "misMatches", "gapopen",
              "qStart", "qEnd", "tStart", "tEnd", "evalue", "bitscore", "qLen", "tLen", "gaps")

  system(paste0(blast.path, "makeblastdb -in Mito-Reference/refMarkers.fa -parse_seqids -dbtype nucl",
                " -out mito-blast_db"), ignore.stdout = quiet, ignore.stderr = quiet)

  for (i in seq_along(spp.samples)){

    contigs = Rsamtools::scanFa(Rsamtools::FaFile(paste0(contig.folder, "/", spp.samples[i], ".fa")))

    system(paste0(blast.path, "blastn -task dc-megablast -db mito-blast_db",
                  " -query ", contig.folder, "/", spp.samples[i], ".fa",
                  " -out ", spp.samples[i], "_match.txt",
                  " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen gaps\" ",
                  " -num_threads 1"),
           ignore.stdout = quiet, ignore.stderr = quiet)

    match.file = paste0(spp.samples[i], "_match.txt")
    if (!file.exists(match.file) || file.info(match.file)$size == 0) {
      message(spp.samples[i], ": no BLAST output. Skipping.")
      next
    }

    match.data = read.table(match.file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    colnames(match.data) = headers
    file.remove(match.file)

    if (nrow(match.data) == 0){
      message(spp.samples[i], ": no matching mitochondrial genes were found.")
      next
    }

    filt.data = match.data[match.data$matches > 8, ]
    filt.data = filt.data[filt.data$evalue <= 0.05, ]

    neg.strand = filt.data$tStart > filt.data$tEnd
    filt.data$qDir = ifelse(neg.strand, "-", "+")
    filt.data[neg.strand, c("tStart", "tEnd")] = filt.data[neg.strand, c("tEnd", "tStart")]

    loci.names = unique(filt.data$tName)
    good.data  = Biostrings::DNAStringSet()
    good.match = c()

    for (j in seq_along(loci.names)){

      sub.data = filt.data[filt.data$tName %in% loci.names[j], ]
      sub.data = sub.data[order(sub.data$tStart, decreasing = FALSE), ]

      if (nrow(sub.data) == 1){
        temp.contig = contigs[names(contigs) %in% sub.data$qName]
        start = sub.data$qStart - (sub.data$tStart - 1)
        if (start <= 0) { start = 1 }
        end = sub.data$qEnd + (sub.data$tLen - sub.data$tEnd)
        if (end > sub.data$qLen) { end = sub.data$qLen }
        new.contig = Biostrings::subseq(temp.contig, start = start, end = end)
        names(new.contig) = loci.names[j]
        good.data  = append(good.data, new.contig)
        good.match = rbind(good.match, sub.data)
        next
      }

      complete.data = sub.data[sub.data$matches >= sub.data$tLen * 0.70, ]

      if (nrow(complete.data) == 1){
        temp.contig = contigs[names(contigs) %in% complete.data$qName]
        start = complete.data$qStart - (complete.data$tStart - 1)
        if (start <= 0) { start = 1 }
        end = complete.data$qEnd + (complete.data$tLen - complete.data$tEnd)
        if (end > complete.data$qLen) { end = complete.data$qLen }
        new.contig = Biostrings::subseq(temp.contig, start = start, end = end)
        names(new.contig) = loci.names[j]
        good.data  = append(good.data, new.contig)
        good.match = rbind(good.match, complete.data)
        next
      }

      if (nrow(complete.data) >= 2){
        temp.data = complete.data[complete.data$bitscore == max(complete.data$bitscore), ][1, ]
        temp.contig = contigs[names(contigs) %in% temp.data$qName]
        start = temp.data$qStart - (temp.data$tStart - 1)
        if (start <= 0) { start = 1 }
        end = temp.data$qEnd + (temp.data$tLen - temp.data$tEnd)
        if (end > temp.data$qLen) { end = temp.data$qLen }
        new.contig = Biostrings::subseq(temp.contig, start = start, end = end)
        names(new.contig) = loci.names[j]
        good.data  = append(good.data, new.contig)
        good.match = rbind(good.match, temp.data)
        next
      }

      cut.contigs = Biostrings::DNAStringSet()
      for (k in seq_len(nrow(sub.data))){
        temp.contig = contigs[names(contigs) %in% sub.data$qName[k]]
        new.contig  = Biostrings::subseq(temp.contig, start = sub.data$qStart[k], end = sub.data$qEnd[k])
        cut.contigs = append(cut.contigs, new.contig)
      }

      cap3.contigs = MitoTrawlR::runCap3(contigs = cut.contigs, read.R = TRUE)

      if (length(cap3.contigs) == 1){
        new.contig = cap3.contigs
        names(new.contig) = loci.names[j]
        good.data  = append(good.data, new.contig)
        good.match = rbind(good.match, sub.data)
        next
      } else {

        sub.data2 = sub.data
        sub.data  = sub.data2
        index = 1
        for (k in seq_len(nrow(sub.data) - 1)){
          if (nrow(sub.data) - 1 == k){ break }
          if (sub.data$tStart[index + 1] <= sub.data$tEnd[index]){
            overlap = sub.data$tEnd[index] - sub.data$tStart[index + 1]
            if (overlap > sub.data$matches[index + 1]){
              sub.data = sub.data[!sub.data$matches == min(sub.data$matches[index + 1], sub.data$matches[index]), ]
              index = index - 1
            }
          }
          if (nrow(sub.data) - 1 == k){ break }
          index = index + 1
        }

        paste.contigs = Biostrings::DNAStringSet()
        for (k in seq_len(nrow(sub.data) - 1)){
          if (sub.data$tStart[k + 1] <= sub.data$tEnd[k]){

            overlap = sub.data$tEnd[k] - sub.data$tStart[k + 1]
            temp.contig1 = contigs[names(contigs) %in% sub.data$qName[k]]
            temp.contig2 = contigs[names(contigs) %in% sub.data$qName[k + 1]]

            if (overlap > sub.data$matches[k + 1]){
              sub.data = sub.data[!sub.data$matches == min(sub.data$matches[k + 1], sub.data$matches[k]), ]
            }
            sub.data$tEnd[k] = sub.data$tEnd[k] - overlap - 1
            sub.data$qEnd[k] = sub.data$qEnd[k] - overlap - 1
            if (sub.data$qStart[k] >= sub.data$qEnd[k]){ next }

            if (length(paste.contigs) == 0){
              new.contig1   = Biostrings::subseq(temp.contig1, start = sub.data$qStart[k],     end = sub.data$qEnd[k])
              new.contig2   = Biostrings::subseq(temp.contig2, start = sub.data$qStart[k + 1], end = sub.data$qEnd[k + 1])
              paste.contigs = Biostrings::DNAStringSet(paste0(as.character(new.contig1), as.character(new.contig2)))
              names(paste.contigs) = loci.names[j]
            } else {
              new.contig2   = Biostrings::subseq(temp.contig2, start = sub.data$qStart[k + 1], end = sub.data$qEnd[k + 1])
              paste.contigs = Biostrings::DNAStringSet(paste0(as.character(paste.contigs), as.character(new.contig2)))
              names(paste.contigs) = loci.names[j]
            }

          } else {

            gap  = sub.data$tStart[k + 1] - sub.data$tEnd[k]
            temp.contig1 = contigs[names(contigs) %in% sub.data$qName[k]]
            temp.contig2 = contigs[names(contigs) %in% sub.data$qName[k + 1]]
            gap.n = paste0(rep("N", gap), collapse = "")
            if (sub.data$qStart[k] >= sub.data$qEnd[k]){ next }

            if (length(paste.contigs) == 0){
              new.contig1   = Biostrings::subseq(temp.contig1, start = sub.data$qStart[k],     end = sub.data$qEnd[k])
              new.contig2   = Biostrings::subseq(temp.contig2, start = sub.data$qStart[k + 1], end = sub.data$qEnd[k + 1])
              paste.contigs = Biostrings::DNAStringSet(paste0(as.character(new.contig1), gap.n, as.character(new.contig2)))
              names(paste.contigs) = loci.names[j]
            } else {
              new.contig2   = Biostrings::subseq(temp.contig2, start = sub.data$qStart[k + 1], end = sub.data$qEnd[k + 1])
              paste.contigs = Biostrings::DNAStringSet(paste0(as.character(paste.contigs), gap.n, as.character(new.contig2)))
              names(paste.contigs) = loci.names[j]
            }
          }
        }

        if (length(paste.contigs) == 0){ next }
        good.data  = append(good.data, paste.contigs)
        good.match = rbind(good.match, sub.data)

      }

    }#end j loop

    names(good.data) = paste0(spp.samples[i], "_|_", names(good.data))
    good.data = good.data[order(names(good.data))]
    write.loci = as.list(as.character(good.data))
    PhyloProcessR::writeFasta(sequences = write.loci, names = names(write.loci),
               paste0("Alignments/sample-fastas/", spp.samples[i], "_sampleFasta.fa"),
               nbchar = 1000000, as.string = TRUE)

  }#end i loop

  file.remove(list.files(pattern = "^mito-blast_db"))

  return(invisible(NULL))

}#end function
