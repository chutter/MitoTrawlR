#' @title annotateMitoContigs
#'
#' @description Annotates draft mitochondrial contigs for each sample by BLASTing
#'   them against the reference marker database, extracting matched gene
#'   regions, merging fragmented hits with CAP3, and running tRNAscan-SE to
#'   locate tRNA genes. Per-sample annotation summaries, individual marker
#'   FASTA files, and contig FASTA files are written to an \code{Annotations/}
#'   directory.
#'
#' @param contig.folder path to a folder of draft contig FASTA files (one .fa
#'   per sample), as produced by \code{mitochondrialCapture}.
#'
#' @param reference.name name of the reference directory created by
#'   \code{buildReference}, which must contain \code{refMarkers.fa}.
#'
#' @param blast.path path to the directory containing the BLAST executables
#'   (\code{blastn}, \code{makeblastdb}). NULL uses the system PATH.
#'
#' @param tRNAscan.path path to the directory containing \code{tRNAscan-SE}.
#'   NULL uses the system PATH.
#'
#' @param cap3.path path to the directory containing \code{cap3}. NULL uses
#'   the system PATH.
#'
#' @param organism.type tRNAscan-SE organism model: one of \code{"mammal"},
#'   \code{"vertebrate"}, or \code{"eukaryotic"}.
#'
#' @param threads number of CPU threads to pass to BLAST.
#'
#' @param overwrite logical; if TRUE, existing output directories are removed
#'   and recreated.
#'
#' @param save.gff logical; if TRUE, a GFF3 annotation file is written
#'   alongside the CSV summary for each sample to
#'   \code{Annotations/sample-summary/<sample>_annotation.gff3}.
#'
#' @param quiet logical; if TRUE, BLAST and tRNAscan-SE screen output is
#'   suppressed.
#'
#' @return Invisibly returns NULL. Writes per-sample annotation CSV files
#'   (and optionally GFF3 files), marker FASTA files, and contig FASTA files
#'   to \code{Annotations/}.
#'
#' @examples
#' \dontrun{
#' annotateMitoContigs(
#'   contig.folder  = "draftContigs",
#'   reference.name = "reference",
#'   blast.path     = "/path/to/blast/bin",
#'   tRNAscan.path  = "/path/to/trnascan/bin",
#'   organism.type  = "vertebrate",
#'   overwrite      = FALSE,
#'   quiet          = TRUE
#' )
#' }
#'
#' @export

#Annotates mitochondrial contigs
annotateMitoContigs = function(contig.folder = NULL,
                               reference.name = "reference",
                               blast.path = NULL,
                               tRNAscan.path = NULL,
                               cap3.path = NULL,
                               organism.type = c("mammal", "vertebrate", "eukaryotic"),
                               threads = 1,
                               overwrite = FALSE,
                               save.gff = FALSE,
                               quiet = TRUE) {

  # # # #Debug
  #  reference.name = "reference"
  #  contig.folder = "draftContigs"
  #  overwrite = TRUE
  #  quiet = FALSE
  #  organism.type = "vertebrate"

  organism.type = match.arg(organism.type)

  if (is.null(blast.path) == FALSE){
    b.string = unlist(strsplit(blast.path, ""))
    if (b.string[length(b.string)] != "/") {
      blast.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { blast.path = "" }

  if (is.null(tRNAscan.path) == FALSE){
    b.string = unlist(strsplit(tRNAscan.path, ""))
    if (b.string[length(b.string)] != "/") {
      tRNAscan.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { tRNAscan.path = "" }

  if (is.null(cap3.path) == FALSE){
    b.string = unlist(strsplit(cap3.path, ""))
    if (b.string[length(b.string)] != "/") {
      cap3.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { cap3.path = "" }

  #Checks for output directories
  if (dir.exists("Annotations") == FALSE) {
    dir.create("Annotations")
  } else if (overwrite == TRUE) {
    unlink("Annotations", recursive = TRUE)
    dir.create("Annotations")
  } else { message("Annotations directory already exists and overwrite = FALSE. Skipping."); return(invisible(NULL)) }

  if (dir.exists("Annotations/sample-contigs") == FALSE) {
    dir.create("Annotations/sample-contigs")
  } else if (overwrite == TRUE) {
    unlink("Annotations/sample-contigs", recursive = TRUE)
    dir.create("Annotations/sample-contigs")
  }#end dir exists

  if (dir.exists("Annotations/sample-markers") == FALSE) {
    dir.create("Annotations/sample-markers")
  } else if (overwrite == TRUE) {
    unlink("Annotations/sample-markers", recursive = TRUE)
    dir.create("Annotations/sample-markers")
  }#end dir exists

  if (dir.exists("Annotations/sample-summary") == FALSE) {
    dir.create("Annotations/sample-summary")
  } else if (overwrite == TRUE) {
    unlink("Annotations/sample-summary", recursive = TRUE)
    dir.create("Annotations/sample-summary")
  }#end dir exists

  #Obtains samples
  spp.samples =  list.files(contig.folder)
  spp.samples = gsub(".fa$", "", spp.samples)

  ref.data = Rsamtools::scanFa(Rsamtools::FaFile(paste0(reference.name, "/refMarkers.fa")))

  headers = c("qName", "tName", "pident", "matches", "misMatches", "gapopen",
              "qStart", "qEnd", "tStart", "tEnd", "evalue", "bitscore", "qLen", "tLen", "gaps")

  system(paste0(blast.path, "makeblastdb -in ", reference.name,
                "/refMarkers.fa -parse_seqids -dbtype nucl",
                " -out ", reference.name, "/ref-blast-db"),
         ignore.stdout = quiet, ignore.stderr = quiet)

  for (i in seq_along(spp.samples)){

    #Load in the data
    contigs = Biostrings::readDNAStringSet(paste0(contig.folder, "/", spp.samples[i], ".fa"))

    #Matches samples to loci
    system(paste0(blast.path, "blastn -task dc-megablast -db ", reference.name, "/ref-blast-db",
                  " -query ", contig.folder, "/", spp.samples[i], ".fa",
                  " -out ", spp.samples[i], "_match.txt",
                  " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen gaps\" ",
                  " -num_threads ", threads),
           ignore.stdout = quiet, ignore.stderr = quiet)

    match.file = paste0(spp.samples[i], "_match.txt")
    if (!file.exists(match.file) || file.info(match.file)$size == 0){
      message(spp.samples[i], ": no matching mitochondrial genes were found.")
      if (file.exists(match.file)) file.remove(match.file)
      next
    }

    match.data = data.table::fread(match.file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    colnames(match.data) = headers
    if (nrow(match.data) == 0){
      message(spp.samples[i], ": no matching mitochondrial genes were found.")
      file.remove(match.file)
      next
    }#end if

    ### Here reorient contigs, merge the same one?


    ##############################################################
    # Part A: Fixes up the match database
    ##############################################################

    #Gets rid of very poor matches
    filt.data = match.data[match.data$matches > 8,]
    filt.data = filt.data[filt.data$evalue <= 0.05,]
    filt.data = filt.data[grep("tRNA", filt.data$tName, invert = TRUE),]

    if (nrow(filt.data) == 0){
      message(spp.samples[i], ": no matching mitochondrial genes were found.")
      file.remove(paste0(spp.samples[i], "_match.txt"))
      next
    }#end if

    #Fixes direction and adds into data (vectorized)
    neg.strand = filt.data$tStart > filt.data$tEnd
    filt.data$qDir = ifelse(neg.strand, "-", "+")
    filt.data[neg.strand, c("tStart", "tEnd")] = filt.data[neg.strand, c("tEnd", "tStart")]

    #Filters and saves
    loci.names = unique(filt.data$tName)
    good.data = Biostrings::DNAStringSet()
    good.match = c()
    for (j in seq_along(loci.names)){
      #pulls out data that matches to multiple contigs
      sub.data = filt.data[filt.data$tName %in% loci.names[j],]
      sub.data = sub.data[order(sub.data$tStart, decreasing = FALSE),]

      #If there is only 1 match
      if (nrow(sub.data) == 1){
        #Subsets to desired contig
        temp.contig = contigs[names(contigs) %in% sub.data$qName]
        start = sub.data$qStart - (sub.data$tStart - 1)
        if (start <= 0) { start = 1 }
        end = sub.data$qEnd + (sub.data$tLen - sub.data$tEnd)
        if (end > sub.data$qLen) {end = sub.data$qLen }
        new.contig = Biostrings::subseq(temp.contig, start = start, end = end)
        names(new.contig) = loci.names[j]
        good.data = append(good.data, new.contig)
        good.match = rbind(good.match, sub.data)
        next
      }

      #Checks for the most complete match if more than 1 match
      complete.data = sub.data[sub.data$matches >= sub.data$tLen * 0.70,]

      if (nrow(complete.data) == 1){
        #Subsets to desired contig
        temp.contig = contigs[names(contigs) %in% complete.data$qName]
        start = complete.data$qStart - (complete.data$tStart - 1)
        if (start <= 0) { start = 1 }
        end = complete.data$qEnd + (complete.data$tLen - complete.data$tEnd)
        if (end > complete.data$qLen) {end = complete.data$qLen }
        new.contig = Biostrings::subseq(temp.contig, start = start, end = end)
        names(new.contig) = loci.names[j]
        good.data = append(good.data, new.contig)
        good.match = rbind(good.match, complete.data)
        next
      }#end if for complete data

      if (nrow(complete.data) >= 2){
        #Gets data with highest bitscore
        temp.data = complete.data[complete.data$bitscore == max(complete.data$bitscore),][1,]
        #Subsets to desired contig
        temp.contig = contigs[names(contigs) %in% temp.data$qName]
        start = temp.data$qStart - (temp.data$tStart - 1)
        if (start <= 0) { start = 1 }
        end = temp.data$qEnd + (temp.data$tLen - temp.data$tEnd)
        if (end > temp.data$qLen) {end = temp.data$qLen }
        new.contig = Biostrings::subseq(temp.contig, start = start, end = end)
        names(new.contig) = loci.names[j]
        good.data = append(good.data, new.contig)
        good.match = rbind(good.match, temp.data)
        next
      }
#
#       if (nrow(complete.data) == 0) {
#         sub.data
#
#       }

      #Fragmented across contigs
      cut.contigs = Biostrings::DNAStringSet()
      for (k in seq_len(nrow(sub.data))){

        temp.contig = contigs[names(contigs) %in% sub.data$qName[k]]
        new.contig = Biostrings::subseq(temp.contig, start = sub.data$qStart[k], end = sub.data$qEnd[k])
        cut.contigs = append(cut.contigs, new.contig)
      }

      cap3.contigs = MitoTrawlR::runCap3(contigs = cut.contigs,
                                         read.R = TRUE,
                                         cap3.path = cap3.path)
      #Saves if there is one contig
      if (length(cap3.contigs) == 1 && length(cap3.contigs) > 0){
        new.contig = cap3.contigs
        names(new.contig) = loci.names[j]
        good.data = append(good.data, new.contig)
        good.match = rbind(good.match, sub.data)
        next
      } else {

        index = 1
        for (k in seq_len(nrow(sub.data)-1)){

          if (nrow(sub.data) == 1){ break }

          #If they are overlapping
          if (sub.data$tStart[index+1] <= sub.data$tEnd[index]){

            overlap = sub.data$tEnd[index] - sub.data$tStart[index+1]
            #Remove smaller overlap
            if (overlap > sub.data$matches[index+1]){
              sub.data = sub.data[!sub.data$matches == min(sub.data$matches[index+1], sub.data$matches[index]),]
              index = index-1
            }

          }#end if
          index = index + 1
        }#end k

        paste.contigs = Biostrings::DNAStringSet()
        for (k in seq_len(nrow(sub.data)-1)){
          #If they are overlapping
          if (sub.data$tStart[k+1] <= sub.data$tEnd[k]){

            overlap = sub.data$tEnd[k] - sub.data$tStart[k+1]
            temp.contig1 = contigs[names(contigs) %in% sub.data$qName[k]]
            temp.contig2 = contigs[names(contigs) %in% sub.data$qName[k+1]]

            #Remove smaller overlap
            if (overlap > sub.data$matches[k+1]){
              sub.data = sub.data[!sub.data$matches == min(sub.data$matches[k+1], sub.data$matches[k]),]
            }

            sub.data$tEnd[k] = sub.data$tEnd[k] - overlap -1
            sub.data$qEnd[k] = sub.data$qEnd[k] - overlap -1

            if (sub.data$qStart[k] >= sub.data$qEnd[k]){ next }

            if (length(paste.contigs) == 0){
              new.contig1 = Biostrings::subseq(temp.contig1, start = sub.data$qStart[k], end = sub.data$qEnd[k])
              new.contig2 = Biostrings::subseq(temp.contig2, start = sub.data$qStart[k+1], end = sub.data$qEnd[k+1])
              paste.contigs = Biostrings::DNAStringSet(paste0(as.character(new.contig1), as.character(new.contig2)))
              names(paste.contigs) = loci.names[j]
            } else {
              new.contig2 = Biostrings::subseq(temp.contig2, start = sub.data$qStart[k+1], end = sub.data$qEnd[k+1])
              paste.contigs = Biostrings::DNAStringSet(paste0(as.character(paste.contigs), as.character(new.contig2)))
              names(paste.contigs) = loci.names[j]
            }#end else paste.contigs

          } else {
            #ADD IN NS
            gap =  sub.data$tStart[k+1] - sub.data$tEnd[k]

            temp.contig1 = contigs[names(contigs) %in% sub.data$qName[k]]
            temp.contig2 = contigs[names(contigs) %in% sub.data$qName[k+1]]
            gap.n = paste0(rep("N", gap), collapse = "")

            if (sub.data$qStart[k] >= sub.data$qEnd[k]){ next }

            if (length(paste.contigs) == 0){
              new.contig1 = Biostrings::subseq(temp.contig1, start = sub.data$qStart[k], end = sub.data$qEnd[k])
              new.contig2 = Biostrings::subseq(temp.contig2, start = sub.data$qStart[k+1], end = sub.data$qEnd[k+1])
              paste.contigs = Biostrings::DNAStringSet(paste0(as.character(new.contig1), as.character(gap.n),
                                                              as.character(new.contig2)))
              names(paste.contigs) = loci.names[j]
            } else {
              new.contig2 = Biostrings::subseq(temp.contig2, start = sub.data$qStart[k+1], end = sub.data$qEnd[k+1])
              paste.contigs = Biostrings::DNAStringSet(paste0(as.character(paste.contigs), as.character(gap.n),
                                                              as.character(new.contig2)))
              names(paste.contigs) = loci.names[j]
            }#end else paste.contigs
          }#end else overlap
        }#end k loop

        if (length(paste.contigs) == 0){ next }

        good.data = append(good.data, paste.contigs)
        good.match = rbind(good.match, sub.data)

      }#end else

    }#end j loop

    file.remove(paste0(spp.samples[i], "_match.txt"))

    # ######## Reposition marker
    # ##############################################
    # ##############################################
    #
    #
    # first.marker = sort(unique(good.match$tName))[1]
    # #if (nrow(filt.data[filt.data$tName == first.marker,]) >= 2){ stop("whoa") }
    #
    # if (unique(good.match[good.match$tName == first.marker,]$qDir) == "+"){
    #
    #   #start.marker = filt.data[filt.data$tName == first.marker,]
    #   #end.marker = filt.data[filt.data$tName == last.marker,]
    #   start.pos = good.match[good.match$tName == first.marker,]$qStart
    #   temp.contig = contigs[names(contigs) == good.match$qName[1]]
    #   first.part = Biostrings::subseq(temp.contig, start = start.pos[1], end = Biostrings::width(temp.contig))
    #   second.part = Biostrings::subseq(temp.contig, start = 1, end = start.pos[1]-1)
    #
    #   #Combine together
    #   combined.contig = Biostrings::DNAStringSet(paste0(as.character(first.part),
    #                                                     as.character(second.part)))
    #
    # }#end if
    #
    # if (unique(good.match[good.match$tName == first.marker,]$qDir) == "-"){
    #
    #   #start.marker = filt.data[filt.data$tName == last.marker,]
    #   #end.marker = filt.data[filt.data$tName == first.marker,]
    #
    #   start.pos = good.match[good.match$tName == first.marker,]$qEnd + 1
    #   temp.contig = contigs[names(contigs) == good.match$qName[1]]
    #   first.part = Biostrings::subseq(temp.contig, start = start.pos[1], end = Biostrings::width(temp.contig))
    #   second.part = Biostrings::subseq(temp.contig, start = 1, end = start.pos[1] - 1)
    #
    #   #Combine together
    #   combined.contig = Biostrings::DNAStringSet(paste0(as.character(first.part),
    #                                                     as.character(second.part)))
    #   #Reverse compliment
    #   #Reverses alignment back to correction orientation
    #   combined.contig = Biostrings::reverseComplement(combined.contig)
    #
    # }#end if
    #
    # #Adds reference locus
    # save.contig = as.list(as.character(combined.contig))
    # names(save.contig) = paste0(sample.names[i], "_sequence-", rep(1:length(save.contig)))
    # #Saves to folder the standard order one already made
    # writeFasta(sequences = save.contig, names = names(save.contig),
    #            paste0(genome.dir, "/final-genomes/", sample.names[i], "_Complete.fa"),
    #            nbchar = 1000000, as.string = TRUE)
    #
    #
    #

    ######## Final step tRNA annotation
    ##############################################
    ##############################################

    ### Annotate tRNAs
    rna.data = tRNAscan(contigs = good.data,
                        tRNAscan.path = tRNAscan.path,
                        organism.type = organism.type,
                        quiet = quiet)

    if (is.null(rna.data) == FALSE){ if (nrow(rna.data) == 0){ rna.data = NULL }}

    if (is.null(rna.data) == FALSE){

      # #Remove duplicate RNAs
      # if (length(unique(good.match$qName)) == 1){
      #   rna.data = rna.data[rna.data$contig %in% unique(good.match$qName),]
      # }

      rna.contigs = Biostrings::DNAStringSet()
      for (j in seq_len(nrow(rna.data))){

        new.contig = Biostrings::subseq(good.data[names(good.data) == rna.data$contig[j]],
                                         start = rna.data$start[j],
                                         end = rna.data$end[j])

        names(new.contig) = paste0("tRNA-", rna.data$type[j])
        rna.contigs = append(rna.contigs, new.contig)
      }#end j loop

      good.data = append(good.data, rna.contigs)
      rna.names = unique(rna.data$type)

      add.rna = data.frame(contig = rna.data$contig,
                           name = paste0("tRNA-", rna.data$type),
                           start = rna.data$start,
                           end = rna.data$end,
                           direction = rna.data$dir)

    } else { add.rna = data.frame() } #end null if

    # Fix later, negative starts when readjsting
    # refine.match = good.match[order(good.match$qStart),]
    #
    # for (x in 1:nrow(refine.match)) {
    #   if (refine.match$qDir[x] == "+"){
    #     refine.match$qStart[x] = refine.match$qStart[x] - (refine.match$tStart[x] - 1)
    #     refine.match$qEnd[x] = refine.match$qEnd[x] + (refine.match$tLen[x] - refine.match$tEnd[x])
    #   } #end + if
    #
    #   if (refine.match$qDir[x] == "-"){
    #     refine.match$qEnd[x] = refine.match$qEnd[x] + (refine.match$tStart[x] - 1)
    #     refine.match$qStart[x] = refine.match$qStart[x] - (refine.match$tLen[x] - refine.match$tEnd[x])
    #   } #end - if
    # }#end x loop
    #
    # #### TO DO:
      #### Fix D-loop, combine into 1 with large range and fix ending and beginning

      # #Sets up new data
      # new.data = data.frame(contig = as.character(),
      #                       feature = as.character(),
      #                       tStart = as.numeric(),
      #                       tEnd = as.numeric(),
      #                       bStart = as.numeric(),
      #                       bEnd = as.numeric(),
      #                       dir = as.numeric(),
      #                       blast = as.logical(),
      #                       tRNAScan = as.logical())
      #
      # for (x in 1:length(rna.names)){
      #
      #   temp.good = refine.match[grep(rna.names[x], refine.match$tName),]
      #   temp.rna = rna.data[rna.data$type == rna.names[x],]
      #
      #   #If there were no reference matches, uses tRNA match
      #   if (nrow(temp.good) == 0){
      #     #Add to the annotation file and check coordinates
      #     new.entry = data.frame(contig = temp.rna$contig,
      #                            feature = paste0("0X_tRNA_", temp.rna$type),
      #                            tStart = temp.rna$start,
      #                            tEnd = temp.rna$end,
      #                            bStart = NA,
      #                            bEnd = NA,
      #                            dir = temp.rna$dir,
      #                            blast = FALSE,
      #                            tRNAScan = TRUE)
      #     new.data = rbind(new.data, new.entry)
      #     next
      #   }#end if
      #
      #   #Deals with the tRNAs with the same name
      #   if (nrow(temp.good) >= 2){
      #     temp.rna = temp.rna[temp.rna$contig %in% temp.good$qName,]
      #     temp.good = temp.good[order(temp.good$qStart),]
      #     temp.rna = temp.rna[order(temp.rna$start),]
      #     #Add to the annotation file and check coordinates
      #     new.entry = data.frame(contig = temp.good$qName,
      #                            feature = temp.good$tName,
      #                            tStart = temp.good$tStart,
      #                            tEnd = temp.good$tEnd,
      #                            bStart = temp.good$qStart,
      #                            bEnd = temp.good$qEnd,
      #                            dir = temp.good$qDir,
      #                            blast = TRUE,
      #                            tRNAScan = TRUE)
      #     new.data = rbind(new.data, new.entry)
      #     next
      #   }
      #
      #   #Duplciate contigs, picks the filtered match set
      #   if (nrow(temp.rna) >= 2){
      #     temp.rna = temp.rna[temp.rna$contig %in% temp.good$qName,]
      #   }
      #
      #   #Even worse! have to pick the best
      #   if (nrow(temp.rna) >= 2){
      #     temp.rna = temp.rna[temp.rna$score == max(temp.rna$score),][1,]
      #     next }
      #
      #   #Add to the annotation file and check coordinates
      #   new.entry = data.frame(contig = temp.rna$contig,
      #                          feature = temp.good$tName,
      #                          tStart = temp.rna$start,
      #                          tEnd = temp.rna$end,
      #                          bStart = temp.good$qStart,
      #                          bEnd = temp.good$qEnd,
      #                          dir = temp.rna$dir,
      #                          blast = TRUE,
      #                          tRNAScan = TRUE)
      #   new.data = rbind(new.data, new.entry)
      #
      # }#End X


    refine.match = good.match
    add.cds = data.frame(contig = refine.match$qName,
                         name = refine.match$tName,
                         start = refine.match$qStart,
                         end = refine.match$qEnd,
                         direction = refine.match$qDir)

    final.sample = rbind(add.rna, add.cds)
    final.sample = final.sample[order(final.sample$contig, final.sample$start),]

    ### Writes the CSV summary
    write.csv(final.sample, file = paste0("Annotations/sample-summary/", spp.samples[i], "_sample-summary.csv"),
              row.names = FALSE)

    ### Optionally writes GFF3 annotation
    if (save.gff == TRUE) {
      gff.rows = data.frame(
        seqname    = final.sample$contig,
        source     = "MitoTrawlR",
        feature    = ifelse(grepl("^tRNA",          final.sample$name), "tRNA",
                     ifelse(grepl("rRNA|12S|16S",   final.sample$name, ignore.case = TRUE), "rRNA",
                     ifelse(grepl("D.loop|D_loop",  final.sample$name, ignore.case = TRUE), "D_loop",
                            "CDS"))),
        start      = final.sample$start,
        end        = final.sample$end,
        score      = ".",
        strand     = final.sample$direction,
        frame      = ifelse(grepl("^tRNA|rRNA|12S|16S|D.loop", final.sample$name, ignore.case = TRUE), ".", "0"),
        attributes = paste0("ID=", final.sample$name, ";Name=", final.sample$name),
        stringsAsFactors = FALSE
      )
      gff.file = paste0("Annotations/sample-summary/", spp.samples[i], "_annotation.gff3")
      writeLines("##gff-version 3", gff.file)
      write.table(gff.rows, file = gff.file, sep = "\t", quote = FALSE,
                  row.names = FALSE, col.names = FALSE, append = TRUE)
    }

    #Writes per-marker extracted sequences
    names(good.data) = paste0(spp.samples[i], "_|_", names(good.data))
    good.data = good.data[order(names(good.data))]
    write.loci = as.list(as.character(good.data))
    PhyloProcessR::writeFasta(sequences = write.loci, names = names(write.loci),
               paste0("Annotations/sample-markers/", spp.samples[i], "_sampleMarkers.fa"), nbchar = 1000000, as.string = TRUE)

    #Writes the original assembly contigs used in annotation
    good.contigs = contigs[names(contigs) %in% unique(good.match$qName)]
    write.loci = as.list(as.character(good.contigs))
    PhyloProcessR::writeFasta(sequences = write.loci, names = names(write.loci),
               paste0("Annotations/sample-contigs/", spp.samples[i], "_sampleContigs.fa"), nbchar = 1000000, as.string = TRUE)

    message("Finished annotation for ", spp.samples[i])

  }#End i loop

  message("Finished annotation for all samples!")
  return(invisible(NULL))

}#end function
