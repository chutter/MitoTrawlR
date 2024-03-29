#' @title prepareAlignments
#'
#' @description Function for running the program spades to assemble short read sequencing data
#'
#' @param contig.folder path to a folder of sequence alignments in phylip format.
#'
#' @param genbank.file contigs are added into existing alignment if algorithm is "add"
#'
#' @param out.dir contigs are added into existing alignment if algorithm is "add"
#'
#' @param blast.path algorithm to use: "add" add sequences with "add.contigs"; "localpair" for local pair align. All others available
#'
#' @param overwrite TRUE to supress screen output
#'
#' @param quiet TRUE to supress screen output
#'
#' @return an alignment of provided sequences in DNAStringSet format. Also can save alignment as a file with save.name
#'
#' @examples
#'
#' your.tree = ape::read.tree(file = "file-path-to-tree.tre")
#' astral.data = astralPlane(astral.tree = your.tree,
#'                           outgroups = c("species_one", "species_two"),
#'                           tip.length = 1)
#'
#'
#' @export

prepareAlignments = function(contig.folder = NULL,
                             genbank.file = NULL,
                             blast.path = "blast",
                             overwrite = FALSE,
                             quiet = TRUE) {

  # #Debug
  # genbank.file = "Crocidura.gb"
  # contig.folder = "draftContigs"
  # overwrite = TRUE
  # quiet = TRUE
  # blast.path = "blast"

  if (dir.exists("Alignments") == FALSE) { dir.create("Alignments") }
  if (dir.exists("Alignments") == TRUE) {
    if (overwrite == TRUE){
      system(paste0("rm -r ", "Alignments"))
      dir.create("Alignments")
    } else { stop("overwrite is false and directory exists. Exiting.") }
  }#end dir exists

  if (dir.exists("Alignments/sample-fastas") == FALSE) { dir.create("Alignments/sample-fastas") }
  if (dir.exists("Alignments/sample-fastas") == TRUE) {
    if (overwrite == TRUE){
      system(paste0("rm -r ", "Alignments/sample-fastas"))
      dir.create("Alignments/sample-fastas")
    } else { stop("overwrite is false and directory exists. Exiting.") }
  }#end dir exists

  #Obtains samples
  spp.samples =  list.files(contig.folder)
  spp.samples = gsub(".fa$", "", spp.samples)

  #Makes reference for blasting
  makeReference(genbank.file = genbank.file,
                overwrite = overwrite)

  headers = c("qName", "tName", "pident", "matches", "misMatches", "gapopen",
              "qStart", "qEnd", "tStart", "tEnd", "evalue", "bitscore", "qLen", "tLen", "gaps")

  system(paste0("makeblastdb -in Mito-Reference/refMarkers.fa -parse_seqids -dbtype nucl",
                " -out mito-blast_db"), ignore.stdout = quiet, ignore.stderr = quiet)

  for (i in 1:length(spp.samples)){

    #Load in the data
    contigs = Rsamtools::scanFa(Rsamtools::FaFile(paste0(contig.folder, "/", spp.samples[i], ".fa")))   # loads up fasta file

    #Matches samples to loci
    system(paste0("blastn -task dc-megablast -db mito-blast_db",
                  " -query ", contig.folder, "/", spp.samples[i], ".fa",
                  " -out ", spp.samples[i], "_match.txt",
                  " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen gaps\" ",
                  " -num_threads 1"))

    #Need to load in transcriptome for each species and take the matching transcripts to the database
    match.data = read.table(paste0(spp.samples[i], "_match.txt"), sep = "\t", header = F, stringsAsFactors = FALSE)
    colnames(match.data) = headers
    if (nrow(match.data) == 0){
      print("No matching mitochondrial genes were found.")
      system(paste0("rm ", spp.samples[i], "_match.txt"))
      next
    }#end if

    ##############################################################
    # Part A: Fixes up the match database
    ##############################################################

    #Gets rid of very poor matches
    filt.data = match.data[match.data$matches > 8,]
    filt.data = filt.data[filt.data$evalue <= 0.05,]

    #Fixes direction and adds into data
    filt.data$qDir = as.character("0")
    #Finds out if they are overlapping
    for (k in 1:nrow(filt.data)){
      if (filt.data$tStart[k] > filt.data$tEnd[k]){
        filt.data$qDir[k]<-"-"
        new.start<-min(filt.data$tStart[k], filt.data$tEnd[k])
        new.end<-max(filt.data$tStart[k], filt.data$tEnd[k])
        filt.data$tStart[k]<-new.start
        filt.data$tEnd[k]<-new.end
      } else { filt.data$qDir[k]<-"+" }
    }#end k loop

    #Filters and saves
    loci.names = unique(filt.data$tName)
    good.data = Biostrings::DNAStringSet()
    good.match = c()
    for (j in 1:length(loci.names)){
      #pulls out data that matches to multiple contigs
      sub.data = filt.data[filt.data$tName %in% loci.names[j],]
      sub.data = sub.data[order(sub.data$tStart, decreasing = F),]

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

      #Fragmented across contigs
      cut.contigs = Biostrings::DNAStringSet()
      for (k in 1:nrow(sub.data)){

        temp.contig = contigs[names(contigs) %in% sub.data$qName[k]]
        new.contig = Biostrings::subseq(temp.contig, start = sub.data$qStart[k], end = sub.data$qEnd[k])
        cut.contigs = append(cut.contigs, new.contig)
      }

      cap3.contigs = runCap3(cut.contigs)
      #Saves if there is one contig
      if (length(cap3.contigs) == 1){
        new.contig = cap3.contigs
        names(new.contig) = loci.names[j]
        good.data = append(good.data, new.contig)
        good.match = rbind(good.match, sub.data)
        next
      } else {

        sub.data2 = sub.data
        sub.data = sub.data2
        index = 1
        for (k in 1:(nrow(sub.data)-1)){

          if (nrow(sub.data)-1 == k){ break }

          #If they are overlapping
          if (sub.data$tStart[index+1] <= sub.data$tEnd[index]){

            overlap = sub.data$tEnd[index] - sub.data$tStart[index+1]
            #Remove smaller overlap
            if (overlap > sub.data$matches[index+1]){
              sub.data = sub.data[!sub.data$matches == min(sub.data$matches[index+1], sub.data$matches[index]),]
              index = index-1
            }

          }#end if
          if (nrow(sub.data)-1 == k){ break }
          index = index + 1
        }#end k

        paste.contigs = Biostrings::DNAStringSet()
        for (k in 1:(nrow(sub.data)-1)){
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
              paste.contigs = Biostrings::DNAStringSet(paste0(as.character(paste.contigs), as.character(new.contig)))
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

    #Writes the full mitochondrial genome file
    names(good.data) = paste0(spp.samples[i], "_|_", names(good.data))
    good.data = good.data[order(names(good.data)),]
    write.loci = as.list(as.character(good.data))
    PhyloProcessR::writeFasta(sequences = write.loci, names = names(write.loci),
               paste0("Alignments/sample-fastas/", spp.samples[i], "_sampleFasta.fa"), nbchar = 1000000, as.string = T)

    system(paste0("rm ", spp.samples[i], "_match.txt"))

  }#End i loop

  system(paste0("rm mito-blast_db*"))

}#end function
