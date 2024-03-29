#' @title annotateMitoContigs
#'
#' @description Function for running the program spades to assemble short read sequencing data
#'
#' @param contig.folder path to a folder of sequence alignments in phylip format.
#'
#' @param genbank.file contigs are added into existing alignment if algorithm is "add"
#'
#' @param blast.path algorithm to use: "add" add sequences with "add.contigs"; "localpair" for local pair align. All others available
#'
#' @param tRNAscan.path algorithm to use: "add" add sequences with "add.contigs"; "localpair" for local pair align. All others available
#'
#' @param organism.type algorithm to use: "add" add sequences with "add.contigs"; "localpair" for local pair align. All others available
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

#Annotates mitochondrial contigs
annotateMitoContigs = function(contig.folder = NULL,
                               reference.name = "reference",
                               blast.path = NULL,
                               tRNAscan.path = NULL,
                               organism.type = c("mammal", "vertebrate", "eukaryotic"),
                               overwrite = FALSE,
                               quiet = TRUE) {

  # # # #Debug
  #  reference.name = "reference"
  #  contig.folder = "draftContigs"
  #  overwrite = TRUE
  #  quiet = FALSE
  #  organism.type = "vertebrate"

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

  #Checks for output directories
  if (dir.exists("Annotations") == FALSE) { dir.create("Annotations") }
  if (dir.exists("Annotations") == TRUE) {
    if (overwrite == TRUE){
      system(paste0("rm -r ", "Annotations"))
      dir.create("Annotations")
    } else { stop("overwrite is false and directory exists. Exiting.") }
  }#end dir exists

  if (dir.exists("Annotations/sample-contigs") == FALSE) { dir.create("Annotations/sample-contigs") }
  if (dir.exists("Annotations/sample-contigs") == TRUE) {
    if (overwrite == TRUE){
      system(paste0("rm -r ", "Annotations/sample-contigs"))
      dir.create("Annotations/sample-contigs")
    } else { stop("overwrite is false and directory exists. Exiting.") }
  }#end dir exists

  if (dir.exists("Annotations/sample-markers") == FALSE) { dir.create("Annotations/sample-markers") }
  if (dir.exists("Annotations/sample-markers") == TRUE) {
    if (overwrite == TRUE){
      system(paste0("rm -r ", "Annotations/sample-markers"))
      dir.create("Annotations/sample-markers")
    } else { stop("overwrite is false and directory exists. Exiting.") }
  }#end dir exists

  if (dir.exists("Annotations/sample-summary") == FALSE) { dir.create("Annotations/sample-summary") }
  if (dir.exists("Annotations/sample-summary") == TRUE) {
    if (overwrite == TRUE){
      system(paste0("rm -r ", "Annotations/sample-summary"))
      dir.create("Annotations/sample-summary")
    } else { stop("overwrite is false and directory exists. Exiting.") }
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

  for (i in 1:length(spp.samples)){

    #Load in the data
    contigs = Biostrings::readDNAStringSet(paste0(contig.folder, "/", spp.samples[i], ".fa"))   # loads up fasta file

    #Matches samples to loci
    system(paste0(blast.path, "blastn -task dc-megablast -db ", reference.name, "/ref-blast-db",
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

    ### Here reorient contigs, merge the same one?


    ##############################################################
    # Part A: Fixes up the match database
    ##############################################################

    #Gets rid of very poor matches
    filt.data = match.data[match.data$matches > 8,]
    filt.data = filt.data[filt.data$evalue <= 0.05,]
    filt.data = filt.data[grep("tRNA", filt.data$tName, invert = T),]

    if (nrow(filt.data) == 0){
      print("No matching mitochondrial genes were found.")
      system(paste0("rm ", spp.samples[i], "_match.txt"))
      next
    }#end if

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
#
#       if (nrow(complete.data) == 0) {
#         sub.data
#
#       }

      #Fragmented across contigs
      cut.contigs = Biostrings::DNAStringSet()
      for (k in 1:nrow(sub.data)){

        temp.contig = contigs[names(contigs) %in% sub.data$qName[k]]
        new.contig = Biostrings::subseq(temp.contig, start = sub.data$qStart[k], end = sub.data$qEnd[k])
        cut.contigs = append(cut.contigs, new.contig)
      }

      cap3.contigs = runCap3(cut.contigs,
                             cap3.path = cap3.path)
      #Saves if there is one contig
      if (length(cap3.contigs) == 1 && cap3.contigs != 0){
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

          if (nrow(sub.data)-1 <= k){ break }
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
          if (nrow(sub.data)-1 == k){ break }
          index = index + 1
        }#end k

        paste.contigs = Biostrings::DNAStringSet()
        for (k in 1:(nrow(sub.data)-1) ){

          if (k >= nrow(sub.data)-1){ next }
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

    system(paste0("rm ", spp.samples[i], "_match.txt"))

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
    #            nbchar = 1000000, as.string = T)
    #
    #
    #

    ######## Final step tRNA annotation
    ##############################################
    ##############################################

    ### Annotate tRNAs
    rna.data = tRNAscan(contigs = good.data,
                        tRNAscan.path = tRNAscan.path,
                        organism.type = "vertebrate",
                        quiet = TRUE)

    if (is.null(rna.data) == FALSE){ if (nrow(rna.data) == 0){ rna.data = NULL }}

    if (is.null(rna.data) == FALSE){

      # #Remove duplicate RNAs
      # if (length(unique(good.match$qName)) == 1){
      #   rna.data = rna.data[rna.data$contig %in% unique(good.match$qName),]
      # }

      rna.contigs = Biostrings::DNAStringSet()
      for (j in 1:nrow(rna.data)){

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

    ### Writes the summary
    write.csv(final.sample, file = paste0("Annotations/sample-summary/", spp.samples[i], "_sample-summary.csv"))

    #Writes the full mitochondrial genome file
    names(good.data) = paste0(spp.samples[i], "_|_", names(good.data))
    good.data = good.data[order(names(good.data))]
    write.loci = as.list(as.character(good.data))
    PhyloProcessR::writeFasta(sequences = write.loci, names = names(write.loci),
               paste0("Annotations/sample-markers/", spp.samples[i], "_sampleMarkers.fa"), nbchar = 1000000, as.string = T)

    #Writes the full mitochondrial genome file
    #good.contigs = contigs[names(contigs) %in% unique(refine.match$qName)]
    #names(good.contigs) = paste0(names(good.contigs))
    write.loci = as.list(as.character(good.data))
    PhyloProcessR::writeFasta(sequences = write.loci, names = names(write.loci),
               paste0("Annotations/sample-contigs/", spp.samples[i], "_sampleContigs.fa"), nbchar = 1000000, as.string = T)

    print(paste0("Finished annotation for ", spp.samples[i]))

  }#End i loop

  print("Finished annotation for all samples!")

}#end function
