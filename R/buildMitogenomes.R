#' @title buildMitogenomes
#'
#' @description Function for batch trimming a folder of alignments, with the various trimming functions available to select from
#'
#' @param annotation.dir path to a folder of sequence alignments in phylip format.
#'
#' @param alignment.folder path to a folder of sequence alignments in phylip format.
#'
#' @param genome.alignment available input alignment formats: fasta or phylip
#'
#' @param genome.dir contigs are added into existing alignment if algorithm is "add"
#'
#' @param output.dir available output formats: phylip
#'
#' @param overwrite TRUE to supress mafft screen output
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

buildMitogenomes = function(annotation.dir = "Annotations",
                            alignment.folder = "Alignments/untrimmed-alignments",
                            genome.alignment = "Genomes/alignments/untrimmed_mitogenome_alignment.phy",
                            genome.dir = "Genomes",
                            output.dir = NULL,
                            overwrite = FALSE) {

  #Debug
  # alignment.folder = "Alignments/untrimmed-alignments"
  # genome.dir = "Genomes"
  # annotation.dir = "Annotations"
  # genome.alignment = "Genomes/alignments/untrimmed_mitogenome_alignment.phy"
  # output.dir = "untrimmed-finished"
  # overwrite = FALSE

  #Prepares output directoires
  output.dir = paste0(genome.dir, "/", output.dir)
  if (dir.exists(output.dir) == FALSE) { dir.create(output.dir) }
  if (dir.exists(output.dir) == TRUE) {
    if (overwrite == TRUE){
      system(paste0("rm -r ", output.dir))
      dir.create(output.dir)
    }
  }#end dir exists

  #Create new direcotires
  dir.create(paste0(genome.dir, "/reference-order"))
  dir.create(paste0(genome.dir, "/sample-markers"))
  dir.create(paste0(genome.dir, "/final-genomes"))

  #Makes reference for blasting
  #makeReference(genbank.file = genbank.file,
  #             overwrite = TRUE,
  #             rep.origin = FALSE)

  #Read in reference and alignment
  ref.data = Rsamtools::scanFa(Rsamtools::FaFile("Mito-Reference/refMarkers.fa"))
  gen.align = Biostrings::DNAStringSet(Biostrings::readAAMultipleAlignment(file = genome.alignment, format = "phylip"))
  sample.names = names(gen.align)
  locus.names = names(ref.data)
  first.marker = names(ref.data)[1]
  last.marker = names(ref.data)[length(ref.data)]

  #Header data for features and whatnot
  header.data = c("Sample", "Circular", "Fragments", "bp_Complete", "Total_Length",
                  "bp_Count", "N_Count", "Marker_Complete", "Marker_Count",
                  "Missing_Markers", "Missing_BP")
  #Gets the new tree fiels it made
  collect.data = data.table(matrix(as.numeric(0), nrow = length(sample.names), ncol = length(header.data)))
  setnames(collect.data, header.data)
  collect.data[, Sample:=as.character(sample.names)]

  #Loops through each sample to assemble genomes
  for (i in 1:length(sample.names)){
    # Read in untrimmed fasta and use with reference
    #n.count = sum(stringr::str_count(as.character(sample.align), "N"))
    # if (n.count >= (length(sample.align) * 0.90)){ next }
    #contigs = Rsamtools::scanFa(Rsamtools::FaFile(paste0(annotation.dir, "/sample-markers/", sample.names[i], "_sampleMarkers.fa")))   # loads up fasta file
    contigs = Rsamtools::scanFa(Rsamtools::FaFile(paste0(annotation.dir, "/sample-contigs/", sample.names[i], "_sampleContigs.fa")))   # loads up fasta file
    n.fragments = length(contigs)

    ##################################
    #Saves reference order sequence
    ##################################
    sample.align = gen.align[names(gen.align) %in% sample.names[i]]
    stand.order = as.list(as.character(sample.align))
    #Saves to folder the standard order one already made
    writeFasta(sequences = stand.order, names = names(stand.order),
               paste0("Genomes/reference-order/", sample.names[i], "_referenceOrder.fa"),
               nbchar = 1000000, as.string = T)

    ##################################
    #Saves original order sequence
    ##################################
    sample.data = Biostrings::DNAStringSet()
    for (j in 1:length(locus.names)){

      align = Biostrings::DNAStringSet(Biostrings::readAAMultipleAlignment(file = paste0(alignment.folder, "/", locus.names[j], ".phy"), format = "phylip"))

      new.align = align[names(align) == sample.names[i]]

      if (length(new.align) == 0){ next }

      names(new.align) = locus.names[j]
      sample.data = append(sample.data, new.align)

    } #end j loop

    #Adds reference locus
    final.save = as.list(as.character(sample.data))
    final.save = lapply(final.save, function (x) gsub("-", "", x) )
    #delete.n = lapply(final.save, function (x) gsub("N", "", x))
    delete.sample = names((final.save[final.save == ""]))
    final.write = final.save[!names(final.save) %in% delete.sample]

    #Saves to folder the standard order one already made
    writeFasta(sequences = final.write, names = names(final.write),
               paste0(genome.dir, "/sample-markers/", sample.names[i], "_Markers.fa"),
               nbchar = 1000000, as.string = T)

    #Make blast database for the probe loci
    #headers
    headers = c("qName", "tName", "pident", "matches", "misMatches", "gapopen",
                "qStart", "qEnd", "tStart", "tEnd", "evalue", "bitscore", "qLen", "tLen", "gaps")

    system(paste0("makeblastdb -in ", genome.dir, "/sample-markers/", sample.names[i], "_Markers.fa -parse_seqids -dbtype nucl ",
                  " -out ", genome.dir, "/sample_blast_db"), ignore.stdout = T, ignore.stderr = T)

    #Matches samples to loci
    system(paste0("blastn -task dc-megablast -db ", genome.dir, "/sample_blast_db",
                  " -query ", annotation.dir, "/sample-contigs/", sample.names[i], "_sampleContigs.fa -out ", sample.names[i], "_match.txt",
                  " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen gaps\" ",
                  " -num_threads 1"))

    #Need to load in transcriptome for each species and take the matching transcripts to the database
    match.data = read.table(paste0(sample.names[i], "_match.txt"), sep = "\t", header = F, stringsAsFactors = FALSE)
    colnames(match.data) = headers

    if (nrow(match.data) == 0){
      print("No matching mitochondrial genes were found.")
      next
    }#end if

    #Check through each marker
    ##############################################################
    # Part A: Fixes up the match database
    ##############################################################

    #Gets rid of very poor matches
    filt.data = match.data[match.data$matches > 10,]
    filt.data = filt.data[filt.data$evalue <= 0.05,]
    filt.data = filt.data[filt.data$pident >= 95,]
    filt.data = filt.data[order(filt.data$qStart),]

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

    #Checks for circularity
    # circular = F
    # n.fragments = c()
    # if (length(unique(filt.data$qName)) == 1){
    #   filt.data = filt.data[order(filt.data$qStart, decreasing = F),]
    #
    #   if (filt.data$tName[1] == filt.data$tName[nrow(filt.data)]) {  circular = T } else {circular = F }
    #   n.fragments = 1
    #
    # } else {  n.fragments = length(unique(filt.data$qName)) }

    #CIRCULAR CHECK

    #If there is only fragment, reorganize to get gene order
    if (length(unique(filt.data$qName)) == 1){

      first.marker = names(final.write)[1]
      #if (nrow(filt.data[filt.data$tName == first.marker,]) >= 2){ stop("whoa") }

      if (unique(filt.data[filt.data$tName == first.marker,]$qDir) == "+"){

        #start.marker = filt.data[filt.data$tName == first.marker,]
        #end.marker = filt.data[filt.data$tName == last.marker,]

        start.pos = filt.data[filt.data$tName == first.marker,]$qStart
        temp.contig = contigs[names(contigs) == filt.data$qName[1]]
        first.part = Biostrings::subseq(temp.contig, start = start.pos[1], end = Biostrings::width(temp.contig))
        second.part = Biostrings::subseq(temp.contig, start = 1, end = start.pos[1]-1)

        #Combine together
        combined.contig = Biostrings::DNAStringSet(paste0(as.character(first.part),
                                                          as.character(second.part)))

      }#end if

      if (unique(filt.data[filt.data$tName == first.marker,]$qDir) == "-"){

        #start.marker = filt.data[filt.data$tName == last.marker,]
        #end.marker = filt.data[filt.data$tName == first.marker,]

        start.pos = filt.data[filt.data$tName == first.marker,]$qEnd + 1
        temp.contig = contigs[names(contigs) == filt.data$qName[1]]
        first.part = Biostrings::subseq(temp.contig, start = start.pos[1], end = Biostrings::width(temp.contig))
        second.part = Biostrings::subseq(temp.contig, start = 1, end = start.pos[1] - 1)

        #Combine together
        combined.contig = Biostrings::DNAStringSet(paste0(as.character(first.part),
                                                          as.character(second.part)))
        #Reverse compliment
        #Reverses alignment back to correction orientation
        combined.contig = Biostrings::reverseComplement(combined.contig)

      }#end if

      #Adds reference locus
      save.contig = as.list(as.character(combined.contig))
      names(save.contig) = paste0(sample.names[i], "_sequence-", rep(1:length(save.contig)))
      #Saves to folder the standard order one already made
      writeFasta(sequences = save.contig, names = names(save.contig),
                 paste0(genome.dir, "/final-genomes/", sample.names[i], "_Complete.fa"),
                 nbchar = 1000000, as.string = T)
    }#end nfragment if here

    if (length(unique(filt.data$qName)) != 1){

      frag.contig = Biostrings::DNAStringSet()
      f.names = unique(filt.data$qName)
      for (k in 1:length(f.names)){
        temp.filt = filt.data[filt.data$qName == f.names[k],]

        if (length(unique(temp.filt$qDir)) !=1){
          temp.contig = contigs[names(contigs) == f.names[k]]
          frag.contig = append(frag.contig, temp.contig)
        } else {

          if (unique(temp.filt$qDir) == "+"){

            temp.contig = contigs[names(contigs) == f.names[k]]
            frag.contig = append(frag.contig, temp.contig)
          }#end if

          if (unique(temp.filt$qDir) == "-"){
            temp.contig = contigs[names(contigs) == f.names[k]]
            #Reverses alignment back to correction orientation
            frag.contig = append(frag.contig, Biostrings::reverseComplement(temp.contig))
          }#end if
        }#end k loop
      }#end out if for qdir

      save.contig = as.list(as.character(frag.contig))
      #Saves to folder the standard order one already made
      names(save.contig) = paste0(sample.names[i], "_sequence-", rep(1:length(save.contig)))
      writeFasta(sequences = save.contig, names = names(save.contig),
                 paste0(genome.dir, "/final-genomes/", sample.names[i], ".fa"),
                 nbchar = 1000000, as.string = T)
    }#end if


    ### Summary data collection
    #Saves location data
    set(collect.data, i = as.integer(i), j = match("Sample", header.data), value =  sample.names[i])

    no.gap = gsub("-", "", stand.order)
    total.nogap = nchar(no.gap)
    total.bp = nchar(gsub("N","", no.gap))
    total.n = total.nogap - total.bp

    #set(collect.data, i = as.integer(i), j = match("Circular", header.data), value =  circular)
    set(collect.data, i = as.integer(i), j = match("Fragments", header.data), value =  n.fragments)
    set(collect.data, i = as.integer(i), j = match("bp_Complete", header.data), value =  round(total.bp/total.nogap, 3))
    set(collect.data, i = as.integer(i), j = match("Total_Length", header.data), value =  total.nogap)
    set(collect.data, i = as.integer(i), j = match("bp_Count", header.data), value =  total.bp)
    set(collect.data, i = as.integer(i), j = match("N_Count", header.data), value =  total.n)

    total.markers = length(final.save)
    total.count = length(final.write)

    set(collect.data, i = as.integer(i), j = match("Marker_Complete", header.data), value =  round(total.count/total.markers, 3) )
    set(collect.data, i = as.integer(i), j = match("Marker_Count", header.data), value =  total.count)
    set(collect.data, i = as.integer(i), j = match("Missing_Markers", header.data), value =  total.markers - total.count)
    set(collect.data, i = as.integer(i), j = match("Missing_BP", header.data), value =  round(total.n/total.nogap, 3))

    system(paste0("rm ", sample.names[i], "_match.txt"))

  }#i loop

  write.csv(collect.data, file = paste0(genome.dir, "/sample_mito_genome_summary.csv"),  row.names = F)

}#end function
