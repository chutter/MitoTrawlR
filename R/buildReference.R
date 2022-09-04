#' @title buildReference
#'
#' @description Function for running the program spades to assemble short read sequencing data
#'
#' @param genbank.file path to a folder of sequence alignments in phylip format.
#'
#' @param overwrite path to a folder of sequence alignments in phylip format.
#'
#' @param rep.origin path to a folder of sequence alignments in phylip format.
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

#Add something to start reference at point of replication

buildReference = function(reference.fasta = NULL,
                          annotation.file = NULL,
                          annotation.type = c("genbank", "table", "gff"),
                          reference.name = "reference",
                          overwrite = TRUE,
                          rep.origin = FALSE) {

  #Debug
  #reference.fasta = "/Users/chutter/Dropbox/Research/1_Main-Projects/1_Collaborative-Projects/Microhylidae_SeqCap/New_Work_2022/MitoCap/Mantella.fa"
  #annotation.file = "/Users/chutter/Dropbox/Research/1_Main-Projects/1_Collaborative-Projects/Microhylidae_SeqCap/New_Work_2022/MitoCap/Mantella.gff"
  # annotation.file = "/Volumes/Rodents/Murinae/Mitochondrial_genomes/Mus_mus.gb"
  # annotation.file = "/Volumes/Rodents/Murinae/Mitochondrial_genomes/Mus_mus.gff"
  # annotation.type = "gff"
   #reference.name = "reference"
   #overwrite = TRUE
   #rep.origin = FALSE

  if (annotation.type != "genbank" && is.null(reference.fasta) == TRUE){
    stop("Please provide a reference fasta file or a single genbank file.")
  }

  if (length(annotation.type) != 1){ stop("Only one annotation type can be selected.") }

  #Checks for directory existing
  if (dir.exists(reference.name) == FALSE) { dir.create(reference.name) }
  if (dir.exists(reference.name) == TRUE) {
    if (overwrite == TRUE){
      system(paste0("rm -r ", reference.name))
      dir.create(reference.name)
    } else { return("overwrite = FALSE and Mito-Reference directory exists") }
  }#end dir exists


  ######################################
  ### PARSE genbank file
  ######################################
  if (annotation.type == "genbank"){

    gb.data = genbankr::readGenBank(file = annotation.file, partial = T)
    ref.genome = gb.data@sequence
    all.data = data.frame(gb.data@genes)

    if (nrow(all.data) == 0){
      all.data = data.frame(gb.data@cds)
    }

    all.data$transcript_id = NULL
    all.data$anticodon = NULL
    all.data$db_xref = NULL
    all.data$name = all.data$gene
    all.data$gene = NULL
    all.data$gene_id = NULL
    all.data$note = NULL
    all.data$codon_start = NULL

    #RNA data stuff
    rna.data = data.frame(gb.data@other_features)
    rna.data$note = NULL
    rna.data$anticodon = NULL
    rna.data$codon_recognized = NULL
    rna.data$product = gsub(" ", "-", rna.data$product)
    rna.data$name = rna.data$product
    rna.data$product = NULL
    rna.data$db_xref = NULL

    if (rep.origin == TRUE){
      rna.data[rna.data$type == "rep_origin",]$name = "rep_origin"
      add.all = cbind(rna.data[rna.data$type == "rep_origin",], gene = NA, gene_id = NA)
      add.all$product = NULL
      all.data = rbind(add.all, all.data)
    } else {
      rna.data = rna.data[rna.data$type != "rep_origin",]
    }

    #Combined tgother
    all.data = rbind(all.data, rna.data)
    all.data = all.data[order(all.data$start),]
    all.data$loctype = NULL
    all.data$width = NULL

  }#end genbank file

  ######################################
  ### PARSE GFF file
  ######################################
  if (annotation.type == "gff"){
    ref.genome = Biostrings::readDNAStringSet(file = reference.fasta)
    names(ref.genome) = reference.name
    #Pulls in GFF
    gtf.data = data.table::data.table(read.delim(annotation.file, sep = "\t", comment.char = "#"))

    header.names = c("seqnames", "method", "type", "start", "end", "un1", "strand", "un2", "info")
    data.table::setnames(gtf.data, header.names)
    gtf.data = gtf.data[gtf.data$type != "exon",]
    gtf.data = gtf.data[gtf.data$type != "gene",]
    gtf.data = gtf.data[gtf.data$type != "region",]
    gtf.data[, un1 := NULL]
    gtf.data[, un2 := NULL]
    gtf.data[, method := NULL]
    gtf.data[, name := gsub(".*;product=", "", gtf.data$info) ]
    gtf.data[, name := gsub(";.*", "", gtf.data$name) ]
    gtf.data[, info := NULL]

    #Fixes D-loop
    gtf.data[gtf.data$type == "D_Loop",]$name = "D_loop"
    gtf.data[gtf.data$type == "D_loop",]$name = "D_loop"
    gtf.data[gtf.data$type == "D-Loop",]$name = "D_loop"
    gtf.data[gtf.data$type == "D-loop",]$name = "D_loop"
    gtf.data[gtf.data$type == "d_Loop",]$name = "D_loop"
    gtf.data[gtf.data$type == "d_loop",]$name = "D_loop"
    gtf.data[gtf.data$type == "d-Loop",]$name = "D_loop"
    gtf.data[gtf.data$type == "d-loop",]$name = "D_loop"
    gtf.data[, name := gsub(" ", "-", gtf.data$name) ]
    all.data = gtf.data
  }#end GFF

  ######################################
  ### Read in simple table
  ######################################
  if (annotation.type == "table"){
    ref.genome = Biostrings::readDNAStringSet(file = reference.fasta)
    names(ref.genome) = reference.name
    #Pulls in GFF
    all.data = data.table::data.table(read.table(annotation.file, sep = "\t", header = T))
  }#END TABLE


  #Goes through each entry and saves
  save.seq = c()
  for (i in 1:nrow(all.data)){

    temp.data = data.frame(all.data[i,])

    #Gets ranges
    start = temp.data$start
    end = temp.data$end

    # #Gets type
    # type.data.e = exon.data[exon.data$name %in% temp.data$name,]
    # type.data.r = rna.data[rna.data$name %in% temp.data$name,]
    #
    # if (nrow(type.data.r) == 0 & nrow(type.data.e) == 0){
    #   type.data.e = exon.data[grep(temp.data$name, exon.data$name),]
    #   type.data.r = rna.data[grep(temp.data$name, rna.data$name),]
    # }#end if
    #
    # if (nrow(type.data.r) == 1 & nrow(type.data.e) == 1){stop("found both rna and exon annotation in genbank file.") }
    # if (nrow(type.data.r) > 1 | nrow(type.data.e) > 1){stop("found both too many annotations in genbank file.") }
    #
    # #If it matches rna
    # if (nrow(type.data.r) == 1){
    #   #Obtains sequence
    #   new.seq = Biostrings::subseq(ref.genome, start = start, end = end)
    #   if (as.character(type.data.r$strand) == "-"){ new.seq = Biostrings::reverseComplement(new.seq) }
    #   names(new.seq) = paste0(sprintf("%02d", i), "_", type.data.r$product)
    #   save.seq = append(save.seq, new.seq)
    # }#end
    #
    # #Exon data
    # if (nrow(type.data.e) == 1){
    #   #Checks strand for codon start
    #   if (as.character(type.data.e$strand) == "+"){
    #    # if (type.data.e$codon_start == 2){ start = start + 1 }
    #   #  if (type.data.e$codon_start == 3){ start = start + 2 }
    #   }#end plus
    #
    #   if (as.character(type.data.e$strand) == "-"){
    #     new.seq = Biostrings::reverseComplement(new.seq) }
    #
    #    # if (type.data.e$codon_start == 2){ end = end - 1 }
    #   #  if (type.data.e$codon_start == 3){ end = end - 2 }
    #   }#end plus

      #Obtains sequence
      new.seq = Biostrings::subseq(ref.genome, start = start, end = end)
      if (as.character(temp.data$strand) == "-"){ new.seq = Biostrings::reverseComplement(new.seq) }
      names(new.seq) = paste0(sprintf("%02d", i), "_", temp.data$type, "_", temp.data$name)

      if (nchar(names(new.seq)) >= 50){
        names(new.seq) = substr(names(new.seq), 1, 50)
      }

      save.seq = append(save.seq, new.seq)
  }#end i loop

  #Writes reference loci
  write.loci = as.list(as.character(save.seq))
  PhyloCap::writeFasta(sequences = write.loci, names = names(write.loci),
                       paste0(reference.name, "/refMarkers.fa"), nbchar = 1000000, as.string = T)

  #WRites reference mito genome
  write.loci = as.list(as.character(ref.genome))
  PhyloCap::writeFasta(sequences = write.loci, names = names(write.loci),
                       paste0(reference.name, "/refGenome.fa"), nbchar = 1000000, as.string = T)

  write.table(all.data, file = paste0(reference.name, "/referenceTable.txt", sep = "\t"), row.names = F, quote = F)

  print(paste0("Mitochondrial reference for ", reference.name, " created in Mito-Reference"))

}#end function
