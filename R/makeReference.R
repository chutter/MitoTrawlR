#' @title makeReference
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

makeReference = function(genbank.file = NULL,
                         overwrite = TRUE,
                         rep.origin = FALSE) {

  #Debug

  #genbank.file = gb.file
  #overwrite = TRUE
  #rep.origin = FALSE

  #Checks for directory existing
  if (dir.exists("Mito-Reference") == FALSE) { dir.create("Mito-Reference") }
  if (dir.exists("Mito-Reference") == TRUE) {
    if (overwrite == TRUE){
      system(paste0("rm -r ", "Mito-Reference"))
      dir.create("Mito-Reference")
    } else { return("overwrite = FALSE and Mito-Reference directory exists") }
  }#end dir exists

  options(stringsAsFactors = FALSE)
  if (is.null(genbank.file) == TRUE){ stop("Please provide a genbank file.") }

  #Gebank setup
  gb.data = genbankr::readGenBank(file = genbank.file, partial = T)
  ref.genome = gb.data@sequence
  all.data = data.frame(gb.data@genes)
  all.data$anticodon = NULL
  all.data$db_xref = NULL
  all.data$name = all.data$gene
  all.data$gene = NULL
  all.data$gene_id = NULL

  #RNA data stuff
  rna.data = data.frame(gb.data@other_features)
  rna.data$note = NULL
  rna.data$anticodon = NULL
  rna.data$codon_recognized = NULL
  rna.data$product = gsub(" ", "-", rna.data$product)
  rna.data$name = rna.data$product
  rna.data$product = NULL

  #Combined tgother
  all.data = rbind(all.data, rna.data)
  all.data = all.data[order(all.data$start),]

  if (rep.origin == TRUE){
    rna.data[rna.data$type == "rep_origin",]$product = "rep_origin"
    rna.data[rna.data$type == "rep_origin",]$db_xref = "rep_origin"
    add.all = cbind(rna.data[rna.data$type == "rep_origin",], gene = NA, gene_id = NA)
    add.all$product = NULL
    all.data = rbind(add.all, all.data)
  }#end rep.origin

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
      names(new.seq) = paste0(sprintf("%02d", i), "_", temp.data$type, "-", temp.data$name)
      save.seq = append(save.seq, new.seq)
  }#end i loop

  #Writes reference loci
  write.loci = as.list(as.character(save.seq))
  writeFasta(sequences = write.loci, names = names(write.loci),
             "Mito-Reference/refMarkers.fa", nbchar = 1000000, as.string = T)

  #WRites reference mito genome
  write.loci = as.list(as.character(ref.genome))
  writeFasta(sequences = write.loci, names = names(write.loci),
             "Mito-Reference/refGenome.fa", nbchar = 1000000, as.string = T)

  print(paste0("Mitochondrial reference for ", genbank.file, " created in Mito-Reference"))

}#end function
