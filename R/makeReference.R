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

  #genbank.file = "Crocidura.gb"
  #overwrite = TRUE

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
  gb.data = genbankr::readGenBank(file = genbank.file)
  ref.genome = gb.data@sequence
  all.data = data.frame(gb.data@genes)
  exon.data = data.frame(gb.data@exons)
  rna.data = data.frame(gb.data@other_features)
  rna.data$product = gsub(" ", "-", rna.data$product)

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

    #Gets type
    type.data.e = exon.data[exon.data$db_xref %in% temp.data$db_xref,]
    type.data.r = rna.data[rna.data$db_xref %in% temp.data$db_xref,]

    if (nrow(type.data.r) == 0 & nrow(type.data.e) == 0){
      type.data.e = exon.data[grep(temp.data$db_xref, exon.data$db_xref),]
      type.data.r = rna.data[grep(temp.data$db_xref, rna.data$db_xref),]
    }#end if

    if (nrow(type.data.r) == 1 & nrow(type.data.e) == 1){stop("found both rna and exon annotation in genbank file.") }
    if (nrow(type.data.r) > 1 | nrow(type.data.e) > 1){stop("found both too many annotations in genbank file.") }

    #If it matches rna
    if (nrow(type.data.r) == 1){
      #Obtains sequence
      new.seq = Biostrings::subseq(ref.genome, start = start, end = end)
      if (as.character(type.data.r$strand) == "-"){ new.seq = Biostrings::reverseComplement(new.seq) }
      names(new.seq) = paste0(sprintf("%02d", i), "_", type.data.r$product)
      save.seq = append(save.seq, new.seq)
    }#end

    #Exon data
    if (nrow(type.data.e) == 1){
      #Checks strand for codon start
      if (as.character(type.data.e$strand) == "+"){
        if (type.data.e$codon_start == 2){ start = start + 1 }
        if (type.data.e$codon_start == 3){ start = start + 2 }
      }#end plus

      if (as.character(type.data.e$strand) == "-"){
        if (type.data.e$codon_start == 2){ end = end - 1 }
        if (type.data.e$codon_start == 3){ end = end - 2 }
      }#end plus

      #Obtains sequence
      new.seq = Biostrings::subseq(ref.genome, start = start, end = end)
      if (as.character(type.data.e$strand) == "-"){ new.seq = Biostrings::reverseComplement(new.seq) }
      names(new.seq) = paste0(sprintf("%02d", i), "_", type.data.e$type, "-", type.data.e$gene)
      save.seq = append(save.seq, new.seq)
    }#end type.data.e
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
