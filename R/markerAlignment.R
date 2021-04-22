#' @title markerAlignment
#'
#' @description Function for running the program spades to assemble short read sequencing data
#'
#' @param input.folder path to a folder of sequence alignments in phylip format.
#'
#' @param genbank.file contigs are added into existing alignment if algorithm is "add"
#'
#' @param out.dir contigs are added into existing alignment if algorithm is "add"
#'
#' @param min.taxa algorithm to use: "add" add sequences with "add.contigs"; "localpair" for local pair align. All others available
#'
#' @param min.prop.coverage TRUE to supress screen output
#'
#' @param threads TRUE to supress screen output
#'
#' @param overwrite TRUE to supress screen output
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

#Aligns all the different markers
markerAlignment = function(input.folder = NULL,
                           reference.name = NULL,
                           threads = 1,
                           overwrite = TRUE){

  #Debug
  # input.folder = "Annotations/sample-fastas"
  # genbank.file = "Crocidura.gb"
  # threads = 4
  # overwrite = TRUE

  #Checks for overwrite
  if (dir.exists("Alignments") == FALSE) { dir.create("Alignments") }
  if (dir.exists("Alignments") == TRUE) {
    if (overwrite == TRUE){
      system(paste0("rm -r ", "Alignments"))
      dir.create("Alignments")
    } else { stop("overwrite is false and directory exists. Exiting.") }
  }#end dir exists

  #Checks and creates directories based on overwrite parameter
  if (dir.exists("Alignments/untrimmed-alignments") == FALSE) { dir.create("Alignments/untrimmed-alignments") }
  if (dir.exists("Alignments/untrimmed-alignments") == TRUE) {
    if (overwrite == TRUE){
      system(paste0("rm -r ", "Alignments/untrimmed-alignments"))
      dir.create("Alignments/untrimmed-alignments")
    } else { stop("overwrite is false and directory exists. Exiting.") }
  }#end dir exists

  if (dir.exists("Alignments/unaligned-markers") == FALSE) { dir.create("Alignments/unaligned-markers") }
  if (dir.exists("Alignments/unaligned-markers") == TRUE) {
    if (overwrite == TRUE){
      system(paste0("rm -r ", "Alignments/unaligned-markers"))
      dir.create("Alignments/unaligned-markers")
    } else { stop("overwrite is false and directory exists. Exiting.") }
  }#end dir exists

  #Sets up the loci to align
  ref.data = Rsamtools::scanFa(Rsamtools::FaFile(paste0(reference.name, "/refMarkers.fa")))

  #Gets the samples
  spp.samples = list.files(input.folder)
  spp.samples = spp.samples[spp.samples != ""]

  #Gets all species data
  final.contigs = Biostrings::DNAStringSet()
  for (i in 1:length(spp.samples)){
    #Loads file
    spp.data =Rsamtools::scanFa(Rsamtools::FaFile(paste0(input.folder, "/", spp.samples[i])))   # loads up fasta file
    #saves
    final.contigs = append(final.contigs, spp.data)
  }#end j loop

  #Aligns each potential locus
  for (i in 1:length(ref.data)){
    ##############
    #STEP 1: Sets up for alignment
    ##############
    #Checks for a minimum length
    sample.markers = final.contigs[grep(paste0(names(ref.data)[i], "$"), names(final.contigs))]

    #Checks for minimum taxa number
    if (length(names(sample.markers)) <= 2){
      print(paste0(names(ref.data)[i], " had too few taxa."))
      next
    }

    #Adds reference locus
    align.data = append(sample.markers, ref.data[i])
    names(align.data)[length(align.data)] = "Reference"
    final.save = as.list(as.character(align.data))

    #Saves to folder to run with mafft
    writeFasta(sequences = final.save, names = names(final.save),
               paste0("Alignments/unaligned-markers/", names(ref.data)[i], ".fa"), nbchar = 1000000, as.string = T)

    ##############
    #STEP 3: Runs MAFFT to align
    ##############
    mafft.cmd = "mafft"
    # if (names(ref.data)[i] == "12S_rRNA" || names(ref.data)[i] == "16S_rRNA"){
    #   if (secondary.structure == TRUE){ mafft.cmd<-"mafft-qinsi" } else { mafft.cmd<-"mafft" }
    # }
    #
    #Runs the mafft command
    system(paste0(mafft.cmd, " --localpair --maxiterate 1000 --adjustdirection --quiet --op 3 --ep 0.123",
                  " --thread ", threads, " ", "Alignments/unaligned-markers/", names(ref.data)[i], ".fa",
                  " > ", "Alignments/unaligned-markers/", names(ref.data)[i], "_align.fa"))

    alignment = Rsamtools::scanFa(Rsamtools::FaFile(paste0("Alignments/unaligned-markers/", names(ref.data)[i], "_align.fa")))   # loads up fasta file

    #Reverses alignment back to correction orientation
    reversed = names(alignment)[grep(pattern = "_R_", names(alignment))]
    if (length(reversed[grep(pattern = "Reference", reversed)]) == 1){ alignment = Biostrings::reverseComplement(alignment) }

    #Renames sequences to get rid of _R_
    names(alignment) = gsub(pattern = "_R_", replacement = "", x = names(alignment))
    names(alignment) = gsub(pattern = "_\\|_.*", replacement = "", x = names(alignment))

    new.align = strsplit(as.character(alignment), "")
    mat.align = lapply(new.align, tolower)
    m.align = as.matrix(ape::as.DNAbin(mat.align))

    #Trims to reference
    ref.aligned = as.character(alignment['Reference'])
    not.gaps = stringr::str_locate_all(ref.aligned, pattern = "[^-]")[[1]][,1]
    ref.start = min(not.gaps)
    ref.finish = max(not.gaps)
    trim.align = Biostrings::subseq(alignment, ref.start, ref.finish)
    trim.align = trim.align[names(trim.align) != "Reference"]

    ### Add in edge Ns
    sample.n = names(trim.align)
    align.list = vector(mode = "list")
    for (x in 1:length(sample.n)){
      ref.aligned = as.character(trim.align[sample.n[x]])
      gaps = stringr::str_locate_all(ref.aligned, pattern = "[-]")[[1]][,1]
      split.align = unlist(strsplit(as.character(ref.aligned), ""))

      #Checks for start gaps to replace with N
      for (y in 1:length(split.align)){
        if (split.align[y] != "-"){ break }
        if (split.align[y] == "-"){ split.align[y] = "N" }
      }#end forward loop

      #Checks for end gaps to replace with N
      for (y in length(split.align):1){
        if (split.align[y] != "-"){ break }
        if (split.align[y] == "-"){ split.align[y] = "N" }
      }#end forward loop

      #Saves in a list
      align.list[[x]] = split.align
      names(align.list)[x] = sample.n[x]
    }#end x

    #readies for saving
    aligned.set = as.matrix(ape::as.DNAbin(align.list) )
    writePhylip(aligned.set, file=paste0("Alignments/untrimmed-alignments/", names(ref.data)[i], ".phy"), interleave = F)
    system(paste0("rm Alignments/unaligned-markers/", names(ref.data)[i], "_align.fa"))
  }#end i loop

} #end function
