#' @title alignMitogenomes
#'
#' @description Function for batch trimming a folder of alignments, with the various trimming functions available to select from
#'
#' @param alignment.folder path to a folder of sequence alignments in phylip format.
#'
#' @param genbank.file available input alignment formats: fasta or phylip
#'
#' @param draft.contigs contigs are added into existing alignment if algorithm is "add"
#'
#' @param output.dir available output formats: phylip
#'
#' @param dataset.name remove samples too divergent from consensus, values 0-1 for proportion similar sites
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

alignMitogenomes = function(alignment.folder = NULL,
                            genbank.file = NULL,
                            draft.contigs = "draftContigs",
                            output.dir = "Genomes",
                            dataset.name = "untrimmed",
                            overwrite = FALSE) {

  #Debug
  # alignment.folder = "Alignments/untrimmed-alignments"
  # genbank.file = "Crocidura.gb"
  # draft.contigs = "draftContigs"
  # output.dir = "Genomes"
  # dataset.name = "untrimmed"
  # overwrite = TRUE

  if (dir.exists(output.dir) == FALSE) { dir.create(output.dir) }
  if (dir.exists(output.dir) == TRUE) {
    if (overwrite == TRUE){
      system(paste0("rm -r ", output.dir))
      dir.create(output.dir)
    }
  }#end dir exists

  #Gets the samples
  locus.names = list.files(alignment.folder)
  locus.names = gsub(".phy$", "", locus.names)
  sample.names = list.files(draft.contigs)
  sample.names = gsub(".fa$", "", sample.names)

  #Create new direcotires
  dir.create(paste0(output.dir, "/alignments"))
  dir.create(paste0(output.dir, "/sample-genomes"))

  #Makes reference for blasting
  makeReference(genbank.file = genbank.file,
                overwrite = TRUE,
                rep.origin = FALSE)

  #Sets up the loci to align
  ref.data = Rsamtools::scanFa(Rsamtools::FaFile("Mito-Reference/refMarkers.fa"))

  #Header data for features and whatnot
  header.data = c("Sample",  locus.names)
  #Gets the new tree fiels it made
  collect.data.bp = data.table(matrix(as.numeric(0), nrow = length(sample.names), ncol = length(header.data)))
  setnames(collect.data.bp, header.data)
  collect.data.bp[, Sample:=as.character(sample.names)]

  #sets up data collection for proportion
  collect.data.pr = data.table(matrix(as.numeric(0), nrow = length(sample.names), ncol = length(header.data)))
  setnames(collect.data.pr, header.data)
  collect.data.pr[, Sample:=as.character(sample.names)]

  #sets up data collection for proportion
  feat.headers = c("Marker", "Start", "End")
  feature.data = data.table(matrix(as.numeric(0), nrow = length(locus.names), ncol = length(feat.headers)))
  setnames(feature.data, feat.headers)
  feature.data[, Marker:=as.character(Marker)]

  draft.genome = c()
  for (i in 1:length(locus.names)){

    # Read in untrimmed fasta and use with reference
    align = Biostrings::DNAStringSet(Biostrings::readAAMultipleAlignment(file = paste0(alignment.folder, "/", locus.names[i], ".phy"), format = "phylip"))

    #Gathers reference
    ref.temp = ref.data[names(ref.data) %in% locus.names[i]]
    align.data = append(align, ref.temp)
    names(align.data)[length(align.data)] = "Reference"

    new.align = runMafft(sequence.data = align.data,
                         add.contigs = ref.seq,
                         adjust.direction = T,
                         threads = 1)

    #Use taxa remove
    mat.align = strsplit(as.character(new.align), "")
    mat.align = lapply(mat.align, tolower)

    # Fill end gaps with Ns
    save.matrix = c()
    for (j in 1:length(mat.align)){

      if (names(mat.align)[j] == "Reference"){ next }

      tar.seq = mat.align[[j]]
      ref.seq = mat.align[[grep("Reference", names(mat.align))]]
      #replace real gap with N
      for (k in 1:length(ref.seq)){
        if (ref.seq[k] != "-" && tar.seq[k] == "-"){ tar.seq[k] = "n" }
      }#end k loop

      save.seq = matrix(tar.seq, nrow = 1)
      rownames(save.seq) = names(mat.align)[j]
      save.matrix = rbind(save.matrix, save.seq)

    } #end j loop

    #Fill in missing samples
    miss.samples = sample.names[!sample.names %in% rownames(save.matrix)]
    if (length(miss.samples) != 0){
      #Adds Ns in for each missing sample
      for (j in 1:length(miss.samples)){
        save.seq = matrix(as.character("n"), nrow = 1, ncol = ncol(save.matrix))
        rownames(save.seq) = miss.samples[j]
        save.matrix = rbind(save.matrix, save.seq)
      }#end j loop
    } #end if

    #Sort and concatenate alignments
    order.matrix = save.matrix[order(rownames(save.matrix)),]

    if (is.null(ncol(draft.genome)) == T){
      start = 1
      end = ncol(order.matrix)
    } else {
      start = ncol(draft.genome) + 1
      end = start + ncol(order.matrix) - 1
    }#end else

    draft.genome = cbind(draft.genome, order.matrix)

    #Saves location data
    set(feature.data, i = as.integer(i), j = match("Marker", feat.headers), value =  locus.names[i])
    set(feature.data, i = as.integer(i), j = match("Start", feat.headers), value =  start)
    set(feature.data, i = as.integer(i), j = match("End", feat.headers), value =  end)

    #Collects stats data
    #per complete
    bp.data = apply(order.matrix, MARGIN = 1, FUN = function (x) length(x[x != "n"]))

    set(collect.data.bp, i = match(collect.data.bp$Sample, names(bp.data)) ,
        j = match(locus.names[i], header.data), value =  bp.data)

    set(collect.data.pr, i = match(collect.data.pr$Sample, names(bp.data)) ,
        j = match(locus.names[i], header.data), value =  round(bp.data/(max(bp.data)), 3))

  }#end i loop

  #Save files
  write.csv(collect.data.bp, file = paste0("logs/", dataset.name, "_mito-alignment_bp-count.csv"),  row.names = F)
  write.csv(collect.data.pr, file = paste0("logs/", dataset.name, "_mito-alignment_bp-prop.csv"),  row.names = F)

  write.genome = as.matrix(ape::as.DNAbin(draft.genome) )
  writePhylip(write.genome, file= paste0(output.dir, "/alignments/", dataset.name, "_mitogenome_alignment.phy"), interleave = F)
  write.table(feature.data, file = paste0(output.dir, "/alignments/", dataset.name, "_alignment_feature_table.txt"),  row.names = F, quote = F)

}#end function

