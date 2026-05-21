#' @title buildReference
#'
#' @description Builds a reference marker database from an annotated
#'   mitochondrial genome. Accepts annotation in GenBank (\code{.gb}), GFF3
#'   (\code{.gff}), or a simple tab-delimited table format. Extracts each
#'   annotated feature as an individual FASTA sequence, writes the per-marker
#'   FASTA (\code{refMarkers.fa}), the whole reference genome FASTA
#'   (\code{refGenome.fa}), and an annotation table
#'   (\code{referenceTable.txt}) to the \code{reference.name} directory.
#'
#' @param reference.fasta path to a FASTA file of the reference mitochondrial
#'   genome. Required when \code{annotation.type} is \code{"gff"} or
#'   \code{"table"}; ignored when \code{annotation.type = "genbank"} (the
#'   sequence is read from the GenBank file).
#'
#' @param annotation.file path to the annotation file. Should be a GenBank
#'   (\code{.gb}), GFF3 (\code{.gff}), or tab-delimited table file depending
#'   on \code{annotation.type}.
#'
#' @param annotation.type format of the annotation file: one of
#'   \code{"genbank"}, \code{"gff"}, or \code{"table"}.
#'
#' @param reference.name name for the output directory that will hold all
#'   reference files.
#'
#' @param overwrite logical; if TRUE, an existing reference directory is
#'   deleted and recreated.
#'
#' @param rep.origin logical; if TRUE, the replication origin feature
#'   (\code{rep_origin}) is included in the reference markers.
#'
#' @return Invisibly returns NULL. Writes \code{refMarkers.fa},
#'   \code{refGenome.fa}, and \code{referenceTable.txt} to
#'   \code{reference.name/}.
#'
#' @examples
#' \dontrun{
#' buildReference(
#'   annotation.file  = "Mus_musculus.gb",
#'   annotation.type  = "genbank",
#'   reference.name   = "reference",
#'   overwrite        = TRUE,
#'   rep.origin       = FALSE
#' )
#' }
#'
#' @export

#Add something to start reference at point of replication

buildReference = function(reference.fasta = NULL,
                          annotation.file = NULL,
                          annotation.type = c("genbank", "table", "gff"),
                          reference.name = "reference",
                          overwrite = FALSE,
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

  annotation.type = match.arg(annotation.type)

  if (annotation.type != "genbank" && is.null(reference.fasta)){
    stop("Please provide a reference fasta file when annotation.type is 'gff' or 'table'.")
  }

  #Checks for directory existing
  if (dir.exists(reference.name) == FALSE) {
    dir.create(reference.name)
  } else if (overwrite == TRUE) {
    unlink(reference.name, recursive = TRUE)
    dir.create(reference.name)
  } else {
    message("Reference directory '", reference.name, "' already exists and overwrite = FALSE. Skipping.")
    return(invisible(NULL))
  }


  ######################################
  ### PARSE genbank file
  ######################################
  if (annotation.type == "genbank"){

    gb.data = genbankr::readGenBank(file = annotation.file, partial = TRUE)
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
    names(ref.genome)[1] = reference.name
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

    # Filter replication origin unless requested
    if (!rep.origin) {
      gtf.data = gtf.data[gtf.data$type != "origin_of_replication",]
    }

    # For each row: prefer gene= (gives short names like ND1, COX1) over product=
    extract.gff.name = function(info.str) {
      gene.m = regmatches(info.str, regexpr("(?<=;gene=)[^;]+", info.str, perl = TRUE))
      if (length(gene.m) > 0 && nzchar(gene.m)) return(gene.m)
      prod.m = regmatches(info.str, regexpr("(?<=product=)[^;]+", info.str, perl = TRUE))
      if (length(prod.m) > 0 && nzchar(prod.m)) return(prod.m)
      return(NA_character_)
    }
    gtf.data[, name := sapply(info, extract.gff.name)]
    gtf.data[, info := NULL]

    #Fixes D-loop (case-insensitive, handles D_loop / D-loop / d-Loop etc.)
    gtf.data[grepl("^d[-_]loop$", gtf.data$type, ignore.case = TRUE), ]$name = "D_loop"
    gtf.data[, name := gsub(" ", "-", gtf.data$name) ]
    all.data = gtf.data
  }#end GFF

  ######################################
  ### Read in simple table
  ######################################
  if (annotation.type == "table"){
    ref.genome = Biostrings::readDNAStringSet(file = reference.fasta)
    names(ref.genome)[1] = reference.name
    #Pulls in GFF
    all.data = data.table::data.table(read.table(annotation.file, sep = "\t", header = TRUE))
  }#END TABLE


  #Goes through each entry and saves
  save.seq = c()
  for (i in seq_len(nrow(all.data))){

    temp.data = data.frame(all.data[i,])

    #Gets ranges
    start = temp.data$start
    end = temp.data$end

    #Obtains sequence
    new.seq = Biostrings::subseq(ref.genome, start = start, end = end)
    if (as.character(temp.data$strand) == "-"){ new.seq = Biostrings::reverseComplement(new.seq) }
    names(new.seq) = paste0(sprintf("%03d", i), "_", temp.data$type, "_", temp.data$name)

    if (nchar(names(new.seq)) >= 50){
      names(new.seq) = substr(names(new.seq), 1, 50)
    }

    save.seq = append(save.seq, new.seq)
  }#end i loop

  #Writes reference loci
  write.loci = as.list(as.character(save.seq))
  PhyloProcessR::writeFasta(sequences = write.loci, names = names(write.loci),
                       paste0(reference.name, "/refMarkers.fa"), nbchar = 1000000, as.string = TRUE)

  #WRites reference mito genome
  write.loci = as.list(as.character(ref.genome))
  PhyloProcessR::writeFasta(sequences = write.loci, names = names(write.loci),
                       paste0(reference.name, "/refGenome.fa"), nbchar = 1000000, as.string = TRUE)

  write.table(all.data, file = paste0(reference.name, "/referenceTable.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

  message("Mitochondrial reference for '", reference.name, "' created successfully.")

  return(invisible(NULL))
}#end function
