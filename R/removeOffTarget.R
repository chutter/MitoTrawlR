#' @title removeOffTarget
#'
#' @description Function for running the program spades to assemble short read sequencing data
#'
#' @param target path to a folder of sequence alignments in phylip format.
#'
#' @param contigs contigs are added into existing alignment if algorithm is "add"
#'
#' @param blast.path contigs are added into existing alignment if algorithm is "add"
#'
#' @param threads algorithm to use: "add" add sequences with "add.contigs"; "localpair" for local pair align. All others available
#'
#' @param quiet TRUE to supress screen output

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

removeOffTarget = function(target = NULL,
                           contigs = NULL,
                           blast.path = "blast",
                           threads = 1,
                           remove.bad = F,
                           quiet = TRUE) {

  # contigs = combined.contigs
  # target = reference
  # blast.path = "/Users/chutter/miniconda3/bin/"
  # threads = 6
  # quiet = F
  # remove.bad = T

  #Writes contigs for cap3
  write.loci = as.list(as.character(contigs))
  PhyloCap::writeFasta(sequences = write.loci, names = names(write.loci),
                       "blast_contigs.fa", nbchar = 1000000, as.string = T)

  #headers
  headers = c("qName", "tName", "pident", "matches", "misMatches", "gapopen",
              "qStart", "qEnd", "tStart", "tEnd", "evalue", "bitscore", "qLen", "tLen", "gaps")

  system(paste0(blast.path, "makeblastdb -in ", target, " -parse_seqids -dbtype nucl",
                " -out blast_db"), ignore.stdout = quiet, ignore.stderr = quiet)

  #Matches samples to loci
  system(paste0(blast.path, "blastn -task dc-megablast -db blast_db",
                " -query blast_contigs.fa -out blast_match.txt",
                " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen gaps\" ",
                " -num_threads ", threads))

  #Need to load in transcriptome for each species and take the matching transcripts to the database
  if (length(readLines("blast_match.txt")) == 0) { return(contigs) }
  match.data = read.table("blast_match.txt", sep = "\t", header = F, stringsAsFactors = FALSE)
  colnames(match.data) = headers
  system(paste0("rm blast_match.txt blast_db* blast_contigs.fa"))

  if (nrow(match.data) == 0){ return(contigs) }#end if

  if (remove.bad == TRUE){
    temp.data = data.frame(contig = match.data$qName, matches = match.data$matches)
    temp.agg = aggregate(temp.data$matches, by = list(temp.data$contig), FUN = "sum")
    colnames(temp.agg) = c("contig", "length")
    temp.agg = temp.agg[temp.agg$length %in% max(temp.agg$length),]
    save.contigs = contigs[names(contigs) %in% temp.agg$contig]
  } else { save.contigs = contigs[names(contigs) %in% match.data$qName] }

  #Grabs only the matches

  return(save.contigs)

}#end function
