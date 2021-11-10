#' @title isCircularGenome
#'
#' @description Function for removing adaptor sequences from raw Illumina sequence data using the program fastp
#'
#' @param contig path to a folder of raw reads in fastq format.
#'
#' @param genbank.reference a csv file with a "File" and "Sample" columns, where "File" is the file name and "Sample" is the desired renamed file
#'
#' @param blast.path the new directory to save the adaptor trimmed sequences
#'
#' @return a new directory of adaptor trimmed reads and a summary of the trimming in logs/
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

### DOESN"T WORK

#Iteratively assembles to reference
isCircularGenome = function(contig = NULL,
                            blast.path = "blastn",
                            cap3.path = "cap3") {

  #contig = contigs

  #Finds probes that match to two or more contigs
  options(stringsAsFactors = FALSE)

  if (is.null(contig) == TRUE){ stop("Please provide a genbank file.") }
  if (length(contig) != 1){ return(FALSE) }

  #Splits in two
  cut.contig1 = Biostrings::subseq(contig,
                                   start = 1,
                                   end = (Biostrings::width(contig)/2) )
  cut.contig2 = Biostrings::subseq(contig,
                                   start = (Biostrings::width(contig)/2) + 1,
                                   end = Biostrings::width(contig) )

  cut.contigs = append(cut.contig1, cut.contig2)
  names(cut.contigs) = c("seq1", "seq2")
  cap.contig = MitoCap::runCap3(contigs = cut.contigs,
                       a = 10,
                       b = 16,
                       c = 6,
                       d = 21,
                       e = 11,
                       o = 16,
                       h = 3,
                       i = 21,
                       j = 40,
                       s = 251,
                       u = 1,
                       v = 1,
                       y = 5,
                       z = 1)

  if (length(cap.contig) == 1){ return(TRUE) } else { return(FALSE) }

} #end function
