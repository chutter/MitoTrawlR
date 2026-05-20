#' @title tRNAscan
#'
#' @description Runs tRNAscan-SE on a set of contigs to identify transfer RNA
#'   genes. Contigs are written to a temporary directory, tRNAscan-SE is
#'   called, the output table is parsed, and strand orientation is inferred
#'   from start/end coordinates. The temporary directory is deleted on
#'   completion. Returns NULL if no tRNAs are found.
#'
#' @param contigs a \code{DNAStringSet} of sequences to scan for tRNA genes.
#'
#' @param tRNAscan.path path to the directory containing \code{tRNAscan-SE}.
#'   NULL uses the system PATH.
#'
#' @param organism.type tRNAscan-SE organism model: one of \code{"mammal"},
#'   \code{"vertebrate"}, or \code{"eukaryotic"}.
#'
#' @param quiet logical; if TRUE, tRNAscan-SE output is suppressed via the
#'   \code{-q} flag.
#'
#' @return A data frame with columns \code{contig}, \code{trna_no},
#'   \code{start}, \code{end}, \code{type}, \code{codon},
#'   \code{codon_start}, \code{codon_end}, \code{score}, and \code{dir}, or
#'   NULL if no tRNAs were detected.
#'
#' @examples
#' \dontrun{
#' library(Biostrings)
#' contigs <- readDNAStringSet("draftContigs/sample1.fa")
#' rna.data <- tRNAscan(
#'   contigs       = contigs,
#'   tRNAscan.path = "/path/to/trnascan/bin",
#'   organism.type = "vertebrate",
#'   quiet         = TRUE
#' )
#' }
#'
#' @export

tRNAscan = function(contigs = NULL,
                    tRNAscan.path = NULL,
                    organism.type = c("mammal", "vertebrate", "eukaryotic"),
                    quiet = TRUE) {

  organism.type = match.arg(organism.type)

  if (!is.null(tRNAscan.path)){
    b.string = unlist(strsplit(tRNAscan.path, ""))
    if (b.string[length(b.string)] != "/") {
      tRNAscan.path = paste0(append(b.string, "/"), collapse = "")
    }
  } else { tRNAscan.path = "" }

  if (organism.type == "eukaryotic"){ org.string = "-E" }
  if (organism.type == "vertebrate"){ org.string = "-M vert" }
  if (organism.type == "mammal"){ org.string = "-M mammal" }
  if (quiet){ org.string = paste0("-q ", org.string) }

  if (dir.exists("temp-trna")){ unlink("temp-trna", recursive = TRUE) }
  dir.create("temp-trna")

  write.loci = as.list(as.character(contigs))
  PhyloProcessR::writeFasta(sequences = write.loci, names = names(write.loci),
             "temp-trna/input_contigs.fa", nbchar = 1000000, as.string = TRUE)

  system(paste0(tRNAscan.path, "tRNAscan-SE ", org.string, " -o temp-trna/output_trnascan.txt",
                " temp-trna/input_contigs.fa"))

  out.file = "temp-trna/output_trnascan.txt"
  if (!file.exists(out.file) || length(readLines(out.file, warn = FALSE)) <= 2) {
    unlink("temp-trna", recursive = TRUE)
    return(NULL)
  }

  scan.headers = c("contig", "trna_no", "start", "end", "type", "codon",
                   "codon_start", "codon_end", "score")
  scan.table = read.table(out.file, header = FALSE, skip = 3)
  colnames(scan.table) = scan.headers

  neg.strand = scan.table$start > scan.table$end
  scan.table$dir = ifelse(neg.strand, "-", "+")
  new.start = ifelse(neg.strand, scan.table$end,   scan.table$start)
  new.end   = ifelse(neg.strand, scan.table$start, scan.table$end)
  scan.table$start = new.start
  scan.table$end   = new.end

  unlink("temp-trna", recursive = TRUE)
  return(scan.table)

}#end function
