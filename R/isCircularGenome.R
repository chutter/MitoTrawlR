#' @title isCircularGenome
#'
#' @description Tests whether a single assembled contig represents a circular
#'   mitochondrial genome. Takes the last 25% (min 200 bp) of the contig as a
#'   tail and the first 25% as a head, then attempts to merge them with CAP3.
#'   If CAP3 produces a single contig the tail and head overlapped at the
#'   junction, indicating circularity.
#'
#' @param contig a \code{DNAStringSet} containing exactly one contig to test.
#'   Returns FALSE immediately if more than one contig is supplied or if the
#'   contig is shorter than 400 bp.
#'
#' @param cap3.path path to the CAP3 executable or its containing directory.
#'   Defaults to searching the system PATH.
#'
#' @return logical; TRUE if the contig appears to be circular, FALSE otherwise.
#'
#' @examples
#' \dontrun{
#' library(Biostrings)
#' contig <- readDNAStringSet("assembled_contig.fa")
#' isCircularGenome(contig = contig)
#' }
#'
#' @export

isCircularGenome = function(contig = NULL,
                            cap3.path = NULL) {

  if (is.null(contig)) { stop("Please provide a contig as a DNAStringSet.") }
  if (length(contig) != 1) { return(FALSE) }

  seq.len = Biostrings::width(contig)
  if (seq.len < 400) { return(FALSE) }

  overlap.size = max(200L, floor(seq.len * 0.25))

  tail.seq = Biostrings::subseq(contig, start = seq.len - overlap.size + 1L, end = seq.len)
  head.seq = Biostrings::subseq(contig, start = 1L, end = overlap.size)

  test.seqs = append(tail.seq, head.seq)
  names(test.seqs) = c("tail", "head")

  cap.result = MitoTrawlR::runCap3(contigs = test.seqs,
                                   read.R = TRUE,
                                   cap3.path = cap3.path)

  if (length(cap.result) == 1) { return(TRUE) } else { return(FALSE) }

} #end function
