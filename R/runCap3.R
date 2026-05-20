#' @title runCap3
#'
#' @description Runs the CAP3 sequence assembly program on a set of contigs to
#'   merge overlapping sequences into longer consensus contigs. Input can be
#'   supplied as a \code{DNAStringSet} (written to a temporary file) or as a
#'   path to an existing FASTA file. When \code{read.R = TRUE}, the assembled
#'   contigs and singlets are read back into R and returned as a
#'   \code{DNAStringSet}; otherwise results are optionally written to
#'   \code{output.name} and temporary files are deleted.
#'
#' @param contigs a \code{DNAStringSet} of sequences to assemble, or a
#'   character string giving the path to an existing FASTA file.
#'
#' @param output.name optional output file path. When provided, the CAP3
#'   contigs and singlets are concatenated and saved here.
#'
#' @param cap3.path path to the directory containing \code{cap3}. NULL uses
#'   the system PATH.
#'
#' @param read.R logical; if TRUE, read the assembled contigs and singlets back
#'   into R and return them as a \code{DNAStringSet}.
#'
#' @param a band expansion size (CAP3 \code{-a}; default 20).
#' @param b base quality cutoff for differences (CAP3 \code{-b}; default 20).
#' @param c base quality cutoff for clipping (CAP3 \code{-c}; default 12).
#' @param d max qscore sum at differences (CAP3 \code{-d}; default 200).
#' @param e clearance between number of differences (CAP3 \code{-e}; default 30).
#' @param f max gap length in any overlap (CAP3 \code{-f}; default 20).
#' @param g gap penalty factor (CAP3 \code{-g}; default 6).
#' @param h max overhang percent length (CAP3 \code{-h}; default 20).
#' @param i segment pair score cutoff (CAP3 \code{-i}; default 40).
#' @param j chain score cutoff (CAP3 \code{-j}; default 80).
#' @param k end clipping flag (CAP3 \code{-k}; default 1).
#' @param m match score factor (CAP3 \code{-m}; default 2).
#' @param n mismatch score factor (CAP3 \code{-n}; default -5).
#' @param o overlap length cutoff (CAP3 \code{-o}; default 40).
#' @param p overlap percent identity cutoff (CAP3 \code{-p}; default 90).
#' @param r reverse orientation value (CAP3 \code{-r}; default 1).
#' @param s overlap similarity score cutoff (CAP3 \code{-s}; default 900).
#' @param t max number of word matches (CAP3 \code{-t}; default 300).
#' @param u min number of constraints for correction (CAP3 \code{-u}; default 3).
#' @param v min number of constraints for linking (CAP3 \code{-v}; default 2).
#' @param y clipping range (CAP3 \code{-y}; default 100).
#' @param z min number of good reads at clip position (CAP3 \code{-z}; default 3).
#'
#' @return When \code{read.R = TRUE}, a \code{DNAStringSet} of assembled
#'   contigs and singlets. Otherwise invisibly returns NULL.
#'
#' @examples
#' \dontrun{
#' library(Biostrings)
#' seqs <- readDNAStringSet("fragments.fa")
#' merged <- runCap3(contigs = seqs, read.R = TRUE)
#' }
#'
#' @export

runCap3 = function(contigs = NULL,
                   output.name = NULL,
                   cap3.path = NULL,
                   read.R = FALSE,
                   a = 20,
                   b = 20,
                   c = 12,
                   d = 200,
                   e = 30,
                   f = 20,
                   g = 6,
                   h = 20,
                   i = 40,
                   j = 80,
                   k = 1,
                   m = 2,
                   n = -5,
                   o = 40,
                   p = 90,
                   r = 1,
                   s = 900,
                   t = 300,
                   u = 3,
                   v = 2,
                   y = 100,
                   z = 3) {

  if (!is.null(cap3.path)){
    b.string = unlist(strsplit(cap3.path, ""))
    if (b.string[length(b.string)] != "/") {
      cap3.path = paste0(append(b.string, "/"), collapse = "")
    }
  } else { cap3.path = "" }

  if (!is.character(contigs)) {
    write.loci = as.list(as.character(contigs))
    PhyloProcessR::writeFasta(sequences = write.loci, names = names(write.loci),
                         "input_contigs.fa", nbchar = 1000000, as.string = TRUE)
    contig.file = "input_contigs.fa"
  } else { contig.file = contigs }

  system(paste0(cap3.path, "cap3 ", contig.file, " -z ", z, " -o ", o, " -e ", e, " -s ", s, " > ",
                "input_contigs.fa.cap.txt"))

  if (!is.null(output.name)){
    system(paste0("cat ", contig.file, ".cap.contigs ",
                  contig.file, ".cap.singlets > ", output.name))
    file.remove(Sys.glob(paste0(contig.file, "*")))
    return(invisible(NULL))
  }

  if (read.R){
    temp.assembled = Biostrings::readDNAStringSet("input_contigs.fa.cap.contigs")
    temp.singlets  = Biostrings::readDNAStringSet("input_contigs.fa.cap.singlets")
    final.save = append(temp.assembled, temp.singlets)
    file.remove(Sys.glob("input_contigs.fa*"))
    return(final.save)
  }

  file.remove(Sys.glob(paste0(contig.file, "*")))
  return(invisible(NULL))

}#end function
