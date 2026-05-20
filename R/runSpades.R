#' @title runSpades
#'
#' @description Runs SPAdes de-novo assembler on one, two, or three sets of
#'   reads. If assembly fails (no \code{contigs.fasta} output) the k-mer list
#'   is shortened by one and SPAdes is retried until either the assembly
#'   succeeds or all k-mer values are exhausted. When \code{read.contigs = TRUE}
#'   the assembled contigs are read back into R; when \code{save.file = TRUE}
#'   the contig FASTA is copied to \code{save.name}.
#'
#' @param read.paths character vector of paths to read files. One element for
#'   single-end reads, two for paired-end, or three for paired-end plus
#'   merged/singleton reads.
#'
#' @param full.path.spades path to the directory containing \code{spades.py}.
#'   NULL uses the system PATH.
#'
#' @param mismatch.corrector logical; if TRUE, passes \code{--careful} to
#'   SPAdes to enable mismatch correction (slower).
#'
#' @param kmer.values integer vector of k-mer sizes to attempt, tried in order
#'   from longest to shortest on failure.
#'
#' @param read.contigs logical; if TRUE, read assembled contigs into R and
#'   return them.
#'
#' @param save.file logical; if TRUE, copy the contig FASTA to \code{save.name}.
#'
#' @param save.name output file path (without extension) used when
#'   \code{save.file = TRUE}.
#'
#' @param threads number of CPU threads to pass to SPAdes.
#'
#' @param memory amount of RAM (GB) to allocate to SPAdes.
#'
#' @param overwrite logical; if TRUE, an existing \code{spades/} output
#'   directory is deleted before running.
#'
#' @param quiet logical; if TRUE, SPAdes stdout is suppressed.
#'
#' @return When \code{read.contigs = TRUE}, a \code{DNAStringSet} of assembled
#'   contigs (empty if assembly failed). When \code{save.file = TRUE} and
#'   \code{read.contigs = FALSE}, returns the string
#'   \code{"Contigs were saved to file."}. Otherwise returns
#'   \code{"Nothing was saved."}.
#'
#' @examples
#' \dontrun{
#' contigs <- runSpades(
#'   read.paths        = c("sample_R1.fastq.gz", "sample_R2.fastq.gz"),
#'   full.path.spades  = "/path/to/spades/bin",
#'   mismatch.corrector = FALSE,
#'   read.contigs      = TRUE,
#'   save.file         = FALSE,
#'   threads           = 4,
#'   memory            = 8
#' )
#' }
#'
#' @export

runSpades = function(read.paths = NULL,
                     full.path.spades = NULL,
                     mismatch.corrector = FALSE,
                     kmer.values = c(33,55,77,99,127),
                     read.contigs = TRUE,
                     save.file = TRUE,
                     save.name = NULL,
                     threads = 1,
                     memory = 4,
                     overwrite = FALSE,
                     quiet = TRUE) {

  if (!is.null(full.path.spades)){
    b.string = unlist(strsplit(full.path.spades, ""))
    if (b.string[length(b.string)] != "/") {
      full.path.spades = paste0(append(b.string, "/"), collapse = "")
    }
  } else { full.path.spades = "" }

  if (!file.exists(read.paths[1])){ stop("Read files not found.") }
  if (overwrite){
    if (dir.exists("spades")){ unlink("spades", recursive = TRUE) }
  }

  if (mismatch.corrector){ mismatch.string = "--careful " } else { mismatch.string = "" }

  k = kmer.values
  k.val = paste(k, collapse = ",")

  while (!file.exists("spades/contigs.fasta")){

    if (length(read.paths) == 1){
      system(paste0(full.path.spades, "spades.py --s1 ", read.paths[1],
                    " -o spades -k ", k.val, " ", mismatch.string, "-t ", threads, " -m ", memory),
             ignore.stdout = quiet)
    }

    if (length(read.paths) == 2){
      system(paste0(full.path.spades, "spades.py --pe1-1 ", read.paths[1], " --pe1-2 ", read.paths[2],
                    " -o spades -k ", k.val, " ", mismatch.string, "-t ", threads, " -m ", memory),
             ignore.stdout = quiet)
    }

    if (length(read.paths) == 3){
      system(paste0(full.path.spades, "spades.py --pe1-1 ", read.paths[1],
                    " --pe1-2 ", read.paths[2], " --merged ", read.paths[3],
                    " -o spades -k ", k.val, " ", mismatch.string, "-t ", threads, " -m ", memory),
             ignore.stdout = quiet)
    }

    k = k[-length(k)]
    if (length(k) == 0) { break }
    k.val = paste(k, collapse = ",")
  }

  if (length(k) == 0) {
    message("k-mer values all used up, cannot assemble!")
    unlink("spades", recursive = TRUE)
    return(Biostrings::DNAStringSet())
  }

  if (read.contigs){
    if (file.exists("spades/contigs.fasta")){
      contigs = Biostrings::readDNAStringSet("spades/contigs.fasta")
    } else {
      contigs = Biostrings::readDNAStringSet("spades/scaffolds.fasta")
    }
  }

  if (save.file){
    if (file.exists("spades/contigs.fasta")){
      system(paste0("cp spades/contigs.fasta ", save.name, ".fa"))
    } else {
      system(paste0("cp spades/scaffolds.fasta ", save.name, ".fa"))
    }
  }

  unlink("spades", recursive = TRUE)

  if (read.contigs){
    if (length(contigs) == 0){ message("No contigs were assembled.") }
    return(contigs)
  }

  if (save.file) { return("Contigs were saved to file.") }

  return("Nothing was saved.")

}# end spades function
