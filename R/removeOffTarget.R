#' @title removeOffTarget
#'
#' @description Filters assembled contigs by BLASTing them against a target
#'   reference FASTA and retaining only contigs that produce a match. When
#'   \code{remove.bad = TRUE}, only the contig with the highest total BLAST
#'   match length is kept; otherwise all contigs with any match are kept. If no
#'   BLAST hits are found the input contigs are returned unchanged.
#'
#' @param target path to a FASTA file used as the BLAST database (e.g., the
#'   reference mitochondrial genome).
#'
#' @param contigs a \code{DNAStringSet} of assembled contigs to filter.
#'
#' @param blast.path path to the directory containing the BLAST executables.
#'   NULL uses the system PATH.
#'
#' @param threads number of CPU threads to pass to \code{blastn}.
#'
#' @param remove.bad logical; if TRUE, keep only the single contig with the
#'   most total BLAST match bases; if FALSE, keep all contigs with any match.
#'
#' @param quiet logical; if TRUE, BLAST screen output is suppressed.
#'
#' @return A \code{DNAStringSet} containing the filtered contigs.
#'
#' @examples
#' \dontrun{
#' filtered <- removeOffTarget(
#'   target     = "reference/refGenome.fa",
#'   contigs    = assembled_contigs,
#'   blast.path = "/path/to/blast/bin",
#'   threads    = 4,
#'   remove.bad = FALSE,
#'   quiet      = TRUE
#' )
#' }
#'
#' @export

removeOffTarget = function(target = NULL,
                           contigs = NULL,
                           blast.path = NULL,
                           threads = 1,
                           remove.bad = FALSE,
                           quiet = TRUE) {

  if (!is.null(blast.path)){
    b.string = unlist(strsplit(blast.path, ""))
    if (b.string[length(b.string)] != "/") {
      blast.path = paste0(append(b.string, "/"), collapse = "")
    }
  } else { blast.path = "" }

  write.loci = as.list(as.character(contigs))
  PhyloProcessR::writeFasta(sequences = write.loci, names = names(write.loci),
                       "blast_contigs.fa", nbchar = 1000000, as.string = TRUE)

  headers = c("qName", "tName", "pident", "matches", "misMatches", "gapopen",
              "qStart", "qEnd", "tStart", "tEnd", "evalue", "bitscore", "qLen", "tLen", "gaps")

  system(paste0(blast.path, "makeblastdb -in ", target, " -parse_seqids -dbtype nucl",
                " -out blast_db"), ignore.stdout = quiet, ignore.stderr = quiet)

  system(paste0(blast.path, "blastn -task dc-megablast -db blast_db",
                " -query blast_contigs.fa -out blast_match.txt",
                " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen gaps\" ",
                " -num_threads ", threads),
         ignore.stdout = quiet, ignore.stderr = quiet)

  cleanup = function() {
    if (file.exists("blast_match.txt")){ file.remove("blast_match.txt") }
    if (file.exists("blast_contigs.fa")){ file.remove("blast_contigs.fa") }
    db.files = list.files(pattern = "^blast_db")
    if (length(db.files) > 0){ file.remove(db.files) }
  }

  if (!file.exists("blast_match.txt") || length(readLines("blast_match.txt", warn = FALSE)) == 0) {
    cleanup()
    return(contigs)
  }

  match.data = read.table("blast_match.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  colnames(match.data) = headers
  cleanup()

  if (nrow(match.data) == 0){ return(contigs) }

  if (remove.bad){
    temp.data = data.frame(contig = match.data$qName, matches = match.data$matches)
    temp.agg  = aggregate(temp.data$matches, by = list(temp.data$contig), FUN = "sum")
    colnames(temp.agg) = c("contig", "length")
    temp.agg  = temp.agg[temp.agg$length %in% max(temp.agg$length), ]
    save.contigs = contigs[names(contigs) %in% temp.agg$contig]
  } else {
    save.contigs = contigs[names(contigs) %in% match.data$qName]
  }

  return(save.contigs)

}#end function
