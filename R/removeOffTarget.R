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
#' @param min.coverage minimum fraction of a contig's length that must be
#'   covered by BLAST hits (summed) for the contig to be retained when
#'   \code{remove.bad = FALSE}. For example, \code{0.25} requires at least
#'   25\% of the contig to have a BLAST match. Default 0 keeps any contig
#'   with any hit (original behaviour).
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
                           min.coverage = 0,
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

  # Summarise total BLAST-matched bases per contig
  temp.agg = aggregate(match.data$matches,
                       by   = list(contig = match.data$qName),
                       FUN  = sum)
  colnames(temp.agg) = c("contig", "total.matches")

  # Attach query length (first occurrence per contig is sufficient)
  qlen.map = match.data[!duplicated(match.data$qName), c("qName", "qLen")]
  colnames(qlen.map) = c("contig", "qLen")
  temp.agg = merge(temp.agg, qlen.map, by = "contig")
  temp.agg$coverage = temp.agg$total.matches / temp.agg$qLen

  if (remove.bad){
    # Keep only the single contig with the most total BLAST match bases
    best = temp.agg[which.max(temp.agg$total.matches), "contig"]
    save.contigs = contigs[names(contigs) %in% best]
  } else {
    # Keep all contigs whose BLAST coverage meets the minimum threshold
    pass = temp.agg[temp.agg$coverage >= min.coverage, "contig"]
    save.contigs = contigs[names(contigs) %in% pass]
  }

  return(save.contigs)

}#end function
