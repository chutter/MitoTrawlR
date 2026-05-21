#' @title tRNAscanBatch
#'
#' @description Runs tRNAscan-SE on a batch of genome/contig FASTA files and
#'   saves the per-sample tRNAscan-SE output tables to \code{out.dir}. A
#'   summary CSV reporting the number of contigs and tRNAs found per sample is
#'   written to \code{logs/tRNAscan-sample-summary.csv}.
#'
#' @param genome.dir path to a folder of genome or contig FASTA files (one
#'   \code{.fa} per sample).
#'
#' @param genbank.file unused legacy parameter; retained for compatibility.
#'
#' @param out.dir path to the output directory where per-sample tRNAscan-SE
#'   result files will be written.
#'
#' @param tRNAscan.path path to the directory containing \code{tRNAscan-SE}.
#'   NULL uses the system PATH.
#'
#' @param organism.type tRNAscan-SE organism model: one of \code{"mammal"},
#'   \code{"vertebrate"}, or \code{"eukaryotic"}.
#'
#' @param overwrite logical; if TRUE, the output directory is deleted and
#'   recreated before processing.
#'
#' @param quiet logical; if TRUE, tRNAscan-SE screen output is suppressed.
#'
#' @return Invisibly returns NULL. Writes per-sample tRNAscan-SE text files to
#'   \code{out.dir} and a summary CSV to \code{logs/}.
#'
#' @examples
#' \dontrun{
#' tRNAscanBatch(
#'   genome.dir    = "draftContigs",
#'   out.dir       = "tRNAscan",
#'   tRNAscan.path = "/path/to/trnascan/bin",
#'   organism.type = "vertebrate",
#'   overwrite     = FALSE,
#'   quiet         = TRUE
#' )
#' }
#'
#' @export

tRNAscanBatch = function(genome.dir = NULL,
                         genbank.file = NULL,
                         out.dir = "tRNAscan",
                         tRNAscan.path = NULL,
                         organism.type = c("mammal", "vertebrate", "eukaryotic"),
                         overwrite = FALSE,
                         quiet = TRUE) {

  organism.type = match.arg(organism.type)

  if (!is.null(tRNAscan.path)){
    b.string = unlist(strsplit(tRNAscan.path, ""))
    if (b.string[length(b.string)] != "/") {
      tRNAscan.path = paste0(append(b.string, "/"), collapse = "")
    }
  } else { tRNAscan.path = "" }

  if (dir.exists(out.dir) == FALSE) {
    dir.create(out.dir)
  } else if (overwrite) {
    unlink(out.dir, recursive = TRUE)
    dir.create(out.dir)
  } else { message("Output directory already exists and overwrite = FALSE. Skipping."); return(invisible(NULL)) }

  if (organism.type == "eukaryotic"){ org.string = "-E" }
  if (organism.type == "vertebrate"){ org.string = "-M vert" }
  if (organism.type == "mammal"){ org.string = "-M mammal" }
  if (quiet){ org.string = paste0("-q ", org.string) }

  spp.samples = list.files(genome.dir)
  spp.samples = gsub(".fa$", "", spp.samples)

  save.data = data.frame(Sample = as.character(),
                         N_contigs = as.numeric(),
                         N_tRNA = as.numeric())

  for (i in seq_along(spp.samples)){

    input.file  = paste0(genome.dir, "/", spp.samples[i], ".fa")
    output.name = paste0(out.dir, "/", spp.samples[i], "_trnascan.txt")

    system(paste0(tRNAscan.path, "tRNAscan-SE ", org.string, " -o ", output.name, " ", input.file))

    if (!file.exists(output.name) || length(readLines(output.name, warn = FALSE)) <= 2) { next }

    scan.headers = c("contig", "trna_no", "start", "end", "type", "codon",
                     "codon_start", "codon_end", "score")
    scan.table = read.table(output.name, header = FALSE, skip = 3)
    colnames(scan.table) = scan.headers

    neg.strand = scan.table$start > scan.table$end
    scan.table$dir = ifelse(neg.strand, "-", "+")
    new.start = ifelse(neg.strand, scan.table$end,   scan.table$start)
    new.end   = ifelse(neg.strand, scan.table$start, scan.table$end)
    scan.table$start = new.start
    scan.table$end   = new.end

    temp.data = data.frame(Sample    = spp.samples[i],
                           N_contigs = length(unique(scan.table$contig)),
                           N_tRNA    = nrow(scan.table))
    save.data = rbind(save.data, temp.data)

  }

  dir.create("logs", recursive = TRUE, showWarnings = FALSE)
  write.csv(save.data, file = "logs/tRNAscan-sample-summary.csv", row.names = FALSE)

  return(invisible(NULL))

}#end function
