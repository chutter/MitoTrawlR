#' @title tRNAscan
#'
#' @description Function for running the program spades to assemble short read sequencing data
#'
#' @param genome.dir path to a folder of sequence alignments in phylip format.
#'
#' @param genbank.file contigs are added into existing alignment if algorithm is "add"
#'
#' @param out.dir contigs are added into existing alignment if algorithm is "add"
#'
#' @param tRNAscan.path algorithm to use: "add" add sequences with "add.contigs"; "localpair" for local pair align. All others available
#'
#' @param organism.type algorithm to use: "add" add sequences with "add.contigs"; "localpair" for local pair align. All others available
#'
#' @param overwrite TRUE to supress screen output
#'
#' @param quiet TRUE to supress screen output
#'
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

### Annotate tRNAs
tRNAscan = function(contigs = NULL,
                    tRNAscan.path = NULL,
                    organism.type = c("mammal", "vertebrate", "eukaryotic"),
                    quiet = TRUE) {

  #Debug
  # contigs = contigs
  # organism.type = "vertebrate"
  # quiet = TRUE
  # tRNAscan.path = trnascan.path

  if (is.null(tRNAscan.path) == FALSE){
    b.string = unlist(strsplit(tRNAscan.path, ""))
    if (b.string[length(b.string)] != "/") {
      tRNAscan.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { tRNAscan.path = "" }

  if (is.null(organism.type) == TRUE){ stop("Please select an organism type.") }

  if (organism.type == "eukaryotic"){ org.string = "-E" }
  if (organism.type == "vertebrate"){ org.string = "-M vert" }
  if (organism.type == "mammal"){ org.string = "-M mammal" }
  if (quiet == TRUE){ org.string = paste0("-q ", org.string) }

  #Set up input and output for tRNAscan
  if (dir.exists("temp-trna") == TRUE) { system("rm -r temp-trna") }
  dir.create("temp-trna")

  #Writes contigs for cap3
  write.loci = as.list(as.character(contigs))
  PhyloCap::writeFasta(sequences = write.loci, names = names(write.loci),
             "temp-trna/input_contigs.fa", nbchar = 1000000, as.string = T)

  #Runs tRNA scan
  system(paste0(tRNAscan.path, "tRNAscan-SE ", org.string, " -o temp-trna/output_trnascan.txt",
                " temp-trna/input_contigs.fa"))

  #Load in tRNAscan data
  if (length(readLines("temp-trna/output_trnascan.txt")) <= 2) { return(NULL) }
  scan.headers = c("contig", "trna_no", "start", "end", "type", "codon","codon_start", "codon_end", "score")
  scan.table = read.table("temp-trna/output_trnascan.txt", header = F, skip = 3)
  colnames(scan.table) = scan.headers

  #Fixes direction and adds into data
  scan.table$dir = as.character("0")
  #Finds out if they are overlapping
  for (k in 1:nrow(scan.table)){
    if (scan.table$start[k] > scan.table$end[k]){
      scan.table$dir[k] = "-"
      new.start = min(scan.table$start[k], scan.table$end[k])
      new.end = max(scan.table$start[k], scan.table$end[k])
      scan.table$start[k] = new.start
      scan.table$end[k] = new.end
    } else { scan.table$dir[k] = "+" }
  }#end k loop

  #Delete temp and return
  system("rm -r temp-trna")
  return(scan.table)

}#end function

