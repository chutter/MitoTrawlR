#' @title tRNAscanBatch
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
tRNAscanBatch = function(genome.dir = NULL,
                         genbank.file = NULL,
                         out.dir = "tRNAscan",
                         tRNAscan.path = "tRNAscan-SE",
                         organism.type = c("mammal", "vertebrate", "eukaryotic"),
                         overwrite = FALSE,
                         quiet = TRUE) {

  #Debug
  # genome.dir = "draftContigs"
  # genbank.file = gb.file
  # out.dir = "tRNAscan"
  # organism.type = "vertebrate"
  # overwrite = TRUE
  # tRNAscan.path = trnascan.path

  #Checks for output directory and overwriting
  if (dir.exists(out.dir) == FALSE) { dir.create(out.dir) }
  if (dir.exists(out.dir) == TRUE) {
    if (overwrite == TRUE){
      system(paste0("rm -r ", out.dir))
      dir.create(out.dir)
    } else { stop("overwrite is false and directory exists. Exiting.") }
  }#end dir exists

  if (is.null(organism.type) == TRUE){ stop("Please select an organism type.") }

  if (organism.type == "eukaryotic"){ org.string = "-E" }
  if (organism.type == "vertebrate"){ org.string = "-M vert" }
  if (organism.type == "mammal"){ org.string = "-M mammal" }
  if (quiet == TRUE){ org.string = paste0("-q ", org.string) }

  #Obtains samples
  spp.samples = list.files(genome.dir)
  spp.samples = gsub(".fa$", "", spp.samples)

  #Empty data.frame for save data
  save.data = data.frame(Sample = as.character(),
                         N_contigs = as.numeric(),
                         N_tRNA = as.numeric())
  for (i in 1:length(spp.samples)){

    #Set up input and output for tRNAscan
    input.file = paste0(genome.dir, "/", spp.samples[i], ".fa")
    output.name = paste0(out.dir, "/", spp.samples[i], "_trnascan.txt")

    #Runs tRNA scan
    system(paste0(tRNAscan.path, " ", org.string, " -o ", output.name, " ",
                  input.file ))

    #Load in tRNAscan data
    scan.headers = c("contig", "trna_no", "start", "end", "type", "codon","codon_start", "codon_end", "score")
    if (length(readLines(output.name)) == 0) { next }
    scan.table = read.table(output.name, header = F, skip = 3)
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

    #Saves summary data
    temp.data = data.frame(Sample = spp.samples[i],
                           N_contigs = length(unique(scan.table$contig)),
                           N_tRNA = nrow(scan.table))

    save.data = rbind(save.data, temp.data)

  }#end i loop

  #Writes the summary data
  write.csv(save.data, file = "logs/tRNAscan-sample-summary.csv")

}#end function

