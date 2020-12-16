#' @title mitochondrialCapture
#'
#' @description Function for removing adaptor sequences from raw Illumina sequence data using the program fastp
#'
#' @param input.reads path to a folder of raw reads in fastq format.
#'
#' @param genbank.file a csv file with a "File" and "Sample" columns, where "File" is the file name and "Sample" is the desired renamed file
#'
#' @param output.dir the new directory to save the adaptor trimmed sequences
#'
#' @param mapper "Sample" to run on a single sample or "Directory" to run on a directory of samples
#'
#' @param min.iterations system path to fastp in case it can't be found
#'
#' @param max.iterations system path to fastp in case it can't be found
#'
#' @param min.length system path to fastp in case it can't be found
#'
#' @param max.length system path to fastp in case it can't be found
#'
#' @param min.ref.id system path to fastp in case it can't be found
#'
#' @param spades.path system path to fastp in case it can't be found
#'
#' @param bbmap.path system path to fastp in case it can't be found
#'
#' @param cap3.path system path to fastp in case it can't be found
#'
#' @param threads number of computation processing threads
#'
#' @param memory amount of system memory to use
#'
#' @param resume TRUE to skip samples already completed
#'
#' @param overwrite TRUE to overwrite a folder of samples with output.dir
#'
#' @param quiet TRUE to supress screen output
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

#Iteratively assembles to reference
mitochondrialCapture = function(input.reads = NULL,
                                genbank.file = NULL,
                                output.dir = "draftContigs",
                                min.iterations = 5,
                                max.iterations = 20,
                                min.length = 17000,
                                max.length = 30000,
                                min.ref.id = 0.75,
                                spades.path = "spades.py",
                                bbmap.path = "bbmap.sh",
                                cap3.path = "cap3",
                                memory = 1,
                                threads = 1,
                                resume = TRUE,
                                overwrite = FALSE,
                                quiet = TRUE) {

  # # #Debug
  # input.reads = "pe-merged-reads"
  # genbank.file = "Crocidura.gb"
  # output.dir = "draftContigs"
  # min.ref.id = 0.80
  # memory = 8
  # threads = 6
  # resume = FALSE
  # overwrite = TRUE
  # max.iterations = 30
  # min.iterations = 5
  # min.length = 15000
  # max.length = 40000
  # spades.path = "/usr/local/Spades/bin/spades.py"
  # bbmap.path = "/usr/local/bin/bbmap.sh"

  #Quick checks
  options(stringsAsFactors = FALSE)
  if (is.null(input.reads) == TRUE){ stop("Please provide input reads.") }
  if (is.null(output.dir) == TRUE){ stop("Please provide an output directory.") }

  #Sets directory and reads in  if (is.null(output.dir) == TRUE){ stop("Please provide an output directory.") }
  if (dir.exists(output.dir) == F){ dir.create(output.dir) } else {
    if (overwrite == TRUE){
      system(paste0("rm -r ", output.dir))
      dir.create(output.dir)
    }
  }#end else

  #Creates output directory
  if (dir.exists("logs") == F){ dir.create("logs") }

  ### Makes the reference
  makeReference(genbank.file = genbank.file,
                overwrite = TRUE,
                rep.origin = FALSE)

  #Sets up the reads
  files = list.files(path = input.reads, full.names = T, recursive = T)
  reads = files[grep(pattern = "fastq|fq|clustS", x = files)]

  samples = gsub(paste0(input.reads, "/"), "", reads)
  samples = unique(gsub("/.*", "", samples))

  # if (sample.folder == T) { samples = list.dirs(read.dir, recursive = F, full.names = F) }
  # if (sample.folder == F) {
  #   samples = gsub("_R1.*", "", reads)
  #   samples = gsub("_R2.*", "", samples)
  #   samples = gsub("_READ1.*", "", samples)
  #   samples = gsub("_READ2.*", "", samples)
  #   samples = gsub("_singletons.*", "", samples)
  #   samples = gsub(".*\\/", "", samples)
  # }#end sample folder if

  #Skips samples already finished
  if (resume == TRUE){
    done.names = list.files(output.dir)
    samples = samples[!samples %in% gsub(".fa$", "", done.names)]
  } else { samples = samples }

  if (length(samples) == 0){ stop("No samples to run or incorrect directory.") }

  #Header data for features and whatnot
  for (i in 1:length(samples)){

    sample.reads = reads[grep(samples[i], reads)]

    mito.contigs = iterativeAssemble(input.reads = sample.reads,
                                     reference = "Mito-Reference/refGenome.fa",
                                     output.name = paste0(output.dir, "/", samples[i]),
                                     min.ref.id = min.ref.id,
                                     memory = memory,
                                     threads = threads,
                                     max.iterations = max.iterations,
                                     min.iterations = min.iterations,
                                     min.length = min.length,
                                     max.length = max.length,
                                     spades.path = spades.path,
                                     bbmap.path = bbmap.path,
                                     mapper = "bbmap")

    #Writes contigs
    names(mito.contigs) = paste0("seq", rep(1:length(mito.contigs), by = 1))
    write.loci = as.list(as.character(mito.contigs))
    writeFasta(sequences = write.loci, names = names(write.loci),
               paste0(output.dir, "/", samples[i], ".fa"), nbchar = 1000000, as.string = T)

    print(paste0(samples[i], " completed!"))

  }#end sample loop

}#end function
