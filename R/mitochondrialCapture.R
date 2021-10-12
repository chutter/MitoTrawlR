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
                                reference.name = "reference",
                                output.dir = "draftContigs",
                                min.iterations = 5,
                                max.iterations = 20,
                                min.length = 17000,
                                max.length = 30000,
                                min.ref.id = 0.75,
                                spades.path = NULL,
                                bbmap.path = NULL,
                                cap3.path = NULL,
                                blast.path = NULL,
                                memory = 1,
                                threads = 1,
                                overwrite = FALSE,
                                quiet = TRUE) {

  # # #Debug
  # setwd("/Volumes/Rodents/Murinae/Mitochondrial_genomes")
  # input.reads = "/Volumes/Rodents/Murinae/processed-reads/adaptor-removed-reads"
  # reference.name = "reference"
  # output.dir = "draftContigs"
  # min.ref.id = 0.8
  # memory = 8
  # threads = 6
  # resume = TRUE
  # overwrite = FALSE
  # max.iterations = 30
  # min.iterations = 5
  # min.length = 15000
  # max.length = 40000
  # spades.path = "/usr/local/Spades/bin/spades.py"
  # bbmap.path = "/usr/local/bin/bbmap.sh"


  if (is.null(spades.path) == FALSE){
    b.string = unlist(strsplit(spades.path, ""))
    if (b.string[length(b.string)] != "/") {
      spades.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { spades.path = "" }

  #Same adds to bbmap path
  if (is.null(bbmap.path) == FALSE){
    b.string = unlist(strsplit(bbmap.path, ""))
    if (b.string[length(b.string)] != "/") {
      bbmap.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { bbmap.path = "" }

  #Same adds to bbmap path
  if (is.null(cap3.path) == FALSE){
    b.string = unlist(strsplit(cap3.path, ""))
    if (b.string[length(b.string)] != "/") {
      cap3.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { cap3.path = "" }

  #Same adds to bbmap path
  if (is.null(blast.path) == FALSE){
    b.string = unlist(strsplit(blast.path, ""))
    if (b.string[length(b.string)] != "/") {
      blast.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { blast.path = "" }

  #Quick checks
  if (is.null(input.reads) == TRUE){ stop("Please provide input reads.") }
  if (is.null(output.dir) == TRUE){ stop("Please provide an output directory.") }

  if (dir.exists(reference.name) == FALSE){
    stop("Please provide a reference directory name made from the function buildReference.")
  }

  #Sets directory and reads in  if (is.null(output.dir) == TRUE){ stop("Please provide an output directory.") }
  if (dir.exists(output.dir) == F){ dir.create(output.dir) } else {
    if (overwrite == TRUE){
      system(paste0("rm -r ", output.dir))
      dir.create(output.dir)
    }
  }#end else

  #Creates output directory
  if (dir.exists("logs") == F){ dir.create("logs") }

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

    #Concatenate together
    read1.reads = sample.reads[grep("_1.f.*|-1.f.*|_R1_.*|-R1_.*|_R1-.*|-R1-.*|READ1.*|_R1.fast.*|-R1.fast.*", sample.reads)]
    read2.reads = sample.reads[grep("_2.f.*|-2.f.*|_R2_.*|-R2_.*|_R2-.*|-R2-.*|READ2.*|_R2.fast.*|-R2.fast.*", sample.reads)]
    read3.reads = sample.reads[grep("_3.f.*|-3.f.*|_R3_.*|-R3_.*|_R3-.*|-R3-.*|READ3.*|_R3.fast.*|-R3.fast.*|_READ3.fast.*|-READ3.fast.*|_singleton.*|-singleton.*|READ-singleton.*|READ_singleton.*|_READ-singleton.*|-READ_singleton.*|-READ-singleton.*|_READ_singleton.*", sample.reads)]

    #Checks for reads
    if (length(read1.reads) == 0){ stop("error: read pairs could not be identified.")}
    if (length(read2.reads) == 0 && length(sample.reads) >= 2){ stop("error: read pairs could not be identified.")}
    if (length(read3.reads) == 0 && length(sample.reads) >= 3){ stop("error: merged set of reads could not be identified.")}

    #Combines duplicate read sets
    it.sample.reads = c()
    if (length(read1.reads) >= 2) {
      system(paste0("cat ", paste0(read1.reads, collapse = " "), " > ", input.reads,
                    "/", samples[i], "_ALL_READ1.fastq.gz"))
      it.sample.reads[1] = paste0(input.reads, "/", samples[i], "_ALL_READ1.fastq.gz")
    } else { it.sample.reads[1] = read1.reads }

    if (length(read2.reads) >= 2) {
      system(paste0("cat ", paste0(read2.reads, collapse = " "), " > ", input.reads,
                    "/", samples[i], "_ALL_READ2.fastq.gz"))
      it.sample.reads[2] = paste0(input.reads, "/", samples[i], "_ALL_READ2.fastq.gz")
    } else { it.sample.reads[2] = read2.reads }

    if (length(read3.reads) >= 2) {
      system(paste0("cat ", paste0(read3.reads, collapse = " "), " > ", input.reads,
                    "/", samples[i], "_ALL_READ3.fastq.gz"))
      it.sample.reads[3] = paste0(input.reads, "/", samples[i], "_ALL_READ3.fastq.gz")
    } else { it.sample.reads[3] = read3.reads }

    #Runs iterative assembly function
    mito.contigs = iterativeAssemble(input.reads = it.sample.reads,
                                     reference = paste0(reference.name, "/refGenome.fa"),
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
                                     blast.path = blast.path,
                                     mapper = "bbmap")

    if (length(mito.contigs) == 0){
      print(paste0(samples[i], " failed: no reads matching to reference."))
      next }

    #Writes contigs
    names(mito.contigs) = paste0("seq", rep(1:length(mito.contigs), by = 1))
    write.loci = as.list(as.character(mito.contigs))
    PhyloCap::writeFasta(sequences = write.loci, names = names(write.loci),
                         paste0(output.dir, "/", samples[i], ".fa"), nbchar = 1000000, as.string = T)

    #Delete combined files
    if (length(sample.reads) > 3) {
      system(paste0("rm ", input.reads, "/", samples[i], "_ALL_READ*"))
    }#end delete

    print(paste0(samples[i], " completed!"))

  }#end sample loop

}#end function


