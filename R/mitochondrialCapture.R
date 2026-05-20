#' @title mitochondrialCapture
#'
#' @description Main entry point for iteratively assembling mitochondrial
#'   genomes from short reads for a batch of samples. For each sample, reads
#'   are matched to the reference genome by \code{iterativeAssemble} and the
#'   resulting draft contigs are saved as FASTA files in \code{output.dir}.
#'   Samples that have already been processed (i.e., whose output FASTA exists)
#'   are skipped when \code{overwrite = FALSE}.
#'
#' @param input.reads path to a folder of processed reads; each sample should
#'   reside in its own subdirectory containing fastq/fq files.
#'
#' @param reference.name name of the reference directory created by
#'   \code{buildReference}, which must contain \code{refGenome.fa}.
#'
#' @param output.dir path to the output directory where per-sample draft contig
#'   FASTA files (\code{<sample>.fa}) will be written.
#'
#' @param min.iterations minimum number of iterative assembly rounds per
#'   sample.
#'
#' @param max.iterations maximum number of iterative assembly rounds per
#'   sample.
#'
#' @param min.length minimum total assembled length (bp) below which off-target
#'   filtering is not applied.
#'
#' @param max.length maximum total assembled length (bp); triggers stricter
#'   read-identity filtering when exceeded.
#'
#' @param min.ref.id minimum fractional read-to-reference identity (0--1) for
#'   BBMap read recruitment.
#'
#' @param spades.path path to the directory containing \code{spades.py}. NULL
#'   uses the system PATH.
#'
#' @param bbmap.path path to the directory containing \code{bbmap.sh}. NULL
#'   uses the system PATH.
#'
#' @param cap3.path path to the directory containing \code{cap3}. NULL uses
#'   the system PATH.
#'
#' @param blast.path path to the directory containing the BLAST executables.
#'   NULL uses the system PATH.
#'
#' @param memory amount of RAM (GB) to allocate to BBMap and SPAdes.
#'
#' @param threads number of CPU threads to use.
#'
#' @param overwrite logical; if TRUE, already-processed samples are
#'   reprocessed and the output directory is recreated.
#'
#' @param quiet logical; if TRUE, external tool screen output is suppressed.
#'
#' @return Invisibly returns NULL. Writes one FASTA file per sample to
#'   \code{output.dir}.
#'
#' @examples
#' \dontrun{
#' mitochondrialCapture(
#'   input.reads    = "processed-reads/adaptor-removed-reads",
#'   reference.name = "reference",
#'   output.dir     = "draftContigs",
#'   min.iterations = 5,
#'   max.iterations = 20,
#'   min.length     = 15000,
#'   max.length     = 30000,
#'   min.ref.id     = 0.75,
#'   memory         = 8,
#'   threads        = 4,
#'   overwrite      = FALSE
#' )
#' }
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
#
#   input.reads = read.directory
#   reference.name = "reference"
#   output.dir = "draftContigs"
#   min.iterations = min.iterations
#   max.iterations = max.iterations
#   min.length = min.length
#   max.length = max.length
#   min.ref.id = min.read.match
#   memory = memory
#   threads = threads
#   spades.path = spades.path
#   bbmap.path = bbmap.path
#   cap3.path = cap3.path
#   blast.path = blast.path
#   overwrite = overwrite


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
  if (file.exists(input.reads) == FALSE){ stop("Input reads directory not found.") }

  if (dir.exists(reference.name) == FALSE){
    stop("Please provide a reference directory name made from the function buildReference.")
  }

  #Sets directory and reads
  if (dir.exists(output.dir) == FALSE) {
    dir.create(output.dir)
  } else {
    if (overwrite == TRUE){
      unlink(output.dir, recursive = TRUE)
      dir.create(output.dir)
    }
  }#end else

  #Creates output directory
  if (dir.exists("logs/sample_logs") == FALSE){ dir.create("logs/sample_logs", recursive = TRUE) }

  #Sets up the reads
  files = list.files(path = input.reads, full.names = TRUE, recursive = TRUE)
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
  if (overwrite == FALSE){
    done.names = list.files(output.dir)
    samples = samples[!samples %in% gsub(".fa$", "", done.names)]
  } else { samples = samples }

  if (length(samples) == 0){ stop("No samples to run or incorrect directory.") }

  #Header data for features and whatnot
  for (i in seq_along(samples)){

    sample.reads = reads[grep(paste0(samples[i], "_"), reads)]
    if (length(sample.reads) == 0){ sample.reads = reads[grep(samples[i], reads)] }

    if (length(sample.reads) == 0){
      warning(samples[i], " does not have any reads present. Skipping.")
      next
    }

    #Skip samples with empty or near-empty read files
    file.sizes = file.info(sample.reads)$size
    if (any(is.na(file.sizes)) || max(file.sizes, na.rm = TRUE) < 1000) {
      warning(samples[i], " read files are empty or near-empty (max file size: ",
              max(file.sizes, na.rm = TRUE), " bytes). Skipping.")
      next
    }

    #Concatenate together
    read1.reads = sample.reads[grep("_1.f.*|-1.f.*|_R1_.*|-R1_.*|_R1-.*|-R1-.*|READ1.*|_R1.fast.*|-R1.fast.*", sample.reads)]
    read2.reads = sample.reads[grep("_2.f.*|-2.f.*|_R2_.*|-R2_.*|_R2-.*|-R2-.*|READ2.*|_R2.fast.*|-R2.fast.*", sample.reads)]
    read3.reads = sample.reads[grep("_3.f.*|-3.f.*|_R3_.*|-R3_.*|_R3-.*|-R3-.*|READ3.*|_R3.fast.*|-R3.fast.*|_READ3.fast.*|-READ3.fast.*|_singleton.*|-singleton.*|READ-singleton.*|READ_singleton.*|_READ-singleton.*|-READ_singleton.*|-READ-singleton.*|_READ_singleton.*", sample.reads)]

    #Checks for reads
    if (length(read1.reads) == 0){
      warning(samples[i], ": read pairs could not be identified. Skipping.")
      next
    }
    if (length(read2.reads) == 0 && length(sample.reads) >= 2){
      warning(samples[i], ": read pairs could not be identified. Skipping.")
      next
    }
    if (length(read3.reads) == 0 && length(sample.reads) >= 3){
      warning(samples[i], ": merged/singleton reads could not be identified. Skipping.")
      next
    }

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
    } else if (length(read3.reads) == 1) {
      it.sample.reads[3] = read3.reads
    }


    #Runs iterative assembly function
    mito.contigs = iterativeAssemble(input.reads = it.sample.reads,
                                     reference = paste0(reference.name, "/refGenome.fa"),
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
                                     cap3.path = cap3.path,
                                     quiet = quiet,
                                     mapper = "bbmap")

    if (length(mito.contigs) == 0){
      message(samples[i], " failed: no reads matching to reference.")
      next
    }

    #Writes contigs
    names(mito.contigs) = paste0("seq", seq_along(mito.contigs))
    write.loci = as.list(as.character(mito.contigs))
    PhyloProcessR::writeFasta(sequences = write.loci, names = names(write.loci),
                         paste0(output.dir, "/", samples[i], ".fa"), nbchar = 1000000, as.string = TRUE)

    #Delete combined files
    if (length(sample.reads) > 3) {
      file.remove(Sys.glob(paste0(input.reads, "/", samples[i], "_ALL_READ*")))
    }#end delete

    message(samples[i], " completed!")

  }#end sample loop

  return(invisible(NULL))

}#end function


