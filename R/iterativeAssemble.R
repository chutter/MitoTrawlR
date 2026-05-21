#' @title iterativeAssemble
#'
#' @description Iteratively assembles a target genomic region (e.g., the
#'   mitochondrial genome) from short reads by repeatedly mapping reads to the
#'   current assembly seed with BBMap and reassembling with SPAdes. Off-target
#'   contigs are removed with BLAST after each iteration. The loop terminates
#'   when the assembly length stops changing, the maximum number of iterations
#'   is reached, or the assembled length exceeds \code{max.length}.
#'
#' @param input.reads character vector of paths to read files. One element for
#'   single-end, two for paired-end, or three for paired-end plus merged/singleton
#'   reads.
#'
#' @param reference path to a FASTA file used as the initial mapping seed.
#'
#' @param mapper read mapper to use; currently only \code{"bbmap"} is
#'   supported.
#'
#' @param min.iterations minimum number of assembly iterations to complete
#'   before checking for convergence.
#'
#' @param max.iterations maximum number of assembly iterations allowed before
#'   the loop is forced to stop.
#'
#' @param min.length minimum total assembled length (bp) below which contigs
#'   are not filtered for off-target content.
#'
#' @param max.length maximum total assembled length (bp); when exceeded, the
#'   mapping identity threshold is raised and small contigs are discarded.
#'
#' @param target.length optional target length (bp); when the total assembled
#'   length reaches or exceeds this value the loop exits early. Useful for
#'   barcode scanning where you only need to cover a short reference region.
#'   \code{NULL} (default) disables early exit.
#'
#' @param min.ref.id minimum fractional read identity (0--1) required by BBMap
#'   for a read to be retained.
#'
#' @param memory amount of RAM (GB) to allocate to BBMap and SPAdes.
#'
#' @param threads number of CPU threads to use.
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
#' @param max.seed.contigs maximum number of contigs to carry in the seed at
#'   the end of each iteration. After BLAST filtering and deduplication the
#'   longest contigs are kept up to this limit, preventing unbounded growth.
#'   Default 30.
#'
#' @param min.blast.coverage minimum fraction (0--1) of a contig that must be
#'   covered by BLAST hits for it to be retained as on-target. Contigs with
#'   sparser matches (e.g. long NUMTs with a single short hit) are discarded.
#'   Default 0.25.
#'
#' @param quiet logical; if TRUE, external tool screen output is suppressed.
#'
#' @return A \code{DNAStringSet} of assembled contigs, or an empty
#'   \code{DNAStringSet} if assembly failed.
#'
#' @examples
#' \dontrun{
#' contigs <- iterativeAssemble(
#'   input.reads    = c("sample_R1.fastq.gz", "sample_R2.fastq.gz"),
#'   reference      = "reference/refGenome.fa",
#'   min.iterations = 5,
#'   max.iterations = 20,
#'   min.length     = 15000,
#'   max.length     = 30000,
#'   min.ref.id     = 0.75,
#'   memory         = 8,
#'   threads        = 4
#' )
#' }
#'
#' @export

#Iteratively assembles to reference
iterativeAssemble = function(input.reads = NULL,
                             reference = NULL,
                             mapper = "bbmap",
                             min.iterations = 5,
                             max.iterations = 20,
                             min.length = 16000,
                             max.length = 30000,
                             target.length = NULL,
                             min.ref.id = 0.75,
                             memory = 1,
                             threads = 1,
                             spades.path = NULL,
                             bbmap.path = NULL,
                             cap3.path = NULL,
                             blast.path = NULL,
                             max.seed.contigs = 30,
                             min.blast.coverage = 0.25,
                             quiet = TRUE) {

  # #Debug
  # input.reads = it.sample.reads
  # reference = paste0(output.directory, "/blast-reference/reference.fa")
  # min.ref.id = 0.7
  # max.iterations = 10
  # min.iterations = 3
  # min.length = Biostrings::width(ref.seq)
  # max.length = Biostrings::width(ref.seq)+(Biostrings::width(ref.seq) * 0.10)
  # spades.path = spades.path
  # bbmap.path = bbmap.path
  # #bbmap.path = "/Users/chutter/Bioinformatics/conda-envs/PhyloCap/bin/java /Users/chutter/Bioinformatics/conda-envs/PhyloCap/bin/"
  # blast.path = blast.path
  # cap3.path = "/Users/chutter/Bioinformatics/conda-envs/PhyloCap/bin"
  # mapper = "bbmap"
  # quiet = FALSE

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
  if (is.null(reference) == TRUE){ stop("Please provide a reference.") }

  #Writes reference to file if its not a file path
  if (file.exists("iterative_temp") == TRUE){ unlink("iterative_temp", recursive = TRUE) }
  if (file.exists("spades") == TRUE){ unlink("spades", recursive = TRUE) }
  dir.create("iterative_temp")
  system(paste0("cp ", reference, " iterative_temp/current_seed.fa"))

  #Sets up variables
  combined.contigs = Biostrings::DNAStringSet()
  new.len  = 0
  new.contigs = 0
  counter = 0
  repeat.counter = 0
  seeding = TRUE

  #############################
  ## While loop start
  #############################
  while (seeding == TRUE){

    #Copy new reference to do recursively
    counter = counter + 1
    prev.len = new.len
    prev.contigs = new.contigs
    #combined.contigs = DNAStringSet()

    message("------------   iteration ", counter, " begin ----------------------")

    #Subsets the reads to the lane
    set.reads = input.reads
    # set.read1 = set.reads[grep("READ1|R1", set.reads)]
    # if (length(input.reads) >= 2){ set.read2 = set.reads[grep("READ2|R2", set.reads)] }
    # if (length(input.reads) >= 3){ set.read3 = set.reads[grep("READ3|MERGE|singletons", set.reads)] }

    #Runs bbmap if selected
    if (mapper == "bbmap"){

      if (length(set.reads) == 1){
        #Pick out matching reads to mt Genomes
        system(paste0(bbmap.path, "bbmap.sh -Xmx", memory, "g ref=iterative_temp/current_seed.fa in1=", set.reads[1],
                      " vslow k=12 minid=", min.ref.id,
                      " outm1=iterative_temp/temp_read1.fq"),
               ignore.stderr = quiet, ignore.stdout = quiet)
        temp.read.path = "iterative_temp/temp_read1.fq"
      }#end 1 read

      if (length(set.reads) >= 2){
        #Pick out matching reads to mt Genomes
        system(paste0(bbmap.path, "bbmap.sh -Xmx", memory, "g ref=iterative_temp/current_seed.fa", " in1=", set.reads[1],
                      " in2=", set.reads[2], " vslow k=12 minid=", min.ref.id,
                      " outm1=iterative_temp/temp_read1.fq outm2=iterative_temp/temp_read2.fq"),
               ignore.stderr = quiet, ignore.stdout = quiet)
        temp.read.path = c("iterative_temp/temp_read1.fq",
                           "iterative_temp/temp_read2.fq")
      }#end 2 read

      if (length(set.reads) == 3){
        #Third set of reads
        system(paste0(bbmap.path, "bbmap.sh -Xmx", memory, "g ref=iterative_temp/current_seed.fa",
                      " in=", set.reads[3], " vslow k=12 minid=", min.ref.id,
                      " outm=iterative_temp/temp_read3.fq"),
               ignore.stderr = quiet, ignore.stdout = quiet)
        if (file.exists("iterative_temp/temp_read3.fq") &&
            file.info("iterative_temp/temp_read3.fq")$size > 0) {
          temp.read.path = c("iterative_temp/temp_read1.fq",
                             "iterative_temp/temp_read2.fq",
                             "iterative_temp/temp_read3.fq")
        }
      }#end set read 3

    }#end bbmap if

    unlink("ref", recursive = TRUE)

    #Runs spades
    spades.contigs = runSpades(read.paths = temp.read.path,
                               full.path.spades = spades.path,
                               mismatch.corrector = FALSE,
                               save.file = FALSE,
                               quiet = TRUE,
                               read.contigs = TRUE,
                               threads = threads,
                               memory = memory)

    spades.contigs = spades.contigs[Biostrings::width(spades.contigs) >= 200]

    #If SPAdes produced contigs, merge with accumulated; otherwise fall back to raw reads via CAP3
    if (length(spades.contigs) > 0) {
      combined.contigs = append(combined.contigs, spades.contigs)
    } else {
      save.seqs = Biostrings::DNAStringSet()
      for (x in seq_along(temp.read.path)){
        if (file.info(temp.read.path[x])$size == 0){ next }
        temp.fastq = Biostrings::readDNAStringSet(temp.read.path[x], format = "fastq")
        save.seqs = append(save.seqs, temp.fastq)
      }#end x
      if (length(save.seqs) >= 1000){ save.seqs = save.seqs[1:200] }
      if (length(save.seqs) != 0){
        save.seqs = append(save.seqs, combined.contigs)
        names(save.seqs) = paste0("seq", seq_along(save.seqs))
        combined.contigs = MitoTrawlR::runCap3(contigs = save.seqs,
                                               read.R = TRUE,
                                               cap3.path = cap3.path)
      } else { combined.contigs = Biostrings::DNAStringSet() }
    }

    #Checks for failure of everything
    if (length(combined.contigs) == 0) {
      seeding = FALSE
      message("mitogenome failed, reads could not be assembled.")
      next
    }

    #Cap3 to merge contigs from this iteration
    if (length(combined.contigs) > 1){
      combined.contigs = MitoTrawlR::runCap3(contigs = combined.contigs,
                                             read.R = TRUE,
                                             cap3.path = cap3.path)
    }

    #Remove contigs with insufficient BLAST coverage against the reference
    combined.contigs = removeOffTarget(target = reference,
                                       contigs = combined.contigs,
                                       blast.path = blast.path,
                                       threads = threads,
                                       quiet = TRUE,
                                       remove.bad = FALSE,
                                       min.coverage = min.blast.coverage)

    #Once we have enough length, keep only the single best-matching contig
    if (length(combined.contigs) > 0 &&
        max(Biostrings::width(combined.contigs)) >= min.length){
      combined.contigs = removeOffTarget(target = reference,
                                         contigs = combined.contigs,
                                         blast.path = blast.path,
                                         threads = threads,
                                         quiet = TRUE,
                                         remove.bad = TRUE)
    }#end if

    #Collapse near-identical / contained contigs to keep the seed manageable
    if (length(combined.contigs) > 1) {
      write.loci = as.list(as.character(combined.contigs))
      PhyloProcessR::writeFasta(sequences = write.loci, names = names(combined.contigs),
                           "iterative_temp/dedup_in.fa", nbchar = 1000000, as.string = TRUE)
      system(paste0(bbmap.path, "dedupe.sh -Xmx", memory, "g",
                    " in=iterative_temp/dedup_in.fa",
                    " out=iterative_temp/dedup_out.fa",
                    " minidentity=97 absorbcontainment=t"),
             ignore.stdout = quiet, ignore.stderr = quiet)
      if (file.exists("iterative_temp/dedup_out.fa") &&
          file.info("iterative_temp/dedup_out.fa")$size > 0) {
        combined.contigs = Biostrings::readDNAStringSet("iterative_temp/dedup_out.fa")
        names(combined.contigs) = paste0("seq", seq_along(combined.contigs))
      }
      file.remove(c("iterative_temp/dedup_in.fa", "iterative_temp/dedup_out.fa"))
    }

    #Hard cap: keep only the longest max.seed.contigs contigs for the next seed
    if (length(combined.contigs) > max.seed.contigs) {
      keep.idx = order(Biostrings::width(combined.contigs), decreasing = TRUE)[1:max.seed.contigs]
      combined.contigs = combined.contigs[keep.idx]
      names(combined.contigs) = paste0("seq", seq_along(combined.contigs))
      message("Seed capped at ", max.seed.contigs, " contigs.")
    }

    # Check for circularity: only when we have a single contig >= min.length
    # and have run at least min.iterations rounds. A circular genome cannot
    # grow further, so there is no reason to keep assembling.
    if (counter >= min.iterations && length(combined.contigs) == 1 &&
        Biostrings::width(combined.contigs) >= min.length) {
      if (MitoTrawlR::isCircularGenome(contig = combined.contigs,
                                       cap3.path = cap3.path)) {
        message("Circular genome detected! Stopping assembly.")
        seeding = FALSE
      }
    }

    #Check size
    new.len = sum(Biostrings::width(combined.contigs))
    new.contigs = length(combined.contigs)

    #If the file gets too large
    if (new.len >= max.length){
      #Merges first
      combined.contigs = combined.contigs[Biostrings::width(combined.contigs) >= 1000]

      new.len = sum(Biostrings::width(combined.contigs))
      new.contigs = length(combined.contigs)

      min.ref.id = 0.98
      #makes sure this doesn't go on forever and ever
      repeat.counter = repeat.counter + 1
      if (repeat.counter >= 5){
        message("repeat counter hit 5, stopping.")
        seeding = FALSE
      }#end if
    }

    #Prints completion info
    message("old length: ", prev.len, " bp; ", prev.contigs, " contigs.")
    message("new length: ", new.len, " bp; ", new.contigs, " contigs.")
    message("iteration ", counter, " complete!")
    message("-------------------------------------------------------")

    ###################################################
    #When target length is reached (e.g. barcode scanning)
    ###################################################
    if (!is.null(target.length) && new.len >= target.length) {
      seeding = FALSE
      message("Target length (", target.length, " bp) reached after ", counter, " iterations.")
    }

    ###################################################
    #When the counter gets too high
    ###################################################
    if (counter >= max.iterations){
      seeding = FALSE
      message("mitogenome complete after ", counter, " iterations!")
    } #end if

    ###################################################
    #When there is no change and min iterations is hit
    ###################################################
    if (new.len == prev.len & counter >= min.iterations){
      seeding = FALSE
      message("mitogenome complete after ", counter, " iterations!")
    } #end if

    #If there is nothing matching at this point
    if (length(combined.contigs) == 0){
      seeding = FALSE
      message("Could not find reads matching reference.")
    } else {
      #Writes contigs for next seed
      names(combined.contigs) = paste0("seq", seq_along(combined.contigs))
      write.loci = as.list(as.character(combined.contigs))
      PhyloProcessR::writeFasta(sequences = write.loci, names = names(write.loci),
                           "iterative_temp/current_seed.fa", nbchar = 1000000, as.string = TRUE)
    }#end else

  }#end while
  #############################
  ## While loop end
  #############################

  # #Saves the raw reads themselves
  # if (new.len <= 200){
  #   #Loops through each set of reads
  #   save.seqs = Biostrings::DNAStringSet()
  #   for (x in 1:length(temp.read.path)){
  #     #Check size
  #     temp.count = scan(file = temp.read.path[x], what = "character")
  #     if (length(temp.count) == 0){ next }
  #     #Sa es if there are reads
  #     temp.fastq = ShortRead::readFastq(temp.read.path[x])
  #     temp.fasta = temp.fastq@sread
  #     temp.seqs = Biostrings::DNAStringSet(unlist(lapply(temp.fasta, FUN = function (x) paste0(x, collapse = "") )))
  #     save.seqs = append(save.seqs, temp.seqs)
  #   }#end x
  #
  #   if (length(save.seqs) >= 1000){
  #     save.seqs = sample(save.seqs, 1000)
  #   }
  #
  #   #Saves them if there are any changes
  #   if (length(save.seqs) != 0){
  #     names(save.seqs) = paste0("seq", rep(1:length(save.seqs), by = 1))
  #     save.seqs = append(save.seqs, combined.contigs)
  #     combined.contigs = runCap3(contigs = save.seqs)
  #   }#end save.seqs if
  # }#end if

  unlink("iterative_temp", recursive = TRUE)

  if (length(combined.contigs) == 0){
    message("Could not find reads matching reference.")
  }

  return(combined.contigs)
  ##########################
}#end function
