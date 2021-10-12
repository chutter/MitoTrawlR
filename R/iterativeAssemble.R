#' @title iterativeAssemble
#'
#' @description Function for removing adaptor sequences from raw Illumina sequence data using the program fastp
#'
#' @param input.reads path to a folder of raw reads in fastq format.
#'
#' @param reference a csv file with a "File" and "Sample" columns, where "File" is the file name and "Sample" is the desired renamed file
#'
#' @param output.name the new directory to save the adaptor trimmed sequences
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
#' @param mem amount of system memory to use
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
iterativeAssemble = function(input.reads = NULL,
                             reference = NULL,
                             output.name = NULL,
                             mapper = "bbmap",
                             min.iterations = 5,
                             max.iterations = 20,
                             min.length = 16000,
                             max.length = 30000,
                             min.ref.id = 0.75,
                             memory = 1,
                             threads = 1,
                             spades.path = NULL,
                             bbmap.path = NULL,
                             cap3.path = NULL,
                             blast.path = NULL,
                             resume = TRUE,
                             overwrite = FALSE,
                             quiet = TRUE) {

  # #Debug
  # input.reads = sample.reads
  # reference = paste0(reference.name, "/refGenome.fa")
  # output.name = paste0(output.dir, "/", samples[i])
  # min.ref.id = 0.80
  # memory = 8
  # threads = 6
  # max.iterations = 20
  # min.iterations = 5
  # min.length = 16000
  # max.length = 40000
  # spades.path = "/usr/local/Spades/bin/spades.py"
  # bbmap.path = "/usr/local/bin/bbmap.sh"
  # mapper = "bbmap"
  #Same adds to bbmap path

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
  options(stringsAsFactors = FALSE)
  if (is.null(input.reads) == TRUE){ stop("Please provide input reads.") }
  if (is.null(reference) == TRUE){ stop("Please provide a reference.") }
  if (is.null(output.name) == TRUE){ stop("Please provide an output name.") }

  #Writes reference to file if its not a file path
  if(file.exists("iterative_temp") == TRUE){ system(paste0("rm -r iterative_temp")) }
  if(file.exists("spades") == TRUE){ system(paste0("rm -r spades")) }
  dir.create("iterative_temp")
  system(paste0("cp ", reference, " iterative_temp/current_seed.fa"))

  #Sets up variables
  combined.contigs = Biostrings::DNAStringSet()
  new.len  = 0
  new.contigs = 0
  counter = 0
  repeat.counter = 0
  seeding = T

  #Makes empty files
  system("touch iterative_temp/read1.fq")

  if (length(input.reads) >= 2){
    system("touch iterative_temp/read2.fq")
  }#end paired end

  if (length(input.reads) >= 3){
    system("touch iterative_temp/read3.fq")
  }#end paired end

  #############################
  ## While loop start
  #############################
  while (seeding == T){

    #Copy new reference to do recursively
    counter = counter + 1
    prev.len = new.len
    prev.contigs = new.contigs
    #combined.contigs = DNAStringSet()

    print(paste0("------------   iteration ", counter, " begin ----------------------"))

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
        #Saves and combines
        system(paste0("cat iterative_temp/read1.fq iterative_temp/temp_read1.fq >> iterative_temp/new_read1.fq"))
        system(paste0("rm iterative_temp/read1.fq"))
        system(paste0("mv iterative_temp/new_read1.fq iterative_temp/read1.fq"))
        system(paste0("rm iterative_temp/temp_read1.fq"))
        temp.read.path = paste0("iterative_temp/read1.fq")
      }#end 1 read

      if (length(set.reads) >= 2){
        #Pick out matching reads to mt Genomes
        system(paste0(bbmap.path, "bbmap.sh -Xmx", memory, "g ref=iterative_temp/current_seed.fa", " in1=", set.reads[1],
                      " in2=", set.reads[2], " vslow k=12 minid=", min.ref.id,
                      " outm1=iterative_temp/temp_read1.fq outm2=iterative_temp/temp_read2.fq"),
               ignore.stderr = quiet, ignore.stdout = quiet)
        #Saves and combines
        system(paste0("cat iterative_temp/read1.fq iterative_temp/temp_read1.fq >> iterative_temp/new_read1.fq"))
        system(paste0("cat iterative_temp/read2.fq iterative_temp/temp_read2.fq >> iterative_temp/new_read2.fq"))
        system(paste0("rm iterative_temp/read1.fq iterative_temp/read2.fq"))
        system(paste0("mv iterative_temp/new_read1.fq iterative_temp/read1.fq"))
        system(paste0("mv iterative_temp/new_read2.fq iterative_temp/read2.fq"))
        system(paste0("rm iterative_temp/temp_read1.fq iterative_temp/temp_read2.fq"))
        temp.read.path = c(paste0("iterative_temp/read1.fq"),
                           paste0("iterative_temp/read2.fq") )

      }#end 2 read

      if (length(set.reads) == 3){
        #Third set of reads
        system(paste0(bbmap.path, "bbmap.sh -Xmx", memory, "g ref=iterative_temp/current_seed.fa",
                      " in=", set.reads[3], " vslow k=12 minid=", min.ref.id,
                      " outm=iterative_temp/temp_read3.fq"),
               ignore.stderr = quiet, ignore.stdout = quiet)
        #Saves and combines
        system(paste0("cat iterative_temp/read3.fq iterative_temp/temp_read3.fq >> iterative_temp/new_read3.fq"))
        system(paste0("rm iterative_temp/read3.fq"))
        system(paste0("mv iterative_temp/new_read3.fq iterative_temp/read3.fq"))
        system(paste0("rm iterative_temp/temp_read3.fq"))
        temp.read.path = c(paste0("iterative_temp/read1.fq"),
                           paste0("iterative_temp/read2.fq"),
                           paste0("iterative_temp/read3.fq"))

        if (file.info("iterative_temp/read3.fq")$size == 0){
          temp.read.path = c(paste0("iterative_temp/read1.fq"),
                             paste0("iterative_temp/read2.fq"))
        }#end size if

      }#end set read 3

    }#end bbmap if

    system("rm -r ref")

    #Runs spades
    spades.contigs = runSpades(read.paths = temp.read.path,
                               full.path.spades = spades.path,
                               mismatch.corrector = FALSE,
                               save.file = F,
                               quiet =T,
                               read.contigs = T,
                               threads = threads,
                               memory = memory)

    spades.contigs = spades.contigs[Biostrings::width(spades.contigs) >= 100]

    #Saves the raw reads themselves
    if (length(spades.contigs) == 0){
      #Loops through each set of reads
      save.seqs = Biostrings::DNAStringSet()
      for (x in 1:length(temp.read.path)){
        #Check size
        temp.count = scan(file = temp.read.path[x], what = "character")
        if (file.info(temp.read.path[x])$size == 0){ next }
        #Sa es if there are reads
        temp.fastq = Biostrings::readDNAStringSet(temp.read.path[x], format = "fastq")
        save.seqs = append(save.seqs, temp.fastq)
      }#end x
      #Removes if too many
      if (length(save.seqs) >= 1000){ save.seqs = save.seqs[1:200] }

      #Saves them if there are any changes
      if (length(save.seqs) != 0){
        save.seqs = append(save.seqs, combined.contigs)
        names(save.seqs) = paste0("seq", rep(1:length(save.seqs), by = 1))
        combined.contigs = MitoCap::runCap3(contigs = save.seqs,
                                            read.R = TRUE)
      } else { combined.contigs = Biostrings::DNAStringSet() }
    } else { combined.contigs = append(combined.contigs, spades.contigs)}

    #Checks for failure of everything
    if (length(combined.contigs) == 0) {
      seeding = F
      print(paste0("mitogenome failed, reads could not be assembled."))
      next
    }

    #Cap3 to combine old and new
    combined.contigs = append(combined.contigs, spades.contigs)
    if (length(combined.contigs) > 1){
      combined.contigs = MitoCap::runCap3(contigs = combined.contigs,
                                          read.R = TRUE)
    }

    #Blasts and remvoes stuff that doesn't match to reference
    if (length(combined.contigs) >= 2){
      combined.contigs = removeOffTarget(target = reference,
                                         contigs = combined.contigs,
                                         blast.path = blast.path,
                                         threads = threads,
                                         quiet = T,
                                         remove.bad = F)
    }#end if

    #Removes even more off target bad matches
    if (max(Biostrings::width(combined.contigs)) >= min.length){
      combined.contigs = removeOffTarget(target = reference,
                                         contigs = combined.contigs,
                                         blast.path = blast.path,
                                         threads = threads,
                                         quiet = T,
                                         remove.bad = T)
    }#end if

    if (max(Biostrings::width(combined.contigs)) >= min.length){
      circ.contigs = combined.contigs[Biostrings::width(combined.contigs) >= min.length]
      circ.genome = isCircularGenome(circ.contigs)

      if (circ.genome == TRUE) {
        print(paste("Circular genome assembled!"))
        seeding = F
        combined.contigs = circ.contigs
      }
    }#end if

    #Check size
    new.len = sum(Biostrings::width(combined.contigs))
    new.contigs = length(combined.contigs)

    #If the file gets too large
    if (new.len >= max.length){
      #Merges first
      combined.contigs = combined.contigs[Biostrings::width(combined.contigs) >= 1000]

      new.len = sum(Biostrings::width(combined.contigs))
      new.contigs = length(combined.contigs)

      min.ref.id = "0.98"
      #makes sure this doesn't go on forever and ever
      repeat.counter = repeat.counter+1
      if (repeat.counter >= 5){
        print(paste("repeat counter hit 5"))
        seeding = F
      }#end if
    } else { min.ref.id = "0.95" }

    #Prints completion info
    print(paste0("old length: ", prev.len, " bp; ", prev.contigs, " contigs."))
    print(paste0("new length: ", new.len, " bp; ", new.contigs, " contigs."))
    print(paste0("iteration ", counter, " complete!"))
    print(paste0("-------------------------------------------------------"))

    ###################################################
    #When the counter gets too high
    ###################################################
    if (counter >= max.iterations){
      seeding = F
      print(paste0("mitogenome complete after ", counter, " iterations!"))
    } #end if

    ###################################################
    #When there is no change and min iterations is hit
    ###################################################
    if (new.len == prev.len & counter >= min.iterations){
      seeding = F
      print(paste0("mitogenome complete after ", counter, " iterations!"))
    } #end if

    #If there is nothing matching at this point
    if (length(combined.contigs) == 0){
      seeding = F
      print(paste0("Could not find reads matching reference."))
    } else {
      #Writes contigs for next seed
      names(combined.contigs) = paste0("seq", rep(1:length(combined.contigs), by = 1))
      write.loci = as.list(as.character(combined.contigs))
      PhyloCap::writeFasta(sequences = write.loci, names = names(write.loci),
                           "iterative_temp/current_seed.fa", nbchar = 1000000, as.string = T)
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

  #If there is nothing matching at this point
  if (length(combined.contigs) == 0){
    system("rm -r iterative_temp")
    print(paste0("Could not find reads matching reference."))
    return(combined.contigs)
  } else {

    system("rm -r iterative_temp")
    return(combined.contigs)
  }#end else
  ##########################
}#end function
