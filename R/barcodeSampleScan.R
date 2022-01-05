#' @title barcodeSampleScan
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
barcodeSampleScan = function(input.reads = NULL,
                             barcode.fasta = NULL,
                             output.directory = "barcodeScan",
                             barcode.database = c("GenBank", "File"),
                             database.file = NULL,
                             hits.per.sample = 5,
                             per.max.length = 0.25,
                             spades.path = NULL,
                             bbmap.path = NULL,
                             blast.path = NULL,
                             cap3.path = NULL,
                             memory = 1,
                             threads = 1,
                             overwrite = FALSE,
                             quiet = TRUE) {

  # #Debug
  # library(PhyloCap)
  # setwd("/Volumes/LaCie/Brygomantis/barcodeScan")
  # input.reads = "/Volumes/LaCie/Brygomantis/processed-reads/pe-merged-reads"
  # barcode.fasta = "/Volumes/LaCie/Brygomantis/barcodeScan/16S_barcode.fa"
  # output.directory = "barcodeScan"
  # memory = 36
  # threads = 8
  # overwrite = TRUE
  # spades.path = "/Users/chutter/Bioinformatics/conda-envs/PhyloCap/bin"
  # bbmap.path = "/Users/chutter/Bioinformatics/conda-envs/PhyloCap/bin/java /Users/chutter/Bioinformatics/conda-envs/PhyloCap/bin"
  # bbmap.path = "/Users/chutter/Bioinformatics/conda-envs/PhyloCap/bin"
  # blast.path = "/Users/chutter/Bioinformatics/conda-envs/PhyloCap/bin"
  # cap3.path = "/Users/chutter/Bioinformatics/conda-envs/PhyloCap/bin"
  # hits.per.sample = 5
  # per.max.length = 0.25
  # barcode.database = "GenBank"

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
  if (is.null(blast.path) == FALSE){
    b.string = unlist(strsplit(blast.path, ""))
    if (b.string[length(b.string)] != "/") {
      blast.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { blast.path = "" }

  #Same adds to bbmap path
  if (is.null(cap3.path) == FALSE){
    b.string = unlist(strsplit(cap3.path, ""))
    if (b.string[length(b.string)] != "/") {
      cap3.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { cap3.path = "" }

  #Quick checks
  if (is.null(input.reads) == TRUE){ stop("Please provide input reads.") }
  if (is.null(output.directory) == TRUE){ stop("Please provide an output directory.") }

  #Sets directory and reads in  if (is.null(output.dir) == TRUE){ stop("Please provide an output directory.") }
  if (dir.exists(output.directory) == F){ dir.create(output.directory) } else {
    if (overwrite == TRUE){
      system(paste0("rm -r ", output.directory))
      dir.create(output.directory)
    }
  }#end else

  #Creates output directory
  if (dir.exists("logs") == F){ dir.create("logs") }

  #Sets up the reads
  files = list.files(path = input.reads, full.names = T, recursive = T)
  reads = files[grep(pattern = "fastq|fq|clustS", x = files)]

  samples = gsub(paste0(input.reads, "/"), "", reads)
  samples = unique(gsub("/.*", "", samples))

  #Skips samples already finished
  if (overwrite == FALSE){
    done.names = list.files(output.directory)
    samples = samples[!samples %in% gsub(".fa$", "", done.names)]
  } else { samples = samples }

  if (length(samples) == 0){ stop("No samples to run or incorrect directory.") }

  #headers
  blast.headers = c("qName", "tName", "pident", "matches", "misMatches", "gapopen",
              "qStart", "qEnd", "tStart", "tEnd", "evalue", "bitscore", "qLen", "tLen", "gaps")

  dir.create(paste0(output.directory, "/sample-barcodes"))
  dir.create(paste0(output.directory, "/blast-reference"))
  system(paste0("cp ", barcode.fasta, " ", output.directory, "/blast-reference/reference.fa"))
  ref.seq = Biostrings::readDNAStringSet(paste0(output.directory, "/blast-reference/reference.fa"))
  system(paste0(blast.path, "makeblastdb -in ", barcode.fasta, " -parse_seqids -dbtype nucl ",
                " -out ", output.directory, "/blast-reference/barcode"))

  #Header data for features and whatnot
  all.barcodes = Biostrings::DNAStringSet()
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
    dir.create(paste0(output.directory, "/sample-barcodes/", samples[i]))
    mito.contigs = iterativeAssemble(input.reads = it.sample.reads,
                                     reference = paste0(output.directory, "/blast-reference/reference.fa"),
                                     min.ref.id = 0.7,
                                     memory = memory,
                                     threads = threads,
                                     max.iterations = 10,
                                     min.iterations = 5,
                                     min.length = Biostrings::width(ref.seq),
                                     max.length = Biostrings::width(ref.seq)+(Biostrings::width(ref.seq) * per.max.length),
                                     spades.path = spades.path,
                                     bbmap.path = bbmap.path,
                                     blast.path = blast.path,
                                     cap3.path = cap3.path,
                                     mapper = "bbmap")

    if (length(mito.contigs) == 0){
      print(paste0(samples[i], " failed: no reads matching to reference."))
      next }

    #Writes contigs
    names(mito.contigs) = paste0("seq", rep(1:length(mito.contigs), by = 1))
    write.loci = as.list(as.character(mito.contigs))
    PhyloCap::writeFasta(sequences = write.loci, names = names(write.loci),
                         paste0(output.directory, "/sample-barcodes/", samples[i], ".fa"), nbchar = 1000000, as.string = T)

    #Delete combined files
    if (length(sample.reads) > 3) {
      system(paste0("rm ", input.reads, "/", samples[i], "_ALL_READ*"))
    }#end delete

    #Matches samples to loci
    system(paste0(blast.path, "blastn -task dc-megablast -db ", output.directory, "/blast-reference/barcode",
                  " -query ", output.directory, "/sample-barcodes/", samples[i], ".fa",
                  " -out ", output.directory, "/", samples[i], "_match.txt",
                  " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen gaps\" ",
                  " -num_threads ", threads))

    #Need to load in transcriptome for each species and take the matching transcripts to the database
    match.data = fread(paste0(output.directory, "/", samples[i], "_match.txt"), sep = "\t", header = F, stringsAsFactors = FALSE)

    if (nrow(match.data) == 0){
      print("No matches found.")
      system(paste0("rm ", output.directory, "/", samples[i], "_match.txt"))
      next
    }#end if

    if (nrow(match.data) == 1){
      system(paste0("rm ", output.directory, "/", samples[i], "_match.txt"))
      #Writes the loci
      write.loci = as.list(as.character(mito.contigs))
      PhyloCap::writeFasta(sequences = write.loci, names = names(write.loci),
                  paste0(output.directory, "/sample-barcodes/", samples[i], ".fa"), nbchar = 1000000, as.string = T)
      #Saves into one file
      names(mito.contigs) = samples[i]
      all.barcodes = append(all.barcodes, mito.contigs)
      next
    }# end if

    #Filters to best
    setnames(match.data, blast.headers)

    #Deals with duplicate matches to same marker
    filt.data = match.data[match.data$bitscore == max(match.data$bitscore)][1]
    save.contigs = mito.contigs[names(mito.contigs) %in% filt.data$qName]

    #saves the file
    system(paste0("rm ", output.directory, "/", samples[i], "_match.txt"))
    #Writes the loci
    write.loci = as.list(as.character(save.contigs))
    PhyloCap::writeFasta(sequences = write.loci, names = names(write.loci),
                         paste0(output.directory, "/sample-barcodes/", samples[i], ".fa"), nbchar = 1000000, as.string = T)
    #Saves into one file
    names(save.contigs) = samples[i]
    all.barcodes = append(all.barcodes, save.contigs)

  }#end i loop

  #Writes the loci
  write.loci = as.list(as.character(all.barcodes))
  PhyloCap::writeFasta(sequences = write.loci, names = names(write.loci),
              paste0(output.directory, "/all-sample_barcodes.fa"), nbchar = 1000000, as.string = T)

  ########## USER PROVIDED DB
  if (barcode.database == "File"){
    dir.create(paste0(output.directory, "/barcode-database"))
    system(paste0("cp ", database.file, " ", output.directory, "/barcode-database/database.fa"))
    #Make blast db
    system(paste0(blast.path, "makeblastdb -in ", database.file, " -parse_seqids -dbtype nucl ",
                  " -out ", output.directory, "/barcode-database/database"))
    #blasting away
    system(paste0(blast.path, "blastn -task dc-megablast -db ", output.directory, "/barcode-database/database",
                  " -query ", output.directory, "/all-sample_barcodes.fa",
                  " -out ", output.directory, "/local-blast-hits.txt",
                  " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen gaps\" ",
                  " -num_threads ", threads))

    #Reads in blast results
    blast.headers = c("qName", "tName", "pident", "matches", "misMatches", "gapopen",
                      "qStart", "qEnd", "tStart", "tEnd", "evalue", "bitscore", "qLen", "tLen", "gaps")
    match.data = fread(paste0(output.directory, "/local-blast-hits.txt"), sep = "\t", header = F, stringsAsFactors = FALSE)
    setnames(match.data, blast.headers)

    ##### TO DO REST WITH EXAMPLE









    #############################

  }################## end Database if

  ########## BLAST to genbank
  if (barcode.database == "GenBank"){

    print("Now using remote blast. This may take some time....")
    #Blast assembled data using remote blast
    system(paste0(blast.path, "blastn -task blastn -db nt -query ",
                  output.directory, "/all-sample_barcodes.fa",
                  " -remote -out ", output.directory, "/remote-blast-hits.txt -max_target_seqs ", hits.per.sample, " -max_hsps 1",
                  " -outfmt \"6 qseqid sseqid sacc pident length mismatch evalue bitscore qlen slen stitle sscinames staxids qcovs\""))

    #Reads in blast results
    blast.headers = c("qseqid", "sseqid", "sacc", "pident", "length", "mismatch", "evalue",
                      "bitscore", "qlen", "slen","stitle","sscinames", "staxids", "qcovs")
    match.data = fread(paste0(output.directory, "/remote-blast-hits.txt"), sep = "\t", header = F, stringsAsFactors = FALSE)
    setnames(match.data, blast.headers)
    print("Remote blast finished!")

    ########## Create nice output of blast results
    header.data = c("Sample", "GenBank_ID", "GenBank_Order", "GenBank_Family", "GenBank_Species",
                    "pident", "length", "mismatch", "evalue", "bitscore", "coverage")

    #CReate empty dataset
    collect.data = data.table(matrix(as.numeric(0), nrow = length(samples)*hits.per.sample, ncol = length(header.data)))
    setnames(collect.data, header.data)
    collect.data[, Sample:=as.character(Sample)]
    collect.data[, GenBank_ID:=as.character(GenBank_ID)]
    collect.data[, GenBank_Species:=as.character(GenBank_Species)]
    collect.data[, GenBank_Family:=as.character(GenBank_Family)]
    collect.data[, GenBank_Order:=as.character(GenBank_Order)]

    index.val = as.integer(1)
    for (i in 1:length(samples)){

      sample.data = match.data[grep(gsub(".fa$", "", samples[i]), match.data$qseqid), ]

      #If no data
      if (nrow(sample.data) == 0){
        print("No GenBanks matches were found or not enough sequence to match. ")
        next
      }#end if

      #Goes through each sample and assigns genbank taxonomy
      for (j in 1:nrow(sample.data)){
        #Looks up genbank taxonomy somehow.
        system(paste0("curl https://taxonomy.jgi-psf.org/accession/", sample.data$sacc[j],
                      " > output.txt"))
        tax = readLines("output.txt", warn = FALSE)

        if (length(grep("Not found.", tax)) != 0) {
          set(collect.data, i = as.integer(index.val), j = match("Sample", header.data), value =  gsub(".fa$", "", samples[i]) )
          set(collect.data, i = as.integer(index.val), j = match("GenBank_ID", header.data), value = "Taxonomy-not-found" )
          index.val = index.val + as.integer(1)
          next
        }

        #Species data
        spp.line = tax[grep("\"species\":", tax)+1]
        spp.line = gsub("\"", "", spp.line)
        spp.line = gsub(",", "", spp.line)
        spp.line = gsub("\\.", "", spp.line)
        spp.line = gsub(".*name: ", "", spp.line)
        spp.line = gsub(" ", "_", spp.line)
        #Family data
        fam.line = tax[grep("\"family\":", tax)+1]
        fam.line = gsub("\"", "", fam.line)
        fam.line = gsub(",", "", fam.line)
        fam.line = gsub("\\.", "", fam.line)
        fam.line = gsub(".*name: ", "", fam.line)
        #Order Data
        ord.line = tax[grep("\"order\":", tax)+1]
        ord.line = gsub("\"", "", ord.line)
        ord.line = gsub(",", "", ord.line)
        ord.line = gsub("\\.", "", ord.line)
        ord.line = gsub(".*name: ", "", ord.line)

        if (length(ord.line) == 0){ord.line = "Not-Found" }

        unlink("output.txt")

        #Saves the data
        set(collect.data, i = as.integer(index.val), j = match("Sample", header.data), value =  gsub(".fa$", "", samples[i]) )
        set(collect.data, i = as.integer(index.val), j = match("GenBank_ID", header.data), value =  sample.data$sacc[j])
        set(collect.data, i = as.integer(index.val), j = match("GenBank_Order", header.data), value =  ord.line)
        set(collect.data, i = as.integer(index.val), j = match("GenBank_Species", header.data), value =  spp.line)
        set(collect.data, i = as.integer(index.val), j = match("GenBank_Family", header.data), value =  fam.line)
        set(collect.data, i = as.integer(index.val), j = match("pident", header.data), value =  sample.data$pident[j])
        set(collect.data, i = as.integer(index.val), j = match("length", header.data), value =  sample.data$length[j])
        set(collect.data, i = as.integer(index.val), j = match("mismatch", header.data), value =  sample.data$mismatch[j])
        set(collect.data, i = as.integer(index.val), j = match("evalue", header.data), value =  sample.data$evalue[j])
        set(collect.data, i = as.integer(index.val), j = match("bitscore", header.data), value =  sample.data$bitscore[j])
        set(collect.data, i = as.integer(index.val), j = match("coverage", header.data), value =  sample.data$qcovs[j])
        index.val = index.val + as.integer(1)
      }#end j loop
    }#end i loop

    #Finds the best match per each sample
    all.data = collect.data[collect.data$Sample != 0,]
    samples = unique(all.data$Sample)

    spp.data = c()
    for (k in 1:length(samples)){

      #Removes duplicate species from the same genbank code (spp data)
      tmp.data = all.data[all.data$Sample %in% samples[k],]
      #tmp.data = tmp.data[!duplicated(tmp.data$GenBank_Species),]
      #spp.data = rbind(spp.data, tmp.data)

      #Filters to the best matches
      best.data = tmp.data[tmp.data$bitscore == max(tmp.data$bitscore),]
      best.data = best.data[best.data$pident == max(best.data$pident),]
      best.data = best.data[best.data$evalue == min(best.data$evalue),]
      best.data = best.data[best.data$coverage == max(best.data$coverage),][1]
      spp.data = rbind(spp.data, best.data)

    }#end i loop

    #Saves the two datasets
    write.csv(all.data, paste0(output.directory, "/All-top-matches_Samples2Genbank.csv"), row.names = F)
    write.csv(spp.data, paste0(output.directory, "/Best-top-match_Samples2GenBank.csv"), row.names = F)

  }################## end GenBank if

}#end function

# END SCRIPT
