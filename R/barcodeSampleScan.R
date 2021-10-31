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
                             output.directory = "draftContigs",
                             min.iterations = 5,
                             max.iterations = 20,
                             min.length = 17000,
                             max.length = 30000,
                             min.ref.id = 0.75,
                             spades.path = NULL,
                             bbmap.path = NULL,
                             blast.path = NULL,
                             memory = 1,
                             threads = 1,
                             overwrite = FALSE,
                             quiet = TRUE) {

  # #Debug
  setwd("/Volumes/Rodents/Mitogenomes")
  input.reads = "/Volumes/Rodents/Mitogenomes/Crocidura/processed-reads/pe-merged-reads"
  barcode.fasta = "/Volumes/Rodents/Mitogenomes/Crocidura/barcode_cox1.fa"
  output.directory = "draftContigs"
  min.ref.id = 0.8
  memory = 8
  threads = 6
  overwrite = FALSE
  max.iterations = 30
  min.iterations = 5
  min.length = 15000
  max.length = 40000
  spades.path = "/usr/local/Spades/bin/spades.py"
  bbmap.path = "/usr/local/bin/bbmap.sh"


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


  #Here

  ##Concateantes all the samples together for an easier blast run
  dir.create("sampleMarkers")

  #Gets samples that succeeded
  samples = list.files(".")
  samples = samples[grep(".fa$", samples)]
  samples = samples[samples != "reference.fa"]
  samples = samples[samples != "sample_blast_pieces.fa"]

  #Make blast database for the probe loci
  #headers
  headers = c("qName", "tName", "pident", "matches", "misMatches", "gapopen",
              "qStart", "qEnd", "tStart", "tEnd", "evalue", "bitscore", "qLen", "tLen", "gaps")

  system(paste0("makeblastdb -in reference.fa -parse_seqids -dbtype nucl ",
                " -out blast_db"))


  #Gets the new tree fiels it made
  header.data = c("Sample", "readPairs", "matchingReads", "numberTargets", "targetPer")
  collect.data = data.table::data.table(matrix(as.numeric(0), nrow = length(samples), ncol = length(header.data)))
  data.table::setnames(collect.data, header.data)
  collect.data[, Sample:=as.character(Sample)]

  #Header data for features and whatnot
  for (i in 1:length(samples)){
    #Gets reads together for this sample
    sample.reads = reads[grep(samples[i], reads)]
    #Gets sets of reads
    sample.sets = gsub(paste0(".*/", samples[i], "/"), "", sample.reads)


    contigs = scanFa(paste0(samples[i]))

    #Matches samples to loci
    system(paste0("blastn -task dc-megablast -db blast_db",
                  " -query ", samples[i], " -out ", samples[i], "_match.txt",
                  " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen gaps\" ",
                  " -num_threads ", threads))








    system(paste0("gzcat ", sample.sets[1], " > ",  gsub(".gz$", "",sample.sets[1]) ))
    system(paste0("gzcat ", sample.sets[2], " > ",  gsub(".gz$", "",sample.sets[2]) ))

    sample.reads = gsub(".gz$", "",sample.sets[1])
    sample.reads[2] = gsub(".gz$", "",sample.sets[2])

    #Magic blast stuff here
    #Make blast database for the probe loci\
    ref.contigs = scanFa(paste0(work.dir, "/", reference))

    system(paste0("makeblastdb -in ", work.dir, "/", reference,
                  " -parse_seqids -dbtype nucl -out transcript_nucl-blast_db"))

    #Matches samples to proteins
    system(paste0("magicblast -db transcript_nucl-blast_db -infmt fastq -reftype transcriptome",
                  " -query ", sample.reads[1], " -query_mate ", sample.reads[2],
                  " -out ", samples[i], "_prot-match.txt -perc_identity 0.65",
                  " -outfmt tabular -num_threads ", threads))

    # system(paste0("magicblast -db ", sample, "_nucl-blast_db -reftype transcriptome",
    #               " -query ", sample, "_dedupe.fa -out ", sample, "_prot-match.txt",
    #               " -outfmt tabular -num_threads ", threads))

    #headers for the blast db
    headers = c("qName", "tName", "pident", "qStart", "qEnd", "tStart", "tEnd", "evalue",
                "qDir", "tDir", "qLen", "BTOP", "numPl", "compart", "lOver", "rOver")

    match.data = fread(paste0(samples[i], "_prot-match.txt"), sep = "\t", header = T, stringsAsFactors = FALSE)
    match.data[,grep("not used", colnames(match.data)) :=  NULL]
    match.data[, c((ncol(match.data)-2):ncol(match.data)):=NULL]
    setnames(match.data, headers)

    #Count stuff
    good.data = match.data[match.data$tName != "-",]
    read.pairs = nrow(match.data)
    no.targets = length(unique(good.data$tName))
    match.reads = nrow(good.data)
    per.targets = no.targets/length(ref.contigs)

    #Saves data
    set(collect.data, i = as.integer(i), j = match("Sample", header.data), value = samples[i] )
    set(collect.data, i = as.integer(i), j = match("readPairs", header.data), value = read.pairs )
    set(collect.data, i = as.integer(i), j = match("matchingReads", header.data), value = match.reads )
    set(collect.data, i = as.integer(i), j = match("numberTargets", header.data), value = no.targets )
    set(collect.data, i = as.integer(i), j = match("targetPer", header.data), value = per.targets )

    #Delete
    system(paste0("rm ", sample.reads[1]))
    system(paste0("rm ", sample.reads[2]))





#### ODL

#Sets up the reads
files = list.files(path = read.dir, full.names = T, recursive = T)
reads = files[grep(pattern = "fastq|fq", x = files)]
if (is.na(sub.folder) != T){ reads = reads[grep(pattern = sub.folder, x = reads)] }

#Creates directories and copies files
system(paste0("cp ", work.dir, "/", reference, " ", work.dir, "/", out.dir, "/reference.fa"))
setwd(paste0(work.dir, "/", out.dir))

if (sample.folder == T) { samples = list.dirs(read.dir, recursive = F, full.names = F) }
if (sample.folder == F) {
  samples = gsub("_R1.*", "", reads)
  samples = gsub("_R2.*", "", samples)
  samples = gsub(".*\\/", "", samples)
}#end sample folder if

samples = unique(gsub("_L7_.*", "", samples))
if (length(samples) == 0){ stop("No samples to run or incorrect directory.") }



#Header data for features and whatnot
for (i in 1:length(samples)){

  #Gets reads together for this sample
  sample.reads = reads[grep(samples[i], reads)]
  #Gets sets of reads
  sample.sets = gsub(paste0(".*/", samples[i], "/"), "", sample.reads)

  system(paste0("gzcat ", sample.sets[1], " > ",  gsub(".gz$", "",sample.sets[1]) ))
  system(paste0("gzcat ", sample.sets[2], " > ",  gsub(".gz$", "",sample.sets[2]) ))

  sample.reads = gsub(".gz$", "",sample.sets[1])
  sample.reads[2] = gsub(".gz$", "",sample.sets[2])

  #Magic blast stuff here
  #Make blast database for the probe loci\
  ref.contigs = scanFa(paste0(work.dir, "/", reference))

  system(paste0("makeblastdb -in ", work.dir, "/", reference,
                " -parse_seqids -dbtype nucl -out transcript_nucl-blast_db"))

  #Matches samples to proteins
  system(paste0("magicblast -db transcript_nucl-blast_db -infmt fastq -reftype transcriptome",
                " -query ", sample.reads[1], " -query_mate ", sample.reads[2],
                " -out ", samples[i], "_prot-match.txt -perc_identity 0.65",
                " -outfmt tabular -num_threads ", threads))

  # system(paste0("magicblast -db ", sample, "_nucl-blast_db -reftype transcriptome",
  #               " -query ", sample, "_dedupe.fa -out ", sample, "_prot-match.txt",
  #               " -outfmt tabular -num_threads ", threads))

  #headers for the blast db
  headers = c("qName", "tName", "pident", "qStart", "qEnd", "tStart", "tEnd", "evalue",
              "qDir", "tDir", "qLen", "BTOP", "numPl", "compart", "lOver", "rOver")

  match.data = fread(paste0(samples[i], "_prot-match.txt"), sep = "\t", header = T, stringsAsFactors = FALSE)
  match.data[,grep("not used", colnames(match.data)) :=  NULL]
  match.data[, c((ncol(match.data)-2):ncol(match.data)):=NULL]
  setnames(match.data, headers)

  #Count stuff
  good.data = match.data[match.data$tName != "-",]
  read.pairs = nrow(match.data)
  no.targets = length(unique(good.data$tName))
  match.reads = nrow(good.data)
  per.targets = no.targets/length(ref.contigs)

  #Saves data
  set(collect.data, i = as.integer(i), j = match("Sample", header.data), value = samples[i] )
  set(collect.data, i = as.integer(i), j = match("readPairs", header.data), value = read.pairs )
  set(collect.data, i = as.integer(i), j = match("matchingReads", header.data), value = match.reads )
  set(collect.data, i = as.integer(i), j = match("numberTargets", header.data), value = no.targets )
  set(collect.data, i = as.integer(i), j = match("targetPer", header.data), value = per.targets )

  #Delete
  system(paste0("rm ", sample.reads[1]))
  system(paste0("rm ", sample.reads[2]))

}#end i loop

write.csv(collect.data, "read-map_Limno.csv")

  #Iteratively assembles mitochondrial genome
  assembled.mt = iterative.assemble(read.names = sample.sets, reference.map = reference,
                                    output.name = samples[i], paired.end = pe.reads,
                                    min.iterations = 5, max.iterations = 20, max.length = 20000,
                                    min.ref.id = min.id, memory = mem, cores = threads)

}#end i

###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################

##Concateantes all the samples together for an easier blast run
setwd(paste0(work.dir, "/", out.dir))
dir.create("sampleMarkers")

#Gets samples that succeeded
samples = list.files(".")
samples = samples[grep(".fa$", samples)]
samples = samples[samples != "reference.fa"]
samples = samples[samples != "sample_blast_pieces.fa"]

#Make blast database for the probe loci
#headers
headers = c("qName", "tName", "pident", "matches", "misMatches", "gapopen",
            "qStart", "qEnd", "tStart", "tEnd", "evalue", "bitscore", "qLen", "tLen", "gaps")

system(paste0("makeblastdb -in reference.fa -parse_seqids -dbtype nucl ",
              " -out blast_db"))

save.contigs = DNAStringSet()
for (i in 1:length(samples)){

  contigs = scanFa(paste0(samples[i]))

  #Matches samples to loci
  system(paste0("blastn -task dc-megablast -db blast_db",
                " -query ", samples[i], " -out ", samples[i], "_match.txt",
                " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen gaps\" ",
                " -num_threads ", threads))

  #Need to load in transcriptome for each species and take the matching transcripts to the database
  match.data = fread(paste0(samples[i], "_match.txt"), sep = "\t", header = F, stringsAsFactors = FALSE)

  if (nrow(match.data) == 0){
    print("No matches found.")
    system(paste0("rm ", samples[i], "_match.txt"))
    next
  }#end if

  setnames(match.data, headers)

  #Deals with duplicate matches to same marker
  marker.names = unique(match.data$tName)
  filt.data = c()
  for (j in 1:length(marker.names)){
    mark.data = match.data[match.data$tName %in% marker.names[j],]

    if(nrow(mark.data) == 0){ next }
    if(nrow(mark.data) != 1){ mark.data = mark.data[mark.data$bitscore == max(mark.data$bitscore)][1] }#end if

    filt.data = rbind(filt.data, mark.data)
  }#end j loop

  #Deals with duplicate matches to same contig
  con.names = names(contigs)
  sample.contigs = DNAStringSet()
  for (j in 1:length(con.names)){

    con.data = filt.data[filt.data$qName %in% con.names[j],]
    temp.contig = contigs[names(contigs) %in% con.names[j]]

    if(nrow(con.data) == 0){ next }
    if(nrow(con.data) == 1){
      names(temp.contig) = paste0(gsub(".fa$", "", samples[i]), "_|_", con.data$tName)
      sample.contigs = append(sample.contigs, temp.contig)
      next
      }#end 1 if

    if (length(unique(con.data$qName)) == length(unique(con.data$tName))){
      con.data = con.data[1,]
      names(temp.contig) = paste0(gsub(".fa$", "", samples[i]), "_|_", con.data$tName)
      sample.contigs = append(sample.contigs, temp.contig)
      next
    }#end 1 if dup

    if (length(unique(con.data$qName)) != length(unique(con.data$tName))){
      for (k in 1:nrow(con.data)){
        cut.contig = subseq(temp.contig, start = con.data$qStart[k], end = con.data$qEnd[k])
        names(cut.contig) = paste0(gsub(".fa$", "", samples[i]), "_|_", con.data$tName[k])
        sample.contigs = append(sample.contigs, cut.contig)
      }#end k loop
      next
    }#end multiple match if

  }#end j loop

  #Saves stuff
  save.contigs = append(save.contigs, sample.contigs)
  system(paste0("rm ", samples[i], "_match.txt"))

  #Writes the loci
  write.loci = as.list(as.character(sample.contigs))
  write.fasta(sequences = write.loci, names = names(write.loci),
              paste0("sampleMarkers/",samples[i]), nbchar = 1000000, as.string = T)

}#end loop

dir.create("markerFastas")
marker.names = unique(gsub(".*_\\|_", "", names(save.contigs)))
for (i in 1:length(marker.names)){

  marker.contigs = save.contigs[grep(marker.names[i], names(save.contigs))]
  names(marker.contigs) = gsub("_\\|_.*", "", names(marker.contigs))
  #Writes the loci
  write.loci = as.list(as.character(marker.contigs))
  write.fasta(sequences = write.loci, names = names(write.loci),
              paste0("markerFastas/", marker.names[i], ".fa"), nbchar = 1000000, as.string = T)
}#end i loop


#Writes the loci
write.loci = as.list(as.character(save.contigs))
write.fasta(sequences = write.loci, names = names(write.loci),
            "sample_blast_pieces.fa", nbchar = 1000000, as.string = T)

#Blast assembled data using remote blast
system(paste0("blastn -task blastn -db nt -query sample_blast_pieces.fa ",
             "-remote -out sample_blast.out -max_target_seqs 8 -max_hsps 5 ",
               "-outfmt \"6 qseqid sseqid sacc pident length mismatch evalue bitscore qlen slen stitle sscinames staxids qcovs\""))

###############################################################################
###############################################################################
###############################################################################
###############################################################################

##Concateantes all the samples together for an easier blast run
setwd(paste0(work.dir, "/", out.dir))

#Gets samples that succeeded
samples = list.files(".")
samples = samples[grep(".fa$", samples)]
samples = samples[samples != "reference.fa"]
samples = samples[samples != "sample_blast_pieces.fa"]

#headers
headers = c("qseqid", "sseqid", "sacc", "pident", "length", "mismatch", "evalue",
            "bitscore", "qlen", "slen","stitle","sscinames", "staxids", "qcovs")

#Need to load in transcriptome for each species and take the matching transcripts to the database
match.data = fread("sample_blast.out", sep = "\t", header = F, stringsAsFactors = FALSE)
setnames(match.data, headers)

header.data = c("Sample", "GenBank_ID", "GenBank_Order", "GenBank_Family", "GenBank_Species",
                "pident", "length", "mismatch", "evalue", "bitscore", "coverage")
#Gets the new tree fiels it made
collect.data = data.table(matrix(as.numeric(0), nrow = length(samples)*1000, ncol = length(header.data)))
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

  #Filters to the best matches
  best.data = sample.data[sample.data$bitscore == max(sample.data$bitscore),]
  best.data = best.data[best.data$pident == max(best.data$pident),]
  best.data = best.data[best.data$evalue == min(best.data$evalue),]
  best.data = best.data[best.data$qcovs == max(best.data$qcovs),]

#  if (nrow(best.data) != 1){stop( "here") }

  for (j in 1:nrow(best.data)){
    #Looks up genbank taxonomy somehow.
    system(paste0("curl https://taxonomy.jgi-psf.org/accession/", best.data$sacc[j],
                             " > output.txt"))
    tax = readLines("output.txt")
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
    set(collect.data, i = as.integer(index.val), j = match("GenBank_ID", header.data), value =  best.data$sacc[j])
    set(collect.data, i = as.integer(index.val), j = match("GenBank_Order", header.data), value =  ord.line)
    set(collect.data, i = as.integer(index.val), j = match("GenBank_Species", header.data), value =  spp.line)
    set(collect.data, i = as.integer(index.val), j = match("GenBank_Family", header.data), value =  fam.line)
    set(collect.data, i = as.integer(index.val), j = match("pident", header.data), value =  best.data$pident[j])
    set(collect.data, i = as.integer(index.val), j = match("length", header.data), value =  best.data$length[j])
    set(collect.data, i = as.integer(index.val), j = match("mismatch", header.data), value =  best.data$mismatch[j])
    set(collect.data, i = as.integer(index.val), j = match("evalue", header.data), value =  best.data$evalue[j])
    set(collect.data, i = as.integer(index.val), j = match("bitscore", header.data), value =  best.data$bitscore[j])
    set(collect.data, i = as.integer(index.val), j = match("coverage", header.data), value =  best.data$qcovs[j])
    index.val = index.val + as.integer(1)
  }#end j loop
}#end i loop

all.data = collect.data[collect.data$Sample != 0,]
samples = unique(all.data$Sample)

spp.data = c()
for (i in 1:length(samples)){

  #Removes duplicate species from the same genbank code (spp data)
  tmp.data = all.data[all.data$Sample %in% samples[i],]
  tmp.data = tmp.data[!duplicated(tmp.data$GenBank_Species),]
  spp.data = rbind(spp.data, tmp.data)

}#end i loop

#Saves the two datasets
write.csv(all.data, "All-top-match_Samples2Genbank.csv", row.names = F)
write.csv(spp.data, "Reduced-top-match_Samples2GenBank.csv", row.names = F)


# END SCRIPT

