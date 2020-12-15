###############################################################################
#Load packages
###############

devtools::install_github("chutter/PHYLOCAP")
library(PHYLOCAP)

devtools::install_github("chutter/MitoCap")
library(MitoCap)

threads = 8
memory = 80

##########################################################################################################
#Step 1: Preprocess
##########################################################################################################

### Example usage
work.dir = "/home/c111h652/scratch/MitoGenomes"
read.dir = "/home/c111h652/scratch/MitoGenomes/raw-reads-frogs"
decontamination.path = "/home/c111h652/scratch/Contamination_Genomes"
sample.file = "/home/c111h652/scratch/MitoGenomes/Mitogenome_study.csv"
dropbox.dir = "/Research/3_Sequence-Database/Raw-Reads"
dropbox.tok = "/home/c111h652/dropbox-token.RDS"
rdrop2::drop_auth(rdstoken = dropbox.tok)

#Sets working directory (where all the stuff will be saved)
setwd(work.dir)

#Run download function
dropboxDownload(sample.spreadsheet = sample.spread,
                dropbox.directory = dropbox.dir,
                out.directory = read.dir,
                dropbox.token = NULL,
                overwrite = TRUE)


### Make file rename function, otherwise use read names
## FILE RENAME

removeAdaptors(raw.reads = read.dir,
               output.dir = "adaptor-removed-reads",
               mode = "directory",
               fastp.path = "fastp",
               threads = threads,
               mem = memory,
               resume = FALSE,
               overwrite = TRUE,
               quiet = TRUE)

## remove external contamination
removeContamination(input.reads = "adaptor-removed-reads",
                    output.dir = "decontaminated-reads",
                    decontamination.path = decontamination.path,
                    mode = "directory",
                    map.match = 0.99,
                    read.mapper = "bwa",
                    mapper.path = NULL,
                    bbmap.path = "/usr/local/bin/bbsplit.sh",
                    samtools.path = "/usr/local/bin/samtools",
                    bwa.path = "/usr/local/bin/bwa",
                    threads = threads,
                    mem = memory,
                    resume = FALSE,
                    overwrite = TRUE,
                    quiet = TRUE)

#merge paired end reads
mergePairedEndReads(input.reads = "decontaminated-reads",
                    output.dir = "pe-merged-reads",
                    mode = "directory",
                    fastp.path = "fastp",
                    threads = threads,
                    mem = memory,
                    resume = FALSE,
                    overwrite = TRUE,
                    quiet = TRUE)


##########################################################################################################
#Step 2: Mitogenome assembly
##########################################################################################################
#
# devtools::install_github("chutter/PHYLOCAP")
# library(PHYLOCAP)
#
# gb.file = "Crocidura.gb"
# trnascan.path = "/Users/chutter/miniconda3/bin/tRNAscan-SE"
# work.dir = "/Volumes/Rodents/Mitogenomes"
# raw.dir = "/Users/chutter/Dropbox/Mammals/Philippine_Shrews_Raw_Data"
# setwd(work.dir)
#
# #Iteratively assembles to reference
# mitochondrialCapture(input.reads = "pe-merged-reads",
#                      genbank.file = gb.file,
#                      output.dir = "draftContigs",
#                      min.iterations = 5,
#                      max.iterations = 30,
#                      min.length = 16000,
#                      max.length = 40000,
#                      min.ref.id = 0.8,
#                      memory = 6,
#                      threads = 6,
#                      spades.path = "/usr/local/Spades/bin/spades.py",
#                      bbmap.path = "/usr/local/bin/bbmap.sh",
#                      resume = TRUE,
#                      overwrite = FALSE)
#
# ### TO DO: Add a redo option for the non-continuous ones, using itself as a reference
#
# #Annotates mitochondrial contigs
# #To do: output format gff and others
# #To do 2: save new contigs (instead of old) but with corrected start point so half of something isn't lost
# #Combine build mtgenomes with this step. difficut.
#
# annotateMitoContigs(contig.folder = "draftContigs",
#                     genbank.file = gb.file,
#                     blast.path = "blast",
#                     tRNAscan.path = "tRNAscan-SE",
#                     organism.type = "vertebrate",
#                     overwrite = TRUE,
#                     quiet = TRUE)
#
# #Aligns all the different markers
# markerAlignment(input.folder = "Annotations/sample-markers",
#                 genbank.file = gb.file,
#                 threads = 6,
#                 overwrite = TRUE)
#
# #Trims the alignments to ready for concatenation or gene tree estimation
# trimMitoAlignments(alignment.dir = "Alignments/untrimmed-alignments",
#                    alignment.format = "phylip",
#                    output.dir = "Alignments/trimmed-alignments",
#                    output.format = "phylip",
#                    sample.similiarity = TRUE,
#                    TrimAl = TRUE,
#                    TrimAl.path = "trimal",
#                    trim.external = TRUE,
#                    min.external.percent = 50,
#                    trim.coverage = TRUE,
#                    min.coverage.percent = 30,
#                    trim.column = TRUE,
#                    min.column.gap.percent = 100,
#                    alignment.assess = TRUE,
#                    min.sample.bp = 10,
#                    min.align.length = 0,
#                    min.taxa.count = 12,
#                    min.gap.percent = 50,
#                    overwrite = TRUE)
#
# #Aligns the mitogenomes and outputs summary stats
# alignMitogenomes(alignment.folder = "Alignments/untrimmed-alignments",
#                  genbank.file = gb.file,
#                  draft.contigs = "draftContigs",
#                  output.dir = "Genomes",
#                  dataset.name = "untrimmed",
#                  overwrite = TRUE)
#
# #Aligns the mitogenomes and outputs summary stats
# alignMitogenomes(alignment.folder = "Alignments/trimmed-alignments",
#                  genbank.file = gb.file,
#                  draft.contigs = "draftContigs",
#                  output.dir = "Genomes",
#                  dataset.name = "trimmed",
#                  overwrite = FALSE)
#
# ### Working here
# #Builds and standardizes final mitogenomes
# #Reverse the whole thing if its backwards
# # Orient each one to begin with 1
# # Remove over-assembled ones next
#
# buildMitogenomes(annotation.dir = "Annotations",
#                  alignment.folder = "Alignments/untrimmed-alignments",
#                  genome.alignment = "Genomes/alignments/untrimmed_mitogenome_alignment.phy",
#                  genome.dir = "Genomes",
#                  output.dir = "untrimmed-finished",
#                  overwrite = FALSE)
#
#
#






# trimMitoCodingAlignments = function(alignment.dir = "Alignments/untrimmed-alignments",
#                                     alignment.format = "phylip",
#                                     output.dir = "Alignments/coding-alignments",
#                                     output.format = "phylip",
#                                     sample.similiarity = TRUE,
#                                     trim.external = TRUE,
#                                     min.external.percent = 50,
#                                     trim.coverage = TRUE,
#                                     min.coverage.percent = 50,
#                                     trim.column = TRUE,
#                                     min.column.gap.percent = 100,
#                                     alignment.assess = TRUE,
#                                     min.sample.bp = 0,
#                                     min.align.length = 0,
#                                     min.taxa.count = 0,
#                                     min.gap.percent = 0,
#                                     overwrite = FALSE) {
#
#
#   # alignment.dir = "Alignments/untrimmed-alignments"
#   # output.dir = "Alignments/trimmed-alignments"
#   # alignment.format = "phylip"
#   # output.format = "phylip"
#   # TrimAl = TRUE
#   # trim.column = TRUE
#   # alignment.assess = TRUE
#   # trim.external = TRUE
#   # trim.coverage = TRUE
#   # min.coverage.percent = 30
#   # min.external.percent = 50
#   # min.column.gap.percent = 100
#   # overwrite = TRUE
#   # min.align.length = 10
#   # min.taxa.count = 12
#   # min.gap.percent = 50
#   # min.sample.bp = 10
#   # sample.similiarity = TRUE
#
#   if (alignment.dir == output.dir){ stop("You should not overwrite the original alignments.") }
#
#   if (dir.exists(output.dir) == FALSE) { dir.create(output.dir) }
#
#   #So I don't accidentally delete everything while testing resume
#   if (overwrite == TRUE){
#     overwrite = FALSE
#     stop("Error: resume = T and overwrite = T, cannot resume if you are going to delete everything!")
#   }
#
#   if (dir.exists(output.dir) == TRUE) {
#     if (overwrite == TRUE){
#       system(paste0("rm -r ", output.dir))
#       dir.create(output.dir)
#     }
#   }#end dir exists
#
#   #Gathers alignments
#   align.files = list.files(alignment.dir)
#
#   #Data to collect
#   header.data = c("Alignment", "Pass", "startSamples", "simSamples", "trimalSamples",
#                   "edgeSamples", "covSamples", "columnSamples",
#                   "startLength", "simLength", "trimalLength",
#                   "edgeLength", "covLength", "columnLength",
#                   "startBasepairs", "simBasepairs", "trimalBasepairs",
#                   "edgeBasepairs", "covBasepairs", "columnBasepairs",
#                   "startGaps","simGaps", "trimalGaps",
#                   "edgeGaps", "covGaps", "columnGaps",
#                   "startPerGaps", "simPerGaps", "trimalPerGaps",
#                   "edgePerGaps", "covPerGaps", "columnPerGaps")
#
#   save.data = data.table::data.table(matrix(as.double(0), nrow = length(align.files), ncol = length(header.data)))
#   data.table::setnames(save.data, header.data)
#   save.data[, Alignment:=as.character(Alignment)]
#   save.data[, Pass:=as.logical(Pass)]
#
#   #Loops through each alignment
#   for (i in 1:length(align.files)){
#     print(paste0(align.files[i], " Starting..."))
#
#     #Load in alignments
#     if (alignment.format == "phylip"){
#       align = Biostrings::readAAMultipleAlignment(file = paste0(alignment.dir, "/", align.files[i]), format = "phylip")
#       align = Biostrings::DNAStringSet(align)
#       save.name = gsub(".phy$", "", align.files[i])
#       save.name = gsub(".phylip$", "", save.name)
#     }#end phylip
#
#     if (alignment.format == "fasta"){
#       align = Rsamtools::scanFa(Rsamtools::FaFile(paste0(alignment.dir, "/", align.files[i])))   # loads up fasta file
#       save.name = gsub(".fa$", "", align.files[i])
#       save.name = gsub(".fasta$", "", save.name)
#     }#end phylip
#
#     # Runs the functions
#     #######
#     # #Step 1: Strip Ns
#     # non.align = replaceAlignmentCharacter(alignment = align,
#     #                                       char.find = "N",
#     #                                       char.replace = "-")
#
#     non.align = align
#     #Summarize all this, no functoin
#     data.table::set(save.data, i = as.integer(i), j = match("Alignment", header.data), value = save.name)
#     data.table::set(save.data, i = as.integer(i), j = match("startSamples", header.data), value = length(non.align))
#     data.table::set(save.data, i = as.integer(i), j = match("startLength", header.data), value = Biostrings::width(non.align)[1] )
#     gap.count = countAlignmentGaps(non.align)
#     data.table::set(save.data, i = as.integer(i), j = match("startBasepairs", header.data), value = gap.count[2] - gap.count[1])
#     data.table::set(save.data, i = as.integer(i), j = match("startGaps", header.data), value = gap.count[1])
#     data.table::set(save.data, i = as.integer(i), j = match("startPerGaps", header.data), value = gap.count[3])
#
#     #Step 2. slice trimming
#     if (sample.similiarity == TRUE){
#
#       sample.align = trimSampleSimilarity(alignment = non.align,
#                                           similarity.threshold = 0.4,
#                                           realign.mafft = TRUE)
#       non.align = sample.align
#       #Saves stat data
#       data.table::set(save.data, i = as.integer(i), j = match("simSamples", header.data), value = length(non.align))
#       data.table::set(save.data, i = as.integer(i), j = match("simLength", header.data), value = Biostrings::width(non.align)[1])
#       gap.count = countAlignmentGaps(non.align)
#       data.table::set(save.data, i = as.integer(i), j = match("simBasepairs", header.data), value = gap.count[2] - gap.count[1])
#       data.table::set(save.data, i = as.integer(i), j = match("simGaps", header.data), value = gap.count[1])
#       data.table::set(save.data, i = as.integer(i), j = match("simPerGaps", header.data), value = gap.count[3])
#     }#end if
#
#     #Step 3. Trimal trimming
#     if (TrimAl == TRUE){
#
#       trimal.align = trimTrimal(alignment = non.align,
#                                 quiet = TRUE)
#       non.align = trimal.align
#
#       #Saves stat data
#       data.table::set(save.data, i = as.integer(i), j = match("trimalSamples", header.data), value = length(non.align))
#       data.table::set(save.data, i = as.integer(i), j = match("trimalLength", header.data), value = Biostrings::width(non.align)[1])
#       gap.count = countAlignmentGaps(non.align)
#       data.table::set(save.data, i = as.integer(i), j = match("trimalBasepairs", header.data), value = gap.count[2] - gap.count[1])
#       data.table::set(save.data, i = as.integer(i), j = match("trimalGaps", header.data), value = gap.count[1])
#       data.table::set(save.data, i = as.integer(i), j = match("trimalPerGaps", header.data), value = gap.count[3])
#     }#end if
#
#     # Step 4. Edge trimming
#     if (trim.external == TRUE){
#       #external trimming function
#       edge.align = trimExternal(alignment = non.align,
#                                 min.n.seq = ceiling(length(non.align) * (min.external.percent/100)),
#                                 codon.trim = F)
#       non.align = edge.align
#       #Saves stat data
#       data.table::set(save.data, i = as.integer(i), j = match("edgeSamples", header.data), value = length(non.align))
#       data.table::set(save.data, i = as.integer(i), j = match("edgeLength", header.data), value = Biostrings::width(non.align)[1])
#       gap.count = countAlignmentGaps(non.align)
#       data.table::set(save.data, i = as.integer(i), j = match("edgeBasepairs", header.data), value = gap.count[2] - gap.count[1])
#       data.table::set(save.data, i = as.integer(i), j = match("edgeGaps", header.data), value = gap.count[1])
#       data.table::set(save.data, i = as.integer(i), j = match("edgePerGaps", header.data), value = gap.count[3])
#     }#end trim external
#
#     #Step 5. Evaluate and cut out each sample
#     if (trim.coverage == TRUE){
#       #sample coverage function
#       cov.align = trimSampleCoverage(alignment = non.align,
#                                      min.coverage.percent = min.coverage.percent,
#                                      min.sample.bp = min.sample.bp)
#       non.align = cov.align
#       #Saves stat data
#       data.table::set(save.data, i = as.integer(i), j = match("covSamples", header.data), value = length(non.align))
#       data.table::set(save.data, i = as.integer(i), j = match("covLength", header.data), value = Biostrings::width(non.align)[1])
#       gap.count = countAlignmentGaps(non.align)
#       data.table::set(save.data, i = as.integer(i), j = match("covBasepairs", header.data), value = gap.count[2] - gap.count[1])
#       data.table::set(save.data, i = as.integer(i), j = match("covGaps", header.data), value = gap.count[1])
#       data.table::set(save.data, i = as.integer(i), j = match("covPerGaps", header.data), value = gap.count[3])
#     }#end trim.external
#
#     if (trim.column == TRUE){
#       #Trim alignment colums
#       col.align = trimAlignmentColumns(alignment = non.align,
#                                        min.gap.percent = min.column.gap.percent)
#       non.align = col.align
#       #Saves stat data
#       data.table::set(save.data, i = as.integer(i), j = match("columnSamples", header.data), value = length(non.align))
#       data.table::set(save.data, i = as.integer(i), j = match("columnLength", header.data), value = Biostrings::width(non.align)[1])
#       gap.count = countAlignmentGaps(non.align)
#       data.table::set(save.data, i = as.integer(i), j = match("columnBasepairs", header.data), value = gap.count[2] - gap.count[1])
#       data.table::set(save.data, i = as.integer(i), j = match("columnGaps", header.data), value = gap.count[1])
#       data.table::set(save.data, i = as.integer(i), j = match("columnPerGaps", header.data), value = gap.count[3])
#     }#end trim column.
#
#     if (alignment.assess == TRUE) {
#       #Assesses the alignment returning TRUE for pass and FALSE for fail
#       test.result = alignmentAssess(alignment = non.align,
#                                     min.gap.percent = min.gap.percent,
#                                     min.taxa.count = min.taxa.count,
#                                     min.align.length = min.align.length)
#
#       data.table::set(save.data, i = as.integer(i), j = match("Pass", header.data), value = test.result)
#
#       if (test.result == FALSE){
#         print(paste0(align.files[i], " Failed and was discarded."))
#       } else {
#         write.temp = strsplit(as.character(non.align), "")
#         aligned.set = as.matrix(ape::as.DNAbin(write.temp) )
#         #readies for saving
#         writePhylip(aligned.set, file= paste0(output.dir, "/", save.name, ".phy"), interleave = F)
#       }#end else test result
#     } else {
#       #If no alignment assessing is done, saves
#       write.temp = strsplit(as.character(non.align), "")
#       aligned.set = as.matrix(ape::as.DNAbin(write.temp) )
#       #readies for saving
#       writePhylip(aligned.set, file= paste0(output.dir, "/", save.name, ".phy"), interleave = F)
#     }#end else
#
#     print(paste0(align.files[i], " Completed."))
#
#   }#end i loop
#
#   #Print and save summary table
#   write.csv(save.data, file = paste0("logs/trimming_sample_stats.csv"), row.names = F)
#
#   #Saves log file of things
#   if (file.exists(paste0("logs/trimming_summary.log")) == TRUE){ system(paste0("rm logs/trimming_summary.log")) }
#   fileConn = file(paste0("logs/trimming_summary.log"), open = "w")
#   writeLines(paste0("Log file for ", output.dir), fileConn)
#   writeLines(paste0("\n"), fileConn)
#   writeLines(paste0("Overall trimming summary:"), fileConn)
#   writeLines(paste0("------------------------------------------------------------------"), fileConn)
#   writeLines(paste0(""), fileConn)
#   writeLines(paste0("Starting alignments: ", length(align.files)), fileConn)
#   writeLines(paste0("Trimmed alignments: ", length(save.data$Pass[save.data$Pass == TRUE])), fileConn)
#   writeLines(paste0("Discarded alignments: ", length(save.data$Pass[save.data$Pass == FALSE])), fileConn)
#   writeLines(paste0(""), fileConn)
#   writeLines(paste0("Mean samples removed per alignment: ",
#                     mean(save.data$startSamples - save.data$columnSamples)), fileConn)
#   writeLines(paste0("Mean alignment length trimmed per alignment: ",
#                     mean(save.data$startLength - save.data$columnLength)), fileConn)
#   writeLines(paste0("Mean basepairs trimmed per alignment: ",
#                     mean(save.data$startBasepairs - save.data$columnBasepairs)), fileConn)
#   writeLines(paste0("Mean gaps trimmed per alignment: ",
#                     mean(save.data$startGaps - save.data$columnGaps)), fileConn)
#   writeLines(paste0("Mean gap percent trimmed per alignment: ",
#                     mean(save.data$startPerGaps - save.data$columnPerGaps)), fileConn)
#   writeLines(paste0(""), fileConn)
#   writeLines(paste0(""), fileConn)
#   writeLines(paste0("Individual trimming step summary:"), fileConn)
#   writeLines(paste0("------------------------------------------------------------------"), fileConn)
#   writeLines(paste0(""), fileConn)
#   writeLines(paste0("Starting alignments:"), fileConn)
#   writeLines(paste0("Mean samples: ",
#                     mean(save.data$startSamples)), fileConn)
#   writeLines(paste0("Mean alignment length: ",
#                     mean(save.data$startLength)), fileConn)
#   writeLines(paste0("Mean basepairs: ",
#                     mean(save.data$startBasepairs)), fileConn)
#   writeLines(paste0("Mean gaps: ",
#                     mean(save.data$startGaps)), fileConn)
#   writeLines(paste0("Mean gap percentage: ",
#                     mean(save.data$startPerGaps)), fileConn)
#   writeLines(paste0(""), fileConn)
#   writeLines(paste0("Trimal:"), fileConn)
#   writeLines(paste0("Mean samples removed: ",
#                     mean(save.data$startSamples - save.data$trimalSamples)), fileConn)
#   writeLines(paste0("Mean alignment length reduction: ",
#                     mean(save.data$startLength - save.data$trimalLength)), fileConn)
#   writeLines(paste0("Mean basepairs trimmed: ",
#                     mean(save.data$startBasepairs - save.data$trimalBasepairs)), fileConn)
#   writeLines(paste0("Mean gap change: ",
#                     mean(save.data$startGaps - save.data$trimalGaps)), fileConn)
#   writeLines(paste0("Mean gap percent change: ",
#                     mean(save.data$startPerGaps - save.data$trimalPerGaps)), fileConn)
#   writeLines(paste0(""), fileConn)
#   writeLines(paste0("External Trimming:"), fileConn)
#   writeLines(paste0("Mean samples removed: ",
#                     mean(save.data$trimalSamples - save.data$edgeSamples)), fileConn)
#   writeLines(paste0("Mean alignment length reduction: ",
#                     mean(save.data$trimalLength - save.data$edgeLength)), fileConn)
#   writeLines(paste0("Mean basepairs trimmed: ",
#                     mean(save.data$trimalBasepairs - save.data$edgeBasepairs)), fileConn)
#   writeLines(paste0("Mean gap change: ",
#                     mean(save.data$trimalGaps - save.data$edgeGaps)), fileConn)
#   writeLines(paste0("Mean gap percent change: ",
#                     mean(save.data$trimalPerGaps - save.data$edgePerGaps)), fileConn)
#   writeLines(paste0(""), fileConn)
#   writeLines(paste0("Sample Coverage Trimming:"), fileConn)
#   writeLines(paste0("Mean samples removed: ",
#                     mean(save.data$edgeSamples - save.data$covSamples)), fileConn)
#   writeLines(paste0("Mean alignment length reduction: ",
#                     mean(save.data$edgeLength - save.data$covLength)), fileConn)
#   writeLines(paste0("Mean basepairs trimmed: ",
#                     mean(save.data$edgeBasepairs - save.data$covBasepairs)), fileConn)
#   writeLines(paste0("Mean gap change: ",
#                     mean(save.data$edgeGaps - save.data$covGaps)), fileConn)
#   writeLines(paste0("Mean gap percent change: ",
#                     mean(save.data$edgePerGaps - save.data$covPerGaps)), fileConn)
#   writeLines(paste0(""), fileConn)
#   writeLines(paste0("Column Coverage Trimming:"), fileConn)
#   writeLines(paste0("Mean samples removed: ",
#                     mean(save.data$covSamples - save.data$columnSamples)), fileConn)
#   writeLines(paste0("Mean alignment length reduction: ",
#                     mean(save.data$covLength - save.data$columnLength)), fileConn)
#   writeLines(paste0("Mean basepairs trimmed: ",
#                     mean(save.data$covBasepairs - save.data$columnBasepairs)), fileConn)
#   writeLines(paste0("Mean gap change: ",
#                     mean(save.data$covGaps - save.data$columnGaps)), fileConn)
#   writeLines(paste0("Mean gap percent change: ",
#                     mean(save.data$covPerGaps - save.data$columnPerGaps)), fileConn)
#   close(fileConn)
#
# } #end function
#
#



#
#
# trimORF = function(input.folder = NULL,
#                    genbank.file = NULL,
#                    min.taxa = 3,
#                    min.prop.coverage = 0.3,
#                    overwrite = TRUE){
#
#   #Debug
#   input.folder = "Alignments/untrimmed-alignments"
#   genbank.file = "Crocidura.gb"
#   min.taxa = 4
#   min.prop.coverage = 0.3
#   trim.len = 100
#
#   #Create directory and loci to trim
#   dir.create("Alignments/trimmed-alignments")
#   locus.names = list.files("Alignments/untrimmed-alignments/.")
#
#   setwd(work.dir)
#
#   slice.size = 100
#   threshold = 0.50
#   min.cov = 0.5
#
#   #Create directory and loci to trim
#   dir.create("Alignments/coding-alignments")
#
#   locus.names = list.files("Alignments/untrimmed-alignments/.")
#   locus.names = locus.names[grep("-CDS", locus.names, invert = F)]
#
#
#   #Loops through each locus and does operations on them
#   for (i in 1:length(locus.names)){
#
#     ##############
#     #STEP 1: Basic steps
#     ##############
#     #Reads in files
#     align = Biostrings::DNAStringSet(Biostrings::readAAMultipleAlignment(file = paste0("Alignments/untrimmed-alignments/", locus.names[i]), format = "phylip"))
#
#     trimmed = trim.ends(align, min.n.seq = ceiling(length(align) * 0.51), codon.trim = T)
#     save.names = names(trimmed)
#
#     #Gets consensus seq for trimming more
#     con.seq = make.consensus(trimmed, method = "majority")
#
#     #Removes the edge gaps
#     ref.aligned = as.character(con.seq)
#     not.gaps = str_locate_all(ref.aligned, pattern = "[^-]")[[1]][,1]
#
#     #Finds weird gaps to fix
#     temp.gaps = as.numeric(1)
#     for (k in 1:length(not.gaps)-1){ temp.gaps = append(temp.gaps, not.gaps[k+1]-not.gaps[k]) }
#     temp.gaps = temp.gaps-1
#     names(temp.gaps) = not.gaps
#     gap.spots = temp.gaps[temp.gaps %% 3 != 0]
#
#     del.col<-c()
#     if (length(gap.spots) != 0){
#       #Loop through each potential bad thing and fix
#       for (k in 1:length(gap.spots)){
#         del.col = append(del.col, (as.numeric(names(gap.spots[k]))-gap.spots[k]):(as.numeric(names(gap.spots[k]))-1))
#       }
#     }#end if
#
#     ##############
#     #STEP 3: Fixes large gaps at ends of alignment
#     ##############
#     #Looks for gaps that clearly wrong and not 3 BP
#     new.align = strsplit(as.character(trimmed), "")
#     x = as.matrix(as.DNAbin(new.align))
#
#     rem.n<-c()
#     for (k in 1:ncol(x)){
#       gaps<-table(as.character(x[,k]))
#       per.gaps<-gaps[names(gaps) == "-"]/nrow(x)
#
#       if (length(per.gaps) == 0){ next }
#
#       #Records column when the gaps exceed this percentage
#       if (per.gaps >= 0.75){ del.col<-append(del.col, k) }
#
#       #Removes gap columns only consisting of Ns
#       n.gaps<-gaps[names(gaps) != "-"]
#       if (length(n.gaps) == 1){
#         if (names(n.gaps) == "n"){ rem.n<-append(rem.n, k)}
#       }
#
#     }#end k loop
#
#     #combines columns to be deleted
#     fin.del<-c(rem.n, del.col)
#     if (length(fin.del) != 0){
#       x<-x[,-fin.del] }
#     #Removes bad columsn and coverts alignment back to DNASTringSet
#     char.align<-as.list(data.frame(t(as.character(x))))
#     temp.align<-lapply(char.align, FUN = function(x) paste(x, collapse = ""))
#     trimmed<-DNAStringSet(unlist(temp.align))
#     names(trimmed)<-save.names
#
#     ##############
#     #STEP 4: Gathers table of best and longest stop codon free frames for each seq
#     ##############
#     #Checks to make sure the codon position is correct
#     save.frame = data.frame()
#     save.all = data.frame()
#     for (j in 1:length(trimmed)){
#       #Finds open reading frames
#       temp.codon = find.orf(trimmed[j], mitochondrial = T, min.size = 50 )
#       if(nrow(temp.codon) == 0){
#         samp.frame<-cbind(Sample = names(trimmed[j]), FrameStart = 0, FrameEnd = 0, Size = 0, sppSize = 0, Frame = "0")
#         save.frame<-rbind(save.frame, samp.frame)
#         next
#       }
#
#       all.frame<-temp.codon[temp.codon$Size >= max(temp.codon$Size) * .70,]
#       big.frame<-temp.codon[temp.codon$Size == max(temp.codon$Size),]
#
#       if (nrow(big.frame) >= 2){
#         #Picks the best from this order of things
#         temp.stop<-big.frame[big.frame$Frame == "F1",]
#         if (nrow(temp.stop) == 0){ temp.stop<-big.frame[big.frame$Frame == "R1",] }
#         if (nrow(temp.stop) == 0){ temp.stop<-big.frame[big.frame$Frame == "F2",] }
#         if (nrow(temp.stop) == 0){ temp.stop<-big.frame[big.frame$Frame == "R2",] }
#         if (nrow(temp.stop) == 0){ temp.stop<-big.frame[big.frame$Frame == "F3",] }
#         if (nrow(temp.stop) == 0){ temp.stop<-big.frame[big.frame$Frame == "R3",] }
#         big.frame<-temp.stop
#       }
#
#       # #Saves teh data
#       samp.frame<-cbind(Sample = names(trimmed[j]), big.frame)
#       temp.size<-unlist(strsplit(as.character(trimmed[j]), ""), use.names = F)
#
#       #Starts from the beginning and end to fill in end gaps
#       sub.size<-0
#       for (q in 1:length(temp.size)){ if (temp.size[q] == "-"){ sub.size<-sub.size+1 } else { break } }
#       for (q in length(temp.size):1){ if (temp.size[q] == "-"){ sub.size<-sub.size+1 } else { break } }
#
#       #Saves final data
#       samp.frame<-cbind(samp.frame, sppSize = length(temp.size)-sub.size)
#       save.frame<-rbind(save.frame, samp.frame)
#       all.frame<-cbind(Sample = names(trimmed[j]), all.frame)
#       all.frame<-cbind(all.frame, sppSize = length(temp.size)-sub.size)
#       save.all<-rbind(save.all, all.frame)
#     }#end j loop
#
#     #Moves on if there are no frames found. Saves to Anon folder?
#     if (unique(save.frame$Frame)[1] == "0"){
#       print(paste(locus.names[i], " sucked. No frames found.", sep = ""))
#       next
#     }
#
#     ##############
#     #STEP 5: Uses previous data to find a consistent frame
#     ##############
#     #Looks at the overall data rather than the best indiv data to find a consistent frame
#     temp.all = save.all
#     frame.names = unique(temp.all$Frame)
#     #Goes through the equally good frames and reduces to frames with the same range
#     very.best = data.frame()
#     for (k in 1:length(frame.names)){
#       temp.best<-temp.all[temp.all$Frame == frame.names[k],]
#       starts<-table(temp.best$FrameStart)[table(temp.best$FrameStart) == max(table(temp.best$FrameStart))]
#       ends<-table(temp.best$FrameEnd)[table(temp.best$FrameEnd) == max(table(temp.best$FrameEnd))]
#
#       #Removes duplicates
#       starts<-starts[as.numeric(names(starts)) == min(as.numeric(names(starts)))]
#       ends<-ends[as.numeric(names(ends)) == min(as.numeric(names(ends)))]
#
#       super.best<-temp.best[temp.best$FrameStart == as.numeric(names(starts)),]
#       super.best<-super.best[super.best$FrameEnd == as.numeric(names(ends)),]
#       very.best<-rbind(very.best, super.best)
#     }#end k loop
#
#     #Moves on if there are no frames found. Saves to Anon folder?
#     if (nrow(very.best) == 0){
#       print(paste(locus.names[i], " sucked. No frames found.", sep = ""))
#       next
#     }
#
#     ##############
#     #STEP 6: Selects the best frame
#     ##############
#     #Picks out the best frame
#     best.frame = table(very.best$Frame)[table(very.best$Frame) == max(table(very.best$Frame))]
#
#     #If there are multiple good frames pick the biggest
#     if (length(best.frame) != 1){
#       temp.fix<-very.best[very.best$Frame %in% names(best.frame),]
#       bigger<-temp.fix[temp.fix$Size == max(temp.fix$Size),]
#       best.frame<-table(bigger$Frame)[table(bigger$Frame) == max(table(bigger$Frame))]
#     }#end if
#
#     #If they are same size just pick from this order
#     if (length(best.frame) != 1){
#       #Picks the best from this order of things
#       temp.stop<-best.frame[names(best.frame) == "F1"]
#       if (length(temp.stop) == 0){ temp.stop<-best.frame[names(best.frame) == "R1"] }
#       if (length(temp.stop) == 0){ temp.stop<-best.frame[names(best.frame) == "F2"] }
#       if (length(temp.stop) == 0){ temp.stop<-best.frame[names(best.frame) == "R2"] }
#       if (length(temp.stop) == 0){ temp.stop<-best.frame[names(best.frame) == "F3"] }
#       if (length(temp.stop) == 0){ temp.stop<-best.frame[names(best.frame) == "R3"] }
#       best.frame<-temp.stop
#     }
#
#     #Checks the remaining ORF size
#     temp.size<-save.all[save.all$Frame == names(best.frame),]
#     samp.spp<-temp.size[temp.size$sppSize == max(temp.size$sppSize),]
#     samp.seq<-trimmed[names(trimmed) == samp.spp$Sample[1]]
#     samp.size<-nchar(gsub("-", "", as.character(samp.seq)))
#
#     if (mean(temp.size$Size) <= samp.size*.5){ print(paste(locus.names[i], " was small.", sep = "")) }
#
#     if (mean(temp.size$Size) <= samp.size*.25){
#       print(paste(locus.names[i], " sucked. Too few sequence left.", sep = ""))
#       next
#     }
#
#     if (best.frame <= length(trimmed) * .5){
#       print(paste(locus.names[i], " sucked. No cosistent frame.", sep = ""))
#       next
#     }
#
#     #Reverses if it needs to
#     if (length(grep("R", names(best.frame))) != 0){
#       new.align<-reverseComplement(trimmed)
#     }else { new.align<-trimmed }
#
#     ##############
#     #STEP 7: Gets start and stop coordinates for each sequence and find best alignment
#     ##############
#     #Gets trimming locations
#     frame.ranges<-save.all[save.all$Frame == names(best.frame),]
#
#     #Gets potential starts and ends
#     starts<-table(frame.ranges$FrameStart)[table(frame.ranges$FrameStart) == max(table(frame.ranges$FrameStart))]
#     ends<-table(frame.ranges$FrameEnd)[table(frame.ranges$FrameEnd) == max(table(frame.ranges$FrameEnd))]
#     starts<-starts[starts == max(starts)]
#     ends<-ends[ends == max(ends)]
#
#     if (length(starts) != 1 || length(ends) != 1){
#       frame.ranges<-save.all[save.all$Frame == names(best.frame),]
#       starts<-table(frame.ranges$FrameStart)[table(frame.ranges$FrameStart) == max(table(frame.ranges$FrameStart))]
#       ends<-table(frame.ranges$FrameEnd)[table(frame.ranges$FrameEnd) == max(table(frame.ranges$FrameEnd))]
#       starts<-starts[starts == max(starts)]
#       ends<-ends[ends == max(ends)]
#     }
#
#     if (length(starts) != 1 || length(ends) != 1){
#       starts<-starts[as.numeric(names(starts)) == min(as.numeric(names(starts)))]
#       ends<-ends[as.numeric(names(ends)) == max(as.numeric(names(ends)))]
#     }
#
#     ###################
#     #STEP 8: Makes sure entire alignment is a multiple of 3
#     ###################
#     anu.start<-as.numeric(names(starts))
#     new.end<-as.numeric(names(ends))
#     new.len<-new.end-(anu.start-1)
#
#     #Gets a new end to keep in multiples of 3 for proteins
#     if (length(new.len[which(new.len %%3==0)]) == 0) {
#       anu.end<-new.end-1
#     } else { anu.end<-new.end }
#
#     new.len<-anu.end-(anu.start-1)
#     if (length(new.len[which(new.len %%3==0)]) == 0) {
#       anu.end<-new.end-2
#     } else { anu.end<-anu.end }
#
#     #Trims sequence with new coords
#     done.seq<-subseq(start = anu.start, end = anu.end, x = new.align)
#
#     ###################
#     #STEP 9: Trim out odd start/end bases
#     ###################
#     codon.seq<-DNAStringSet()
#     for (k in 1:length(done.seq)){
#       ref.aligned<-as.character(done.seq[k])
#
#       #Chcecks at beginning of sequence
#       not.gaps<-str_locate_all(ref.aligned, pattern = "[^-]")[[1]][,1]
#
#       if (length(not.gaps) <= min.len){ next }
#
#       #Checks if its odd, delete 1 base
#       if ( (not.gaps[1]-1) %%3 == 2){ substr(ref.aligned, not.gaps[1], not.gaps[1])<-"-" }
#
#       #Deletes 2 bases its off by
#       if ( (not.gaps[1]-1) %%3 == 1){
#         substr(ref.aligned, not.gaps[1], not.gaps[1])<-"-"
#         substr(ref.aligned, not.gaps[1]+1, not.gaps[1]+1)<-"-"
#       }#end if
#
#       #checks for end of sequence
#       not.gaps<-str_locate_all(ref.aligned, pattern = "[^-]")[[1]][,1]
#       #counts characters
#       char.len<-(not.gaps[length(not.gaps)]-not.gaps[1])+1
#       end.pos<-not.gaps[length(not.gaps)]
#       #removes odd characters at ends
#       if ( char.len %%3 == 1){ substr(ref.aligned, end.pos, end.pos)<-"-" }
#
#       if ( char.len %%3 == 2){
#         substr(ref.aligned, end.pos-1, end.pos-1)<-"-"
#         substr(ref.aligned, end.pos, end.pos)<-"-"
#       } #end if
#
#       #Saves final seqs
#       save.seq<-DNAStringSet(ref.aligned)
#       names(save.seq)<-names(done.seq[k])
#       codon.seq<-append(codon.seq, save.seq)
#
#     } #END K
#
#     ###################
#     #STEP 10: Change stop codons to N
#     ###################
#     #Finds stop codons to replace
#     n.seq = DNAStringSet(gsub("-", "N", as.character(codon.seq)))
#     stop.seq = DNAStringSet()
#     del.taxa = c()
#     for (k in 1:length(n.seq)){
#       stop.data = find.codon(n.seq[k], genetic.code = "mitochondrial", reverse = T)
#       stop.data = stop.data[stop.data$frame == "F1",]
#       #If there is a stop codon take bold action
#       if (stop.data$start[1] == 0){
#         stop.seq = append(stop.seq, n.seq[k])
#       } else {
#         #Goes through each codon
#         #  for (y in 1:nrow(stop.data)){
#         #Saves final seqs
#         #    substr(ref.aligned, stop.data$start[y], stop.data$start[y])<-"N"
#         #    substr(ref.aligned, stop.data$start[y]+1, stop.data$start[y]+1)<-"N"
#         #    substr(ref.aligned, stop.data$start[y]+2, stop.data$start[y]+2)<-"N"
#         #  }#end Y LOOP
#         del.taxa = append(del.taxa, names(ref.aligned) )
#         #Saves final data
#         #save.seq = DNAStringSet(ref.aligned)
#         #names(save.seq) = names(n.seq[k])
#         #stop.seq = append(stop.seq, save.seq)
#       }#end if state
#     }# END K loop
#
#     if (length(stop.seq) <= min.taxa){ next }
#
#     ###################
#     #FINAL STEP: Save everything after some final spp and length filtering
#     ###################
#     #Removes sequences that are less than a certain coverage\
#     t.align<-strsplit(as.character(stop.seq), "")
#     len.loci<-lapply(t.align, function (x) x[x != "N"])
#     spp.len<-unlist(lapply(len.loci, function (x) length(x)))
#     spp.rem<-spp.len[spp.len <= max(spp.len) * as.numeric(min.cov)]
#     spp.rem<-append(spp.rem, spp.len[spp.len <= as.numeric(min.len)]  )
#
#     if (length(spp.rem) > 0){
#       red.align<-t.align[!names(t.align) %in% unique(names(spp.rem))]
#     } else { red.align<-t.align }
#
#     #Removes if all removed
#     if (length(red.align) == 0){
#       print(paste0("not enough alignment remains for ", locus.names[i]))
#       next
#     }
#
#     #removes loci with too few taxa
#     if (length(names(red.align)) <= as.numeric(min.taxa)){
#       system(paste("rm ", input.file, sep = ""))
#       print(paste(input.file, "deleted. Too few taxa after trimming."))
#       next
#     }
#
#     #writes alignment
#     mat.align<-lapply(red.align, tolower)
#     write.align<-as.matrix(as.DNAbin(mat.align))
#
#     #readies for saving
#     write.phy(write.align, file=paste0("Alignments/coding-alignments/", locus.names[i]), interleave = F)
#
#     #readies for saving
#     write.phy(write.align, file=paste0("Alignments/trimmed-alignments/", locus.names[i]), interleave = F)
#
#
#   }#end i loop final #######
#
#   # *** trimming trims out good alignment that should be included in the mtgenome.
#   # have an alternative where it saves slice trim above?
#
#   # END SCRIPT
#
# }
#


#numpt detecting?

#But first, need to match raw reads again to the mtgenomes
#Call SNPS, figure out

#Another weird thing, look at when pseudogenes introgress into the nuc genome
#Ie some are closer maybe because they introgressed to the ref species a long time ago

#Stats at this point:
#Percent complete
#Number of numpts?
#Number genes
#BP
#Ns
# Circular = Yes/No

