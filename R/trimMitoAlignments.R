#' @title trimMitoAlignments
#'
#' @description Function for batch trimming a folder of alignments, with the various trimming functions available to select from
#'
#' @param alignment.dir path to a folder of sequence alignments in phylip format.
#'
#' @param alignment.format available input alignment formats: fasta or phylip
#'
#' @param output.dir contigs are added into existing alignment if algorithm is "add"
#'
#' @param output.format available output formats: phylip
#'
#' @param sample.similarity remove samples too divergent from consensus, values 0-1 for proportion similar sites
#'
#' @param TrimAl if a file name is provided, save.name will be used to save aligment to file as a fasta
#'
#' @param TrimAl.path path to a folder of sequence alignments in phylip format.
#'
#' @param trim.external give a save name if you wnat to save the summary to file.
#'
#' @param min.external.percent TRUE to supress mafft screen output
#'
#' @param trim.coverage path to a folder of sequence alignments in phylip format.
#'
#' @param min.coverage.percent contigs are added into existing alignment if algorithm is "add"
#'
#' @param trim.column algorithm to use: "add" add sequences with "add.contigs"; "localpair" for local pair align. All others available
#'
#' @param min.column.gap.percent TRUE applies the adjust sequence direction function of MAFFT
#'
#' @param alignment.assess if a file name is provided, save.name will be used to save aligment to file as a fasta
#'
#' @param min.sample.bp path to a folder of sequence alignments in phylip format.
#'
#' @param min.align.length give a save name if you wnat to save the summary to file.
#'
#' @param min.taxa.count TRUE to supress mafft screen output
#'
#' @param min.gap.percent if a file name is provided, save.name will be used to save aligment to file as a fasta
#'
#' @param overwrite TRUE to supress mafft screen output
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

trimMitoAlignments = function(alignment.dir = NULL,
                              alignment.format = "phylip",
                              output.dir = NULL,
                              output.format = "phylip",
                              TAPER = FALSE,
                              TAPER.path = NULL,
                              julia.path = NULL,
                              TrimAl = FALSE,
                              TrimAl.path = NULL,
                              trim.external = TRUE,
                              min.external.percent = 50,
                              trim.coverage = TRUE,
                              min.coverage.percent = 50,
                              trim.column = TRUE,
                              min.column.gap.percent = 100,
                              convert.ambiguous.sites = FALSE,
                              alignment.assess = TRUE,
                              min.coverage.bp = 0,
                              min.alignment.length = 0,
                              min.taxa.alignment = 0,
                              max.alignment.gap.percent = 0,
                              threads = 1,
                              memory = 1,
                              overwrite = FALSE) {

  # #devtools::install_github("chutter/PHYLOCAP", upgrade = "never")
  # library(PhyloCap)
  # library(foreach)
  #
  #
  # #Save directory
  # #work.dir = "/Volumes/Rodents/Murinae/Trimming"
  # #align.dir = "/Volumes/Rodents/Murinae/Trimming/genes-untrimmed"
  # #feat.gene.names = "/Volumes/Rodents/Murinae/Selected_Transcripts/Mus_gene_metadata.csv"
  # #work.dir = "/home/c111h652/scratch/Rodents/Trimming"
  #
  # work.dir = "/home/c111h652/scratch/Rodents/Trimming"
  # align.dir = "/home/c111h652/scratch/Rodents/Trimming/Emily/genes_untrimmed"
  #
  # work.dir = "/Volumes/Rodents/Murinae/Trimming"
  # align.dir =  "/Volumes/Rodents/Murinae/Trimming/Emily/genes_untrimmed"
  # out.name = "Emily"
  #
  #
  # setwd(work.dir)
  # alignment.dir = align.dir
  # output.dir = paste0(out.name, "/genes_trimmed2")
  # alignment.format = "phylip"
  # output.format = "phylip"
  # TAPER = TRUE
  # TAPER.path = "/usr/local/bin"
  # julia.path = "/Applications/Julia-1.6.app/Contents/Resources/julia/bin"
  # #TAPER.path = "/home/c111h652/programs"
  # #julia.path = "/programs/julia/bin"
  # TrimAl = TRUE
  # TrimAl.path = "/Users/chutter/miniconda3/bin"
  # trim.external = TRUE
  # min.external.percent = 50
  # trim.coverage = TRUE
  # min.coverage.percent = 35
  # trim.column = TRUE
  # min.column.gap.percent = 50
  # convert.ambiguous.sites = FALSE
  # alignment.assess = TRUE
  # min.coverage.bp = 60
  # min.alignment.length = 100
  # min.taxa.alignment = 5
  # max.alignment.gap.percent = 50
  # threads = 2
  # memory = 6
  # overwrite = FALSE
  # resume = TRUE

  if (is.null(TrimAl.path) == FALSE){
    b.string = unlist(strsplit(TrimAl.path, ""))
    if (b.string[length(b.string)] != "/") {
      TrimAl.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { TrimAl.path = NULL }

  if (alignment.dir == output.dir){ stop("You should not overwrite the original alignments.") }

  # if (dir.exists(output.dir) == FALSE) { dir.create(output.dir) }

  if (dir.exists(output.dir) == TRUE) {
    if (overwrite == TRUE){
      system(paste0("rm -r ", output.dir))
      dir.create(output.dir)
    }
  } else { dir.create(output.dir) }

  #Gathers alignments
  align.files = list.files(alignment.dir)

  if (length(align.files) == 0) { stop("alignment files could not be found.") }

  #Skips files done already if resume = TRUE
  if (overwrite == FALSE){
    done.files = list.files(output.dir)
    align.files = align.files[!gsub("\\..*", "", align.files) %in% gsub("\\..*", "", done.files)]
  }

  if (length(align.files) == 0) { stop("All alignments have already been completed and overwrite = FALSE.") }

  #Data to collect
  header.data = c("Alignment", "Pass", "startSamples", "trimalSamples",
                  "edgeSamples", "columnSamples", "covSamples",
                  "startLength",  "trimalLength",
                  "edgeLength", "columnLength", "covLength",
                  "startBasepairs", "tapirBasepairs", "trimalBasepairs",
                  "edgeBasepairs", "columnBasepairs", "covBasepairs",
                  "startGaps", "trimalGaps",
                  "edgeGaps", "columnGaps", "covGaps",
                  "startPerGaps", "trimalPerGaps",
                  "edgePerGaps", "columnPerGaps", "covPerGaps")

  save.data = data.table::data.table(matrix(as.double(0), nrow = length(align.files), ncol = length(header.data)))
  data.table::setnames(save.data, header.data)
  save.data[, Alignment:=as.character(Alignment)]
  save.data[, Pass:=as.logical(Pass)]

  #Sets up multiprocessing
  #cl = parallel::makeCluster(threads, outfile = "")
  #doParallel::registerDoParallel(cl)
  #mem.cl = floor(memory/threads)

  #Loops through each locus and does operations on them
  #save.data = foreach::foreach(i=1:length(align.files), .combine = rbind, .packages = c("PhyloCap", "foreach", "Biostrings","Rsamtools", "ape", "stringr", "data.table")) %dopar% {
  for (i in 1:length(align.files)){
    print(paste0(align.files[i], " Starting..."))

    #Load in alignments
    if (alignment.format == "phylip"){
      align = Biostrings::readAAMultipleAlignment(file = paste0(alignment.dir, "/", align.files[i]), format = "phylip")

      #  align = Biostrings::readDNAStringSet(file = paste0(alignment.dir, "/", align.files[i]), format = "phylip")
      #  align = readLines(paste0(alignment.dir, "/", align.files[i]))[-1]
      #  align = gsub(".*\\ ", "", align)
      #  char.count = nchar(align)

      align = Biostrings::DNAStringSet(align)
      save.name = gsub(".phy$", "", align.files[i])
      save.name = gsub(".phylip$", "", save.name)
    }#end phylip

    if (alignment.format == "fasta"){
      align = Biostrings::readDNAStringSet(paste0(alignment.dir, "/", align.files[i]) )
      save.name = gsub(".fa$", "", align.files[i])
      save.name = gsub(".fasta$", "", save.name)
    }#end phylip

    # Runs the functions
    #######
    #Step 1: Strip Ns
    non.align = PhyloCap::replaceAlignmentCharacter(alignment = align,
                                          char.find = "N",
                                          char.replace = "-")

    #Convert ambiguous sites
    if (convert.ambiguous.sites == TRUE){
      non.align = PhyloCap::convertAmbiguousConsensus(alignment = non.align)
    }#end phylip


    #Summarize all this, no functoin
    data.table::set(save.data, i = as.integer(1), j = match("Alignment", header.data), value = save.name)
    data.table::set(save.data, i = as.integer(1), j = match("startSamples", header.data), value = length(non.align))
    data.table::set(save.data, i = as.integer(1), j = match("startLength", header.data), value = Biostrings::width(non.align)[1] )
    gap.count = PhyloCap::countAlignmentGaps(non.align)
    data.table::set(save.data, i = as.integer(1), j = match("startBasepairs", header.data), value = gap.count[2] - gap.count[1])
    data.table::set(save.data, i = as.integer(1), j = match("startGaps", header.data), value = gap.count[1])
    data.table::set(save.data, i = as.integer(1), j = match("startPerGaps", header.data), value = gap.count[3])

    #Step 3. Trimal trimming
    if (TrimAl == TRUE && length(non.align) != 0){

      trimal.align = PhyloCap::trimTrimal(alignment = non.align,
                                trimal.path = TrimAl.path,
                                quiet = TRUE)
      non.align = trimal.align
      #Saves stat data
      data.table::set(save.data, i = as.integer(1), j = match("trimalSamples", header.data), value = length(trimal.align))
      data.table::set(save.data, i = as.integer(1), j = match("trimalLength", header.data), value = Biostrings::width(trimal.align)[1])
      gap.count = PhyloCap::countAlignmentGaps(non.align)
      data.table::set(save.data, i = as.integer(1), j = match("trimalBasepairs", header.data), value = gap.count[2] - gap.count[1])
      data.table::set(save.data, i = as.integer(1), j = match("trimalGaps", header.data), value = gap.count[1])
      data.table::set(save.data, i = as.integer(1), j = match("trimalPerGaps", header.data), value = gap.count[3])
    }#end if

    # Step 4. Edge trimming
    if (trim.external == TRUE && length(non.align) != 0){
      #external trimming function
      edge.align = PhyloCap::trimExternal(alignment = non.align,
                                min.n.seq = ceiling(length(non.align) * (min.external.percent/100)),
                                codon.trim = F)

      if (class(edge.align) == "numeric") { edge.align = DNAStringSet() }
      non.align = edge.align

      #Saves stat data
      data.table::set(save.data, i = as.integer(1), j = match("edgeSamples", header.data), value = length(edge.align))
      data.table::set(save.data, i = as.integer(1), j = match("edgeLength", header.data), value = Biostrings::width(edge.align)[1])
      gap.count = PhyloCap::countAlignmentGaps(non.align)
      data.table::set(save.data, i = as.integer(1), j = match("edgeBasepairs", header.data), value = gap.count[2] - gap.count[1])
      data.table::set(save.data, i = as.integer(1), j = match("edgeGaps", header.data), value = gap.count[1])
      data.table::set(save.data, i = as.integer(1), j = match("edgePerGaps", header.data), value = gap.count[3])
    }#end trim external

    if (trim.column == TRUE && length(non.align) != 0){
      #Trim alignment colums
      col.align = PhyloCap::trimAlignmentColumns(alignment = non.align,
                                       min.gap.percent = min.column.gap.percent)
      non.align = col.align
      #Saves stat data
      data.table::set(save.data, i = as.integer(1), j = match("columnSamples", header.data), value = length(col.align))
      data.table::set(save.data, i = as.integer(1), j = match("columnLength", header.data), value = Biostrings::width(col.align)[1])
      gap.count = PhyloCap::countAlignmentGaps(non.align)
      data.table::set(save.data, i = as.integer(1), j = match("columnBasepairs", header.data), value = gap.count[2] - gap.count[1])
      data.table::set(save.data, i = as.integer(1), j = match("columnGaps", header.data), value = gap.count[1])
      data.table::set(save.data, i = as.integer(1), j = match("columnPerGaps", header.data), value = gap.count[3])
    }#end trim column.

    #Step 5. Evaluate and cut out each sample
    if (trim.coverage == TRUE && length(non.align) != 0){
      #sample coverage function
      cov.align = PhyloCap::trimSampleCoverage(alignment = non.align,
                                     min.coverage.percent = min.coverage.percent,
                                     min.coverage.bp = min.coverage.bp,
                                     relative.width = "sample")

      non.align = cov.align
      #Saves stat data
      data.table::set(save.data, i = as.integer(1), j = match("covSamples", header.data), value = length(cov.align))
      data.table::set(save.data, i = as.integer(1), j = match("covLength", header.data), value = Biostrings::width(cov.align)[1])
      gap.count = PhyloCap::countAlignmentGaps(non.align)
      data.table::set(save.data, i = as.integer(1), j = match("covBasepairs", header.data), value = gap.count[2] - gap.count[1])
      data.table::set(save.data, i = as.integer(1), j = match("covGaps", header.data), value = gap.count[1])
      data.table::set(save.data, i = as.integer(1), j = match("covPerGaps", header.data), value = gap.count[3])
    }#end trim.external

    #Step 6
    if (alignment.assess == TRUE) {
      #Assesses the alignment returning TRUE for pass and FALSE for fail
      test.result = PhyloCap::alignmentAssess(alignment = non.align,
                                    max.alignment.gap.percent = max.alignment.gap.percent,
                                    min.taxa.alignment = min.taxa.alignment,
                                    min.alignment.length = min.alignment.length)

      data.table::set(save.data, i = as.integer(1), j = match("Pass", header.data), value = test.result)

      if (test.result == FALSE){
        print(paste0(align.files[i], " Failed and was discarded."))
      } else {
        write.temp = strsplit(as.character(non.align), "")
        aligned.set = as.matrix(ape::as.DNAbin(write.temp) )
        #readies for saving
        PhyloCap::writePhylip(aligned.set, file= paste0(output.dir, "/", save.name, ".phy"), interleave = F)
      }#end else test result
    } else {
      #If no alignment assessing is done, saves
      write.temp = strsplit(as.character(non.align), "")
      aligned.set = as.matrix(ape::as.DNAbin(write.temp) )
      #readies for saving
      PhyloCap::writePhylip(aligned.set, file= paste0(output.dir, "/", save.name, ".phy"), interleave = F)
    }#end else

    print(paste0(align.files[i], " Completed."))

    #rm()
    #gc()

  }#end i loop

 # parallel::stopCluster(cl)

  #Print and save summary table
  write.csv(save.data, file = paste0(output.dir, "_trimming-summary.csv"), row.names = F)
  #Saves log file of things
  if (file.exists(paste0(output.dir, ".log")) == TRUE){ system(paste0("rm ", output.dir, ".log")) }
  fileConn = file(paste0(output.dir, ".log"), open = "w")
  writeLines(paste0("Log file for ", output.dir), fileConn)
  writeLines(paste0("\n"), fileConn)
  writeLines(paste0("Overall trimming summary:"), fileConn)
  writeLines(paste0("------------------------------------------------------------------"), fileConn)
  writeLines(paste0(""), fileConn)
  writeLines(paste0("Starting alignments: ", length(align.files)), fileConn)
  writeLines(paste0("Trimmed alignments: ", length(save.data$Pass[save.data$Pass == TRUE])), fileConn)
  writeLines(paste0("Discarded alignments: ", length(save.data$Pass[save.data$Pass == FALSE])), fileConn)
  writeLines(paste0(""), fileConn)
  writeLines(paste0("Mean samples removed per alignment: ",
                    mean(save.data$startSamples - save.data$columnSamples)), fileConn)
  writeLines(paste0("Mean alignment length trimmed per alignment: ",
                    mean(save.data$startLength - save.data$columnLength)), fileConn)
  writeLines(paste0("Mean basepairs trimmed per alignment: ",
                    mean(save.data$startBasepairs - save.data$columnBasepairs)), fileConn)
  writeLines(paste0("Mean gaps trimmed per alignment: ",
                    mean(save.data$startGaps - save.data$columnGaps)), fileConn)
  writeLines(paste0("Mean gap percent trimmed per alignment: ",
                    mean(save.data$startPerGaps - save.data$columnPerGaps)), fileConn)
  writeLines(paste0(""), fileConn)
  writeLines(paste0(""), fileConn)
  writeLines(paste0("Individual trimming step summary:"), fileConn)
  writeLines(paste0("------------------------------------------------------------------"), fileConn)
  writeLines(paste0(""), fileConn)
  writeLines(paste0("Starting alignments:"), fileConn)
  writeLines(paste0("Mean samples: ",
                    mean(save.data$startSamples)), fileConn)
  writeLines(paste0("Mean alignment length: ",
                    mean(save.data$startLength)), fileConn)
  writeLines(paste0("Mean basepairs: ",
                    mean(save.data$startBasepairs)), fileConn)
  writeLines(paste0("Mean gaps: ",
                    mean(save.data$startGaps)), fileConn)
  writeLines(paste0("Mean gap percentage: ",
                    mean(save.data$startPerGaps)), fileConn)
  writeLines(paste0(""), fileConn)
  writeLines(paste0("hmmCleaner:"), fileConn)
  writeLines(paste0("Mean samples removed: ",
                    mean(save.data$startSamples - save.data$hmmSamples)), fileConn)
  writeLines(paste0("Mean alignment length reduction: ",
                    mean(save.data$startLength - save.data$hmmLength)), fileConn)
  writeLines(paste0("Mean basepairs trimmed: ",
                    mean(save.data$startBasepairs - save.data$hmmBasepairs)), fileConn)
  writeLines(paste0("Mean gap change: ",
                    mean(save.data$startGaps - save.data$hmmGaps)), fileConn)
  writeLines(paste0("Mean gap percent change: ",
                    mean(save.data$startPerGaps - save.data$hmmPerGaps)), fileConn)
  writeLines(paste0(""), fileConn)
  writeLines(paste0("Trimal:"), fileConn)
  writeLines(paste0("Mean samples removed: ",
                    mean(save.data$hmmSamples - save.data$trimalSamples)), fileConn)
  writeLines(paste0("Mean alignment length reduction: ",
                    mean(save.data$hmmLength - save.data$trimalLength)), fileConn)
  writeLines(paste0("Mean basepairs trimmed: ",
                    mean(save.data$hmmBasepairs - save.data$trimalBasepairs)), fileConn)
  writeLines(paste0("Mean gap change: ",
                    mean(save.data$hmmGaps - save.data$trimalGaps)), fileConn)
  writeLines(paste0("Mean gap percent change: ",
                    mean(save.data$hmmPerGaps - save.data$trimalPerGaps)), fileConn)
  writeLines(paste0(""), fileConn)
  writeLines(paste0("External Trimming:"), fileConn)
  writeLines(paste0("Mean samples removed: ",
                    mean(save.data$trimalSamples - save.data$edgeSamples)), fileConn)
  writeLines(paste0("Mean alignment length reduction: ",
                    mean(save.data$trimalLength - save.data$edgeLength)), fileConn)
  writeLines(paste0("Mean basepairs trimmed: ",
                    mean(save.data$trimalBasepairs - save.data$edgeBasepairs)), fileConn)
  writeLines(paste0("Mean gap change: ",
                    mean(save.data$trimalGaps - save.data$edgeGaps)), fileConn)
  writeLines(paste0("Mean gap percent change: ",
                    mean(save.data$trimalPerGaps - save.data$edgePerGaps)), fileConn)
  writeLines(paste0(""), fileConn)
  writeLines(paste0("Sample Coverage Trimming:"), fileConn)
  writeLines(paste0("Mean samples removed: ",
                    mean(save.data$edgeSamples - save.data$covSamples)), fileConn)
  writeLines(paste0("Mean alignment length reduction: ",
                    mean(save.data$edgeLength - save.data$covLength)), fileConn)
  writeLines(paste0("Mean basepairs trimmed: ",
                    mean(save.data$edgeBasepairs - save.data$covBasepairs)), fileConn)
  writeLines(paste0("Mean gap change: ",
                    mean(save.data$edgeGaps - save.data$covGaps)), fileConn)
  writeLines(paste0("Mean gap percent change: ",
                    mean(save.data$edgePerGaps - save.data$covPerGaps)), fileConn)
  writeLines(paste0(""), fileConn)
  writeLines(paste0("Column Coverage Trimming:"), fileConn)
  writeLines(paste0("Mean samples removed: ",
                    mean(save.data$covSamples - save.data$columnSamples)), fileConn)
  writeLines(paste0("Mean alignment length reduction: ",
                    mean(save.data$covLength - save.data$columnLength)), fileConn)
  writeLines(paste0("Mean basepairs trimmed: ",
                    mean(save.data$covBasepairs - save.data$columnBasepairs)), fileConn)
  writeLines(paste0("Mean gap change: ",
                    mean(save.data$covGaps - save.data$columnGaps)), fileConn)
  writeLines(paste0("Mean gap percent change: ",
                    mean(save.data$covPerGaps - save.data$columnPerGaps)), fileConn)
  close(fileConn)

} #end function
