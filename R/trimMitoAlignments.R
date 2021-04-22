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

trimMitoAlignments = function(alignment.dir = "Alignments/untrimmed-alignments",
                              alignment.format = "phylip",
                              output.dir = "Alignments/trimmed-alignments",
                              output.format = "phylip",
                              sample.similiarity = TRUE,
                              TrimAl = TRUE,
                              TrimAl.path = "trimal",
                              trim.external = TRUE,
                              min.external.percent = 50,
                              trim.coverage = TRUE,
                              min.coverage.percent = 50,
                              trim.column = TRUE,
                              min.column.gap.percent = 100,
                              alignment.assess = TRUE,
                              min.sample.bp = 0,
                              min.align.length = 0,
                              min.taxa.count = 0,
                              min.gap.percent = 0,
                              overwrite = FALSE) {


  # alignment.dir = "Alignments/untrimmed-alignments"
  # output.dir = "Alignments/trimmed-alignments"
  # alignment.format = "phylip"
  # output.format = "phylip"
  # TrimAl = TRUE
  # trim.column = TRUE
  # alignment.assess = TRUE
  # trim.external = TRUE
  # trim.coverage = TRUE
  # min.coverage.percent = 30
  # min.external.percent = 50
  # min.column.gap.percent = 100
  # overwrite = TRUE
  # min.align.length = 10
  # min.taxa.count = 12
  # min.gap.percent = 50
  # min.sample.bp = 10
  # sample.similiarity = TRUE

  if (alignment.dir == output.dir){ stop("You should not overwrite the original alignments.") }

  if (dir.exists(output.dir) == FALSE) { dir.create(output.dir) }
  if (dir.exists(output.dir) == TRUE) {
    if (overwrite == TRUE){
      system(paste0("rm -r ", output.dir))
      dir.create(output.dir)
    }
  }#end dir exists

  #Gathers alignments
  align.files = list.files(alignment.dir)

  #Data to collect
  header.data = c("Alignment", "Pass", "startSamples", "simSamples", "trimalSamples",
                  "edgeSamples", "covSamples", "columnSamples",
                  "startLength", "simLength", "trimalLength",
                  "edgeLength", "covLength", "columnLength",
                  "startBasepairs", "simBasepairs", "trimalBasepairs",
                  "edgeBasepairs", "covBasepairs", "columnBasepairs",
                  "startGaps","simGaps", "trimalGaps",
                  "edgeGaps", "covGaps", "columnGaps",
                  "startPerGaps", "simPerGaps", "trimalPerGaps",
                  "edgePerGaps", "covPerGaps", "columnPerGaps")

  save.data = data.table::data.table(matrix(as.double(0), nrow = length(align.files), ncol = length(header.data)))
  data.table::setnames(save.data, header.data)
  save.data[, Alignment:=as.character(Alignment)]
  save.data[, Pass:=as.logical(Pass)]

  #Loops through each alignment
  for (i in 1:length(align.files)){
    print(paste0(align.files[i], " Starting..."))

    #Load in alignments
    if (alignment.format == "phylip"){
      align = Biostrings::readAAMultipleAlignment(file = paste0(alignment.dir, "/", align.files[i]), format = "phylip")
      align = Biostrings::DNAStringSet(align)
      save.name = gsub(".phy$", "", align.files[i])
      save.name = gsub(".phylip$", "", save.name)
    }#end phylip

    if (alignment.format == "fasta"){
      align = Rsamtools::scanFa(Rsamtools::FaFile(paste0(alignment.dir, "/", align.files[i])))   # loads up fasta file
      save.name = gsub(".fa$", "", align.files[i])
      save.name = gsub(".fasta$", "", save.name)
    }#end phylip

    # Runs the functions
    #######
    # #Step 1: Strip Ns
    # non.align = replaceAlignmentCharacter(alignment = align,
    #                                       char.find = "N",
    #                                       char.replace = "-")

    non.align = align
    #Summarize all this, no functoin
    data.table::set(save.data, i = as.integer(i), j = match("Alignment", header.data), value = save.name)
    data.table::set(save.data, i = as.integer(i), j = match("startSamples", header.data), value = length(non.align))
    data.table::set(save.data, i = as.integer(i), j = match("startLength", header.data), value = Biostrings::width(non.align)[1] )
    gap.count = countAlignmentGaps(non.align)
    data.table::set(save.data, i = as.integer(i), j = match("startBasepairs", header.data), value = gap.count[2] - gap.count[1])
    data.table::set(save.data, i = as.integer(i), j = match("startGaps", header.data), value = gap.count[1])
    data.table::set(save.data, i = as.integer(i), j = match("startPerGaps", header.data), value = gap.count[3])

    #Step 2. slice trimming
    if (sample.similiarity == TRUE){

      sample.align = trimSampleSimilarity(alignment = non.align,
                                          similarity.threshold = 0.4,
                                          realign.mafft = TRUE)
      non.align = sample.align
      #Saves stat data
      data.table::set(save.data, i = as.integer(i), j = match("simSamples", header.data), value = length(non.align))
      data.table::set(save.data, i = as.integer(i), j = match("simLength", header.data), value = Biostrings::width(non.align)[1])
      gap.count = countAlignmentGaps(non.align)
      data.table::set(save.data, i = as.integer(i), j = match("simBasepairs", header.data), value = gap.count[2] - gap.count[1])
      data.table::set(save.data, i = as.integer(i), j = match("simGaps", header.data), value = gap.count[1])
      data.table::set(save.data, i = as.integer(i), j = match("simPerGaps", header.data), value = gap.count[3])
    }#end if

    #Step 3. Trimal trimming
    if (TrimAl == TRUE){

      trimal.align = trimTrimal(alignment = non.align,
                                quiet = TRUE)
      non.align = trimal.align

      #Saves stat data
      data.table::set(save.data, i = as.integer(i), j = match("trimalSamples", header.data), value = length(non.align))
      data.table::set(save.data, i = as.integer(i), j = match("trimalLength", header.data), value = Biostrings::width(non.align)[1])
      gap.count = countAlignmentGaps(non.align)
      data.table::set(save.data, i = as.integer(i), j = match("trimalBasepairs", header.data), value = gap.count[2] - gap.count[1])
      data.table::set(save.data, i = as.integer(i), j = match("trimalGaps", header.data), value = gap.count[1])
      data.table::set(save.data, i = as.integer(i), j = match("trimalPerGaps", header.data), value = gap.count[3])
    }#end if

    # Step 4. Edge trimming
    if (trim.external == TRUE){
      #external trimming function
      edge.align = trimExternal(alignment = non.align,
                                min.n.seq = ceiling(length(non.align) * (min.external.percent/100)),
                                codon.trim = F)
      non.align = edge.align
      #Saves stat data
      data.table::set(save.data, i = as.integer(i), j = match("edgeSamples", header.data), value = length(non.align))
      data.table::set(save.data, i = as.integer(i), j = match("edgeLength", header.data), value = Biostrings::width(non.align)[1])
      gap.count = countAlignmentGaps(non.align)
      data.table::set(save.data, i = as.integer(i), j = match("edgeBasepairs", header.data), value = gap.count[2] - gap.count[1])
      data.table::set(save.data, i = as.integer(i), j = match("edgeGaps", header.data), value = gap.count[1])
      data.table::set(save.data, i = as.integer(i), j = match("edgePerGaps", header.data), value = gap.count[3])
    }#end trim external

    #Step 5. Evaluate and cut out each sample
    if (trim.coverage == TRUE){
      #sample coverage function
      cov.align = trimSampleCoverage(alignment = non.align,
                                     min.coverage.percent = min.coverage.percent,
                                     min.sample.bp = min.sample.bp,
                                     relative.width = "sample")
      non.align = cov.align
      #Saves stat data
      data.table::set(save.data, i = as.integer(i), j = match("covSamples", header.data), value = length(non.align))
      data.table::set(save.data, i = as.integer(i), j = match("covLength", header.data), value = Biostrings::width(non.align)[1])
      gap.count = countAlignmentGaps(non.align)
      data.table::set(save.data, i = as.integer(i), j = match("covBasepairs", header.data), value = gap.count[2] - gap.count[1])
      data.table::set(save.data, i = as.integer(i), j = match("covGaps", header.data), value = gap.count[1])
      data.table::set(save.data, i = as.integer(i), j = match("covPerGaps", header.data), value = gap.count[3])
    }#end trim.external

    if (trim.column == TRUE){
      #Trim alignment colums
      col.align = trimAlignmentColumns(alignment = non.align,
                                       min.gap.percent = min.column.gap.percent)
      non.align = col.align
      #Saves stat data
      data.table::set(save.data, i = as.integer(i), j = match("columnSamples", header.data), value = length(non.align))
      data.table::set(save.data, i = as.integer(i), j = match("columnLength", header.data), value = Biostrings::width(non.align)[1])
      gap.count = countAlignmentGaps(non.align)
      data.table::set(save.data, i = as.integer(i), j = match("columnBasepairs", header.data), value = gap.count[2] - gap.count[1])
      data.table::set(save.data, i = as.integer(i), j = match("columnGaps", header.data), value = gap.count[1])
      data.table::set(save.data, i = as.integer(i), j = match("columnPerGaps", header.data), value = gap.count[3])
    }#end trim column.

    if (alignment.assess == TRUE) {
      #Assesses the alignment returning TRUE for pass and FALSE for fail
      test.result = alignmentAssess(alignment = non.align,
                                    min.gap.percent = min.gap.percent,
                                    min.taxa.count = min.taxa.count,
                                    min.align.length = min.align.length)

      data.table::set(save.data, i = as.integer(i), j = match("Pass", header.data), value = test.result)

      if (test.result == FALSE){
        print(paste0(align.files[i], " Failed and was discarded."))
      } else {
        write.temp = strsplit(as.character(non.align), "")
        aligned.set = as.matrix(ape::as.DNAbin(write.temp) )
        #readies for saving
        writePhylip(aligned.set, file= paste0(output.dir, "/", save.name, ".phy"), interleave = F)
      }#end else test result
    } else {
      #If no alignment assessing is done, saves
      write.temp = strsplit(as.character(non.align), "")
      aligned.set = as.matrix(ape::as.DNAbin(write.temp) )
      #readies for saving
      writePhylip(aligned.set, file= paste0(output.dir, "/", save.name, ".phy"), interleave = F)
    }#end else

    print(paste0(align.files[i], " Completed."))

  }#end i loop

  #Print and save summary table
  write.csv(save.data, file = paste0("logs/trimming_sample_stats.csv"), row.names = F)

  #Saves log file of things
  if (file.exists(paste0("logs/trimming_summary.log")) == TRUE){ system(paste0("rm logs/trimming_summary.log")) }
  fileConn = file(paste0("logs/trimming_summary.log"), open = "w")
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
  writeLines(paste0("Trimal:"), fileConn)
  writeLines(paste0("Mean samples removed: ",
                    mean(save.data$startSamples - save.data$trimalSamples)), fileConn)
  writeLines(paste0("Mean alignment length reduction: ",
                    mean(save.data$startLength - save.data$trimalLength)), fileConn)
  writeLines(paste0("Mean basepairs trimmed: ",
                    mean(save.data$startBasepairs - save.data$trimalBasepairs)), fileConn)
  writeLines(paste0("Mean gap change: ",
                    mean(save.data$startGaps - save.data$trimalGaps)), fileConn)
  writeLines(paste0("Mean gap percent change: ",
                    mean(save.data$startPerGaps - save.data$trimalPerGaps)), fileConn)
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

