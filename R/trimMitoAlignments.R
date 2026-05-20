#' @title trimMitoAlignments
#'
#' @description Batch-trims a folder of mitochondrial marker alignments using
#'   a pipeline of optional steps: N-to-gap conversion, TrimAl automated
#'   trimming, edge (external gap) trimming, column gap trimming, and sample
#'   coverage trimming. Each alignment is optionally assessed for minimum
#'   quality thresholds before being written to the output directory in phylip
#'   format. A CSV summary of trimming statistics and a plain-text log file are
#'   also written.
#'
#' @param alignment.dir path to a folder of untrimmed alignments (phylip or
#'   fasta format).
#'
#' @param alignment.format format of input alignments: \code{"phylip"} or
#'   \code{"fasta"}.
#'
#' @param output.dir path to the output directory for trimmed alignments.
#'   Must differ from \code{alignment.dir}.
#'
#' @param output.format output alignment format; currently only
#'   \code{"phylip"} is written.
#'
#' @param TAPER logical; if TRUE, run TAPER masking (requires \code{julia} and
#'   TAPER to be installed).
#'
#' @param TAPER.path path to the directory containing the TAPER script.
#'
#' @param julia.path path to the directory containing \code{julia}.
#'
#' @param TrimAl logical; if TRUE, run TrimAl automated column trimming.
#'
#' @param TrimAl.path path to the directory containing \code{trimal}.
#'
#' @param trim.external logical; if TRUE, trim edge (leading/trailing) gaps
#'   from the alignment.
#'
#' @param min.external.percent minimum proportion of samples (0--100) that must
#'   have sequence at a column before it is kept during edge trimming.
#'
#' @param trim.coverage logical; if TRUE, remove samples with insufficient
#'   sequence coverage.
#'
#' @param min.coverage.percent minimum percentage of alignment columns a sample
#'   must cover to be retained.
#'
#' @param trim.column logical; if TRUE, remove alignment columns with too many
#'   gaps.
#'
#' @param min.column.gap.percent maximum proportion of gaps (0--100) allowed in
#'   a column before it is removed.
#'
#' @param convert.ambiguous.sites logical; if TRUE, convert ambiguous IUPAC
#'   bases to the consensus character.
#'
#' @param alignment.assess logical; if TRUE, apply final quality assessment
#'   filters and discard failing alignments.
#'
#' @param min.coverage.bp minimum number of unambiguous base pairs a sample
#'   must have to be retained (passed to \code{trimSampleCoverage}).
#'
#' @param min.alignment.length minimum alignment length (columns) required for
#'   an alignment to pass assessment.
#'
#' @param min.taxa.alignment minimum number of taxa required for an alignment
#'   to pass assessment.
#'
#' @param max.alignment.gap.percent maximum overall gap percentage allowed in
#'   an alignment for it to pass assessment.
#'
#' @param threads number of CPU threads (reserved for future parallel use).
#'
#' @param memory amount of RAM in GB (reserved for future parallel use).
#'
#' @param overwrite logical; if TRUE, already-trimmed files in \code{output.dir}
#'   are re-trimmed; if FALSE, existing output files are skipped.
#'
#' @return Invisibly returns NULL. Writes trimmed phylip alignments to
#'   \code{output.dir}, a trimming summary CSV, and a plain-text log file.
#'
#' @examples
#' \dontrun{
#' trimMitoAlignments(
#'   alignment.dir          = "Alignments/untrimmed-alignments",
#'   output.dir             = "Alignments/trimmed-alignments",
#'   alignment.format       = "phylip",
#'   TrimAl                 = TRUE,
#'   TrimAl.path            = "/path/to/trimal/bin",
#'   trim.external          = TRUE,
#'   min.external.percent   = 50,
#'   trim.coverage          = TRUE,
#'   min.coverage.percent   = 50,
#'   trim.column            = TRUE,
#'   min.column.gap.percent = 100,
#'   alignment.assess       = TRUE,
#'   min.alignment.length   = 100,
#'   min.taxa.alignment     = 4,
#'   overwrite              = FALSE
#' )
#' }
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

  # #devtools::install_github("chutter/PhyloProcessR", upgrade = "never")
  # library(PhyloProcessR)
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
      unlink(output.dir, recursive = TRUE)
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
  #save.data = foreach::foreach(i=1:length(align.files), .combine = rbind, .packages = c("PhyloProcessR", "foreach", "Biostrings","Rsamtools", "ape", "stringr", "data.table")) %dopar% {
  for (i in seq_along(align.files)){
    message(align.files[i], " Starting...")

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
    }#end fasta

    # Runs the functions
    #######
    #Step 1: Strip Ns
    non.align = PhyloProcessR::replaceAlignmentCharacter(alignment = align,
                                          char.find = "N",
                                          char.replace = "-")

    #Convert ambiguous sites
    if (convert.ambiguous.sites == TRUE){
      non.align = PhyloProcessR::convertAmbiguousConsensus(alignment = non.align)
    }#end phylip


    #Summarize all this, no functoin
    data.table::set(save.data, i = as.integer(i), j = match("Alignment", header.data), value = save.name)
    data.table::set(save.data, i = as.integer(i), j = match("startSamples", header.data), value = length(non.align))
    data.table::set(save.data, i = as.integer(i), j = match("startLength", header.data), value = Biostrings::width(non.align)[1] )
    gap.count = PhyloProcessR::countAlignmentGaps(non.align)
    data.table::set(save.data, i = as.integer(i), j = match("startBasepairs", header.data), value = gap.count[2] - gap.count[1])
    data.table::set(save.data, i = as.integer(i), j = match("startGaps", header.data), value = gap.count[1])
    data.table::set(save.data, i = as.integer(i), j = match("startPerGaps", header.data), value = gap.count[3])

    #Step 3. Trimal trimming
    if (TrimAl == TRUE && length(non.align) != 0){

      trimal.align = PhyloProcessR::trimTrimal(alignment = non.align,
                                trimal.path = TrimAl.path,
                                quiet = TRUE)
      non.align = trimal.align
      #Saves stat data
      data.table::set(save.data, i = as.integer(i), j = match("trimalSamples", header.data), value = length(trimal.align))
      data.table::set(save.data, i = as.integer(i), j = match("trimalLength", header.data), value = Biostrings::width(trimal.align)[1])
      gap.count = PhyloProcessR::countAlignmentGaps(non.align)
      data.table::set(save.data, i = as.integer(i), j = match("trimalBasepairs", header.data), value = gap.count[2] - gap.count[1])
      data.table::set(save.data, i = as.integer(i), j = match("trimalGaps", header.data), value = gap.count[1])
      data.table::set(save.data, i = as.integer(i), j = match("trimalPerGaps", header.data), value = gap.count[3])
    }#end if

    # Step 4. Edge trimming
    if (trim.external == TRUE && length(non.align) != 0){
      #external trimming function
      edge.align = PhyloProcessR::trimExternal(alignment = non.align,
                                min.n.seq = ceiling(length(non.align) * (min.external.percent/100)),
                                codon.trim = FALSE)

      if (is.numeric(edge.align)) { edge.align = Biostrings::DNAStringSet() }
      non.align = edge.align

      #Saves stat data
      data.table::set(save.data, i = as.integer(i), j = match("edgeSamples", header.data), value = length(edge.align))
      data.table::set(save.data, i = as.integer(i), j = match("edgeLength", header.data), value = Biostrings::width(edge.align)[1])
      gap.count = PhyloProcessR::countAlignmentGaps(non.align)
      data.table::set(save.data, i = as.integer(i), j = match("edgeBasepairs", header.data), value = gap.count[2] - gap.count[1])
      data.table::set(save.data, i = as.integer(i), j = match("edgeGaps", header.data), value = gap.count[1])
      data.table::set(save.data, i = as.integer(i), j = match("edgePerGaps", header.data), value = gap.count[3])
    }#end trim external

    if (trim.column == TRUE && length(non.align) != 0){
      #Trim alignment colums
      col.align = PhyloProcessR::trimAlignmentColumns(alignment = non.align,
                                       min.gap.percent = min.column.gap.percent)
      non.align = col.align
      #Saves stat data
      data.table::set(save.data, i = as.integer(i), j = match("columnSamples", header.data), value = length(col.align))
      data.table::set(save.data, i = as.integer(i), j = match("columnLength", header.data), value = Biostrings::width(col.align)[1])
      gap.count = PhyloProcessR::countAlignmentGaps(non.align)
      data.table::set(save.data, i = as.integer(i), j = match("columnBasepairs", header.data), value = gap.count[2] - gap.count[1])
      data.table::set(save.data, i = as.integer(i), j = match("columnGaps", header.data), value = gap.count[1])
      data.table::set(save.data, i = as.integer(i), j = match("columnPerGaps", header.data), value = gap.count[3])
    }#end trim column.

    #Step 5. Evaluate and cut out each sample
    if (trim.coverage == TRUE && length(non.align) != 0){
      #sample coverage function
      cov.align = PhyloProcessR::trimSampleCoverage(alignment = non.align,
                                     min.coverage.percent = min.coverage.percent,
                                     min.coverage.bp = min.coverage.bp,
                                     relative.width = "sample")

      non.align = cov.align
      #Saves stat data
      data.table::set(save.data, i = as.integer(i), j = match("covSamples", header.data), value = length(cov.align))
      data.table::set(save.data, i = as.integer(i), j = match("covLength", header.data), value = Biostrings::width(cov.align)[1])
      gap.count = PhyloProcessR::countAlignmentGaps(non.align)
      data.table::set(save.data, i = as.integer(i), j = match("covBasepairs", header.data), value = gap.count[2] - gap.count[1])
      data.table::set(save.data, i = as.integer(i), j = match("covGaps", header.data), value = gap.count[1])
      data.table::set(save.data, i = as.integer(i), j = match("covPerGaps", header.data), value = gap.count[3])
    }#end trim.external

    #Step 6
    if (alignment.assess == TRUE) {
      #Assesses the alignment returning TRUE for pass and FALSE for fail
      test.result = PhyloProcessR::alignmentAssess(alignment = non.align,
                                    max.alignment.gap.percent = max.alignment.gap.percent,
                                    min.taxa.alignment = min.taxa.alignment,
                                    min.alignment.length = min.alignment.length)

      data.table::set(save.data, i = as.integer(i), j = match("Pass", header.data), value = test.result)

      if (test.result == FALSE){
        message(align.files[i], " Failed and was discarded.")
      } else {
        write.temp = strsplit(as.character(non.align), "")
        aligned.set = as.matrix(ape::as.DNAbin(write.temp) )
        #readies for saving
        PhyloProcessR::writePhylip(aligned.set, file= paste0(output.dir, "/", save.name, ".phy"), interleave = FALSE)
      }#end else test result
    } else {
      #If no alignment assessing is done, saves
      write.temp = strsplit(as.character(non.align), "")
      aligned.set = as.matrix(ape::as.DNAbin(write.temp) )
      #readies for saving
      PhyloProcessR::writePhylip(aligned.set, file= paste0(output.dir, "/", save.name, ".phy"), interleave = FALSE)
    }#end else

    message(align.files[i], " Completed.")

    #rm()
    #gc()

  }#end i loop

 # parallel::stopCluster(cl)

  #Print and save summary table
  write.csv(save.data, file = paste0(output.dir, "_trimming-summary.csv"), row.names = FALSE)
  #Saves log file of things
  if (file.exists(paste0(output.dir, ".log")) == TRUE){ file.remove(paste0(output.dir, ".log")) }
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
                    mean(save.data$startSamples - save.data$covSamples)), fileConn)
  writeLines(paste0("Mean alignment length trimmed per alignment: ",
                    mean(save.data$startLength - save.data$covLength)), fileConn)
  writeLines(paste0("Mean basepairs trimmed per alignment: ",
                    mean(save.data$startBasepairs - save.data$covBasepairs)), fileConn)
  writeLines(paste0("Mean gaps trimmed per alignment: ",
                    mean(save.data$startGaps - save.data$covGaps)), fileConn)
  writeLines(paste0("Mean gap percent trimmed per alignment: ",
                    mean(save.data$startPerGaps - save.data$covPerGaps)), fileConn)
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
  writeLines(paste0("Column Trimming:"), fileConn)
  writeLines(paste0("Mean samples removed: ",
                    mean(save.data$edgeSamples - save.data$columnSamples)), fileConn)
  writeLines(paste0("Mean alignment length reduction: ",
                    mean(save.data$edgeLength - save.data$columnLength)), fileConn)
  writeLines(paste0("Mean basepairs trimmed: ",
                    mean(save.data$edgeBasepairs - save.data$columnBasepairs)), fileConn)
  writeLines(paste0("Mean gap change: ",
                    mean(save.data$edgeGaps - save.data$columnGaps)), fileConn)
  writeLines(paste0("Mean gap percent change: ",
                    mean(save.data$edgePerGaps - save.data$columnPerGaps)), fileConn)
  writeLines(paste0(""), fileConn)
  writeLines(paste0("Sample Coverage Trimming:"), fileConn)
  writeLines(paste0("Mean samples removed: ",
                    mean(save.data$columnSamples - save.data$covSamples)), fileConn)
  writeLines(paste0("Mean alignment length reduction: ",
                    mean(save.data$columnLength - save.data$covLength)), fileConn)
  writeLines(paste0("Mean basepairs trimmed: ",
                    mean(save.data$columnBasepairs - save.data$covBasepairs)), fileConn)
  writeLines(paste0("Mean gap change: ",
                    mean(save.data$columnGaps - save.data$covGaps)), fileConn)
  writeLines(paste0("Mean gap percent change: ",
                    mean(save.data$columnPerGaps - save.data$covPerGaps)), fileConn)
  close(fileConn)

  return(invisible(NULL))

} #end function
