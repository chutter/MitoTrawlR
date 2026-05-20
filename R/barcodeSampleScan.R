#' @title barcodeSampleScan
#'
#' @description Assembles a target barcode marker (e.g. 16S rRNA, COI) from
#'   cleaned reads for each sample using iterative BBMap read recruitment
#'   followed by SPAdes / CAP3 assembly, then identifies the best-matching
#'   species via BLAST. The iterative approach starts permissive (min.ref.id =
#'   0.70) and accumulates reads across rounds, assembling progressively better
#'   sequences even from low-coverage or divergent samples. Assembly stops as
#'   soon as the barcode region length is covered (target.length = reference
#'   barcode length) so the function does not attempt to build a full genome.
#'
#'   A per-sample summary is appended to logs/barcodeSampleScan_summary.csv
#'   after every sample so the file grows safely across successive runs.
#'
#' @param input.reads path to a directory of cleaned reads. Samples may be in
#'   per-sample sub-directories or identified by a shared filename prefix.
#'
#' @param barcode.fasta path to a FASTA file of target barcode reference
#'   sequence(s) (e.g. a single 16S or COI representative). Used as the
#'   initial assembly seed and to filter off-target contigs during assembly.
#'
#' @param output.directory path where per-sample assemblies and BLAST results
#'   are written. Default: \code{"barcodeScan"}.
#'
#' @param database.fasta path to a FASTA file of named barcode sequences for
#'   local BLAST identification. \code{NULL} (the default) queries the NCBI
#'   \code{nt} database remotely — no local file needed, but an internet
#'   connection is required and queries will be slower.
#'
#' @param hits.per.sample number of top BLAST hits to keep per sample.
#'   Default: \code{5}.
#'
#' @param per.max.length fraction above the reference barcode length that
#'   triggers the max-length guard during assembly. Default: \code{0.50}.
#'
#' @param min.iterations minimum iterative assembly rounds before convergence
#'   is tested. Default: \code{3}.
#'
#' @param max.iterations maximum iterative assembly rounds. Default: \code{10}.
#'
#' @param min.ref.id starting BBMap minimum identity for read recruitment.
#'   Permissive on the first round (0.70), auto-tightened to 0.95 thereafter.
#'   Default: \code{0.70}.
#'
#' @param spades.path system path to the directory containing
#'   \code{spades.py}; NULL searches the system PATH.
#'
#' @param bbmap.path system path to the directory containing \code{bbmap.sh};
#'   NULL searches the system PATH.
#'
#' @param blast.path system path to the directory containing \code{blastn} and
#'   \code{makeblastdb}; NULL searches the system PATH.
#'
#' @param cap3.path system path to the directory containing \code{cap3}; NULL
#'   searches the system PATH.
#'
#' @param memory RAM in GB to allocate to BBMap and SPAdes. Default: \code{1}.
#'
#' @param threads CPU threads for BBMap, SPAdes, and BLAST. Default: \code{1}.
#'
#' @param overwrite logical; if TRUE the output directory is deleted and
#'   recreated. Default: \code{FALSE}.
#'
#' @param quiet logical; if TRUE tool screen output is suppressed.
#'   Default: \code{TRUE}.
#'
#' @return invisibly; per-sample FASTA assemblies and BLAST result tables are
#'   written to output.directory, and a cross-sample summary is appended to
#'   logs/barcodeSampleScan_summary.csv.
#'
#' @export

barcodeSampleScan = function(input.reads = NULL,
                             barcode.fasta = NULL,
                             output.directory = "barcodeScan",
                             database.fasta = NULL,
                             hits.per.sample = 5,
                             per.max.length = 0.50,
                             min.iterations = 3,
                             max.iterations = 10,
                             min.ref.id = 0.70,
                             spades.path = NULL,
                             bbmap.path = NULL,
                             blast.path = NULL,
                             cap3.path = NULL,
                             memory = 1,
                             threads = 1,
                             overwrite = FALSE,
                             quiet = TRUE) {

  # #Debug
  # setwd("/Volumes/LaCie/Brygomantis")
  # input.reads = "processed-reads/cleaned-reads"
  # barcode.fasta = "barcodeScan/16S_barcode.fa"
  # output.directory = "barcodeScan"
  # database.fasta = NULL   # NULL = remote NCBI nt
  # memory = 36
  # threads = 8
  # overwrite = TRUE
  # spades.path = "/Users/chutter/Bioinformatics/conda-envs/PhyloCap/bin"
  # bbmap.path  = "/Users/chutter/Bioinformatics/conda-envs/PhyloCap/bin"
  # blast.path  = "/Users/chutter/Bioinformatics/conda-envs/PhyloCap/bin"
  # cap3.path   = "/Users/chutter/Bioinformatics/conda-envs/PhyloCap/bin"
  # hits.per.sample = 5
  # per.max.length  = 0.50
  # min.iterations  = 3
  # max.iterations  = 10
  # min.ref.id      = 0.70

  #################################################
  ### Tool path normalisation
  #################################################

  if (!is.null(spades.path)) {
    b.string = unlist(strsplit(spades.path, ""))
    if (b.string[length(b.string)] != "/") {
      spades.path = paste0(append(b.string, "/"), collapse = "")
    }
  } else { spades.path = "" }

  if (!is.null(bbmap.path)) {
    b.string = unlist(strsplit(bbmap.path, ""))
    if (b.string[length(b.string)] != "/") {
      bbmap.path = paste0(append(b.string, "/"), collapse = "")
    }
  } else { bbmap.path = "" }

  if (!is.null(blast.path)) {
    b.string = unlist(strsplit(blast.path, ""))
    if (b.string[length(b.string)] != "/") {
      blast.path = paste0(append(b.string, "/"), collapse = "")
    }
  } else { blast.path = "" }

  if (!is.null(cap3.path)) {
    b.string = unlist(strsplit(cap3.path, ""))
    if (b.string[length(b.string)] != "/") {
      cap3.path = paste0(append(b.string, "/"), collapse = "")
    }
  } else { cap3.path = "" }

  #################################################
  ### Quick checks
  #################################################

  if (is.null(input.reads))   { stop("Please provide input reads.") }
  if (is.null(barcode.fasta)) { stop("Please provide a barcode reference FASTA.") }
  if (!file.exists(barcode.fasta)) { stop("Barcode reference FASTA not found: ", barcode.fasta) }

  use.remote.blast = is.null(database.fasta)
  if (!use.remote.blast && !file.exists(database.fasta)) {
    stop("Barcode database FASTA not found: ", database.fasta)
  }

  #################################################
  ### Output and log directories
  #################################################

  if (!dir.exists(output.directory)) {
    dir.create(output.directory)
  } else if (overwrite) {
    unlink(output.directory, recursive = TRUE)
    dir.create(output.directory)
  }

  if (!dir.exists("logs")) { dir.create("logs") }

  dir.create(paste0(output.directory, "/sample-barcodes"),  showWarnings = FALSE)
  dir.create(paste0(output.directory, "/blast-reference"),  showWarnings = FALSE)

  #################################################
  ### One-time setup: reference blast db + identification db
  #################################################

  # Copy barcode reference and build BLAST db used internally by iterativeAssemble
  # for off-target filtering, and for the final per-sample identification BLAST.
  ref.fa = paste0(output.directory, "/blast-reference/reference.fa")
  if (!file.exists(ref.fa)) {
    system(paste0("cp ", barcode.fasta, " ", ref.fa))
  }
  ref.db = paste0(output.directory, "/blast-reference/barcode")
  if (!file.exists(paste0(ref.db, ".nhr"))) {
    system(paste0(blast.path, "makeblastdb -in ", ref.fa,
                  " -parse_seqids -dbtype nucl -out ", ref.db),
           ignore.stdout = quiet, ignore.stderr = quiet)
  }

  # Read reference to get expected barcode length (used as early-exit target
  # so iterativeAssemble stops once the barcode region is covered, rather than
  # continuing to assemble a full genome).
  ref.seq = Biostrings::readDNAStringSet(ref.fa)
  ref.len  = max(Biostrings::width(ref.seq))
  min.len  = ref.len
  max.len  = ref.len + (ref.len * per.max.length)

  # Local identification database (built once if database.fasta is provided)
  id.blast.db = NULL
  if (!use.remote.blast) {
    id.blast.dir = paste0(output.directory, "/barcode-database")
    dir.create(id.blast.dir, showWarnings = FALSE)
    id.db.fa = paste0(id.blast.dir, "/database.fa")
    if (!file.exists(id.db.fa)) {
      system(paste0("cp ", database.fasta, " ", id.db.fa))
      system(paste0(blast.path, "makeblastdb -in ", id.db.fa,
                    " -parse_seqids -dbtype nucl -out ", id.blast.dir, "/database"),
             ignore.stdout = quiet, ignore.stderr = quiet)
    }
    id.blast.db = paste0(id.blast.dir, "/database")
  } else {
    cat(" Barcode BLAST mode: remote NCBI nt (internet connection required)\n")
  }

  #################################################
  ### Sample discovery and resume
  #################################################

  files   = list.files(path = input.reads, full.names = TRUE, recursive = TRUE)
  reads   = files[grep(pattern = "fastq|fq", x = files)]

  # Prefer sub-directory-per-sample layout; fall back to flat files
  samples = list.dirs(input.reads, recursive = FALSE, full.names = FALSE)
  if (length(samples) == 0) {
    flat = list.files(input.reads, recursive = FALSE, full.names = FALSE)
    samples = unique(gsub("_L00.*|_R[12][._].*|_READ[123][._].*|\\.fastq.*|\\.fq.*", "", flat))
    samples = samples[nchar(samples) > 0]
  }

  if (!overwrite) {
    done = list.files(paste0(output.directory, "/sample-barcodes"), full.names = FALSE)
    done = gsub("\\.fa$", "", done)
    samples = samples[!samples %in% done]
  }

  if (length(samples) == 0) { return("No samples remain to analyze.") }

  blast.headers = c("qName", "tName", "pident", "matches", "misMatches", "gapopen",
                    "qStart", "qEnd", "tStart", "tEnd", "evalue", "bitscore",
                    "qLen", "tLen", "gaps")

  # Helper: append a result row to the rolling per-sample summary CSV
  append.summary = function(row) {
    csv.path = "logs/barcodeSampleScan_summary.csv"
    if (file.exists(csv.path)) {
      existing = read.csv(csv.path, stringsAsFactors = FALSE)
      existing = existing[!existing$Sample %in% row$Sample, ]
      row = rbind(existing, row)
    }
    write.csv(row, file = csv.path, row.names = FALSE)
  }

  #################################################
  ### Main loop — one sample at a time
  #################################################

  for (i in seq_along(samples)) {

    cat("\n [barcodeSampleScan] sample", i, "of", length(samples), ":", samples[i], "\n")

    #------------------------------------------------------
    # Identify and concatenate reads across lanes
    #------------------------------------------------------
    sample.reads = reads[grep(samples[i], reads)]

    read1.reads = sample.reads[grep("_1\\.f.*|-1\\.f.*|_R1_.*|-R1_.*|_R1-.*|-R1-.*|READ1.*|_R1\\.fast.*|-R1\\.fast.*", sample.reads)]
    read2.reads = sample.reads[grep("_2\\.f.*|-2\\.f.*|_R2_.*|-R2_.*|_R2-.*|-R2-.*|READ2.*|_R2\\.fast.*|-R2\\.fast.*", sample.reads)]
    read3.reads = sample.reads[grep("_3\\.f.*|-3\\.f.*|_R3_.*|-R3_.*|READ3.*|_singleton.*|-singleton.*", sample.reads)]

    if (length(read1.reads) == 0) {
      warning(samples[i], " read1 files not found. Skipping.")
      next
    }

    it.sample.reads = c()

    if (length(read1.reads) >= 2) {
      combined.r1 = paste0(input.reads, "/", samples[i], "_ALL_READ1.fastq.gz")
      system(paste0("cat ", paste(read1.reads, collapse = " "), " > ", combined.r1))
      it.sample.reads[1] = combined.r1
    } else {
      it.sample.reads[1] = read1.reads
    }

    if (length(read2.reads) >= 2) {
      combined.r2 = paste0(input.reads, "/", samples[i], "_ALL_READ2.fastq.gz")
      system(paste0("cat ", paste(read2.reads, collapse = " "), " > ", combined.r2))
      it.sample.reads[2] = combined.r2
    } else if (length(read2.reads) == 1) {
      it.sample.reads[2] = read2.reads
    }

    if (length(read3.reads) >= 2) {
      combined.r3 = paste0(input.reads, "/", samples[i], "_ALL_READ3.fastq.gz")
      system(paste0("cat ", paste(read3.reads, collapse = " "), " > ", combined.r3))
      it.sample.reads[3] = combined.r3
    } else if (length(read3.reads) == 1) {
      it.sample.reads[3] = read3.reads
    }

    it.sample.reads = it.sample.reads[!is.na(it.sample.reads)]

    #------------------------------------------------------
    # Iterative BBMap + SPAdes / CAP3 assembly.
    # target.length = ref.len so the loop exits as soon as
    # the barcode region is covered — no full-genome assembly.
    #------------------------------------------------------
    dir.create(paste0(output.directory, "/sample-barcodes/", samples[i]),
               showWarnings = FALSE)

    barcode.contigs = iterativeAssemble(
      input.reads    = it.sample.reads,
      reference      = ref.fa,
      mapper         = "bbmap",
      min.ref.id     = min.ref.id,
      memory         = memory,
      threads        = threads,
      min.iterations = min.iterations,
      max.iterations = max.iterations,
      min.length     = min.len,
      max.length     = max.len,
      target.length  = ref.len,
      spades.path    = spades.path,
      bbmap.path     = bbmap.path,
      blast.path     = blast.path,
      cap3.path      = cap3.path,
      quiet          = quiet
    )

    # Delete any multi-lane concatenated files created above
    combined.files = grep("_ALL_READ", c(it.sample.reads), value = TRUE)
    if (length(combined.files) > 0) { unlink(combined.files) }

    if (length(barcode.contigs) == 0) {
      cat(" ", samples[i], ": no barcode assembled — no reads matching reference.\n")
      temp.row = data.frame(Sample = samples[i],
                            ContigLength = 0L, ContigCount = 0L,
                            BestMatch = "No-match",
                            Pident = NA_real_, AlignLength = NA_integer_,
                            Evalue = NA_real_, Bitscore = NA_real_,
                            stringsAsFactors = FALSE)
      append.summary(temp.row)
      next
    }

    # Save assembled barcode contigs
    names(barcode.contigs) = paste0("seq", seq_along(barcode.contigs))
    out.fa = paste0(output.directory, "/sample-barcodes/", samples[i], ".fa")
    write.loci = as.list(as.character(barcode.contigs))
    PhyloProcessR::writeFasta(sequences = write.loci, names = names(write.loci),
                              out.fa, nbchar = 1000000, as.string = TRUE)

    contig.len   = sum(Biostrings::width(barcode.contigs))
    contig.count = length(barcode.contigs)

    #------------------------------------------------------
    # BLAST assembled contigs for species identification
    #------------------------------------------------------
    blast.out = paste0(output.directory, "/sample-barcodes/", samples[i], "_blast.txt")

    if (use.remote.blast) {
      cat(" BLASTing", samples[i], "against NCBI nt (remote)...\n")
      system(paste0(blast.path, "blastn -task dc-megablast -remote -db nt",
                    " -query ", out.fa,
                    " -out ", blast.out,
                    " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen gaps\"",
                    " -max_target_seqs ", hits.per.sample),
             ignore.stdout = quiet, ignore.stderr = quiet)
    } else {
      system(paste0(blast.path, "blastn -task dc-megablast",
                    " -db ", id.blast.db,
                    " -query ", out.fa,
                    " -out ", blast.out,
                    " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen gaps\"",
                    " -max_target_seqs ", hits.per.sample,
                    " -num_threads ", threads),
             ignore.stdout = quiet, ignore.stderr = quiet)
    }

    if (!file.exists(blast.out) || file.info(blast.out)$size == 0) {
      cat(" ", samples[i], ": assembly succeeded but no BLAST hits found.\n")
      temp.row = data.frame(Sample = samples[i],
                            ContigLength = contig.len, ContigCount = contig.count,
                            BestMatch = "No-hit",
                            Pident = NA_real_, AlignLength = NA_integer_,
                            Evalue = NA_real_, Bitscore = NA_real_,
                            stringsAsFactors = FALSE)
    } else {
      match.data = read.table(blast.out, sep = "\t", header = FALSE,
                              stringsAsFactors = FALSE)
      colnames(match.data) = blast.headers
      best = match.data[which.max(match.data$bitscore), ]
      temp.row = data.frame(Sample = samples[i],
                            ContigLength = contig.len,
                            ContigCount  = contig.count,
                            BestMatch    = best$tName,
                            Pident       = round(best$pident, 2),
                            AlignLength  = best$matches,
                            Evalue       = best$evalue,
                            Bitscore     = best$bitscore,
                            stringsAsFactors = FALSE)
    }

    #------------------------------------------------------
    # Append to rolling summary CSV
    #------------------------------------------------------
    append.summary(temp.row)
    print(paste0(samples[i], " barcode scan complete!"))

  }  # end for i

}  # end function

# END SCRIPT
