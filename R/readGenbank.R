#' @title readGenBank
#'
#' @description Function for running the program spades to assemble short read sequencing data
#'
#' @param genbank.file path to a folder of sequence alignments in phylip format.
#'
#' @param overwrite path to a folder of sequence alignments in phylip format.
#'
#' @param rep.origin path to a folder of sequence alignments in phylip format.
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

#Add something to start reference at point of replication
#
# readGenBank = function (genbank.file = NULL,
#                         partial = NA,
#                         ret.seq = TRUE,
#                         verbose = FALSE) {
#
#   file = "/Volumes/Rodents/Mitogenomes/Crocidura.gb"
#
#   #Reads in text file
#   text = readLines(file)
#
#   if (missing(text)) {
#     if (missing(file))
#       stop("One of text or file must be specified.")
#   }
#
#   #PARSE GENBANK FUNCTION HERE
#   fldlines = grepl("^[[:upper:]]+[[:space:]]+", text)
#   fldfac = cumsum(fldlines)
#   fldnames = gsub("^([[:upper:]]+).*", "\\1", text[fldlines])[fldfac]
#   spl = split(text, fldnames)
#
#   #Read locus function here
#   result.data = list(LOCUS = strsplit(spl[["LOCUS"]], "[\\t]+", spl[["LOCUS"]])[[1]])
#   feat.data = spl[["FEATURES"]]
#
#
#
#
#
#
#
#
#
#   if (substr(feat.data[1], 1, 8) == "FEATURES")
#     feat.data = feat.data[-1]
#   fttypelins = grepl("^( {5}|\\t)[[:alnum:]'_]+[[:space:]]+(complement|join|order|[[:digit:]<,])", feat.data)
#   featfactor = cumsum(fttypelins)
#
#   chr = "unk"
#   totsources = length(grep("[[:space:]]+source[[:space:]]+[<[:digit:]]",
#                            feat.data[which(fttypelins)]))
#   numsources = 0
#   everhadchr = FALSE
#   type = gsub("[[:space:]]+([[:alnum:]_']+).*", "\\1",
#               feat.data[1])
#   attrstrts = cumsum(grepl("^[[:space:]]+/[^[:space:]]+($|=([[:digit:]]|\"))",
#                            feat.data))
#
#   feat.data = tapply(feat.data, attrstrts, function(x) {
#     paste(gsub("^[[:space:]]+", "", x), collapse = "")
#   }, simplify = TRUE)
#
#   rawlocstring = feat.data[1]
#   rngstr = gsub("^[[:space:]]*[[:alnum:]_']+[[:space:]]+((complement\\(|join\\(|order\\(|[[:digit:]<]+).*)", "\\1", rawlocstring)
#
#   feat.data = feat.data[-1]
#   if (length(feat.data)) {
#
#     num =  grepl('=[[:digit:]]+(\\.[[:digit:]]+){0,1}$', feat.data)
#     val = gsub('[[:space:]]*/[^=]+($|="{0,1}([^"]*)"{0,1})', "\\2", feat.data)
#     mapply(function(val, num) {
#       if(nchar(val)==0)
#         TRUE
#       else if(num)
#         as.numeric(val)
#       else
#         val
#     }, val = val, num = num, SIMPLIFY=FALSE)
#
#     attrs = val
#     names(attrs) = gsub("^[[:space:]]*/([^=]+)($|=[^[:space:]].*$)",
#                         "\\1", feat.data)
#   }
#
#
#
#
#
#
#
#
#           else if ("strain" %in% names(attrs)) {
#             chr <<- if (totsources == 1)
#               attrs$strain
#             else paste(attrs$strain, numsources, sep = ":")
#           }
#           else {
#             chr <<- if (totsources == 1)
#               attrs$organism
#             else paste(attrs$organism, numsources, sep = ":")
#           }
#         }
#       }
#       else {
#         attrs = list()
#       }
#       make_feat_gr(str = rngstr, chr = chr, ats = c(type = type,
#                                                     attrs), partial = partial)
#     }
#     if (verbose)
#       message(Sys.time(), " Starting feature parsing")
#     resgrs = tapply(lines, featfactor, do_readfeat, simplify = FALSE,
#                     partial = partial)
#     if (verbose)
#       message(Sys.time(), " Done feature parsing")
#     resgrs
#   }
#
#
#
#   seqtype = .seqTypeFromLocus(resthang$LOCUS)
#   resthang$ORIGIN = if (ret.seq)
#     readOrigin(spl[["ORIGIN"]], seqtype = seqtype)
#
#
#
#
#
#   ret = make_gbrecord(rawgbk = prsed, verbose = verbose)
#   if (!ret.seq)
#     sequence(ret) = NULL
#   ret
# }
#
#
# function (file, text = readLines(file), partial = NA, verbose = FALSE,
#           ret.anno = TRUE, ret.seq = TRUE)
# {
#   if (!ret.anno && !ret.seq)
#     stop("Must return at least one of annotations or sequence.")
#   bf = proc.time()["elapsed"]
#   if (missing(text) && !file.exists(file))
#     stop("No text provided and file does not exist or was not specified. Either an existing file or text to parse must be provided.")
#   if (length(text) == 1)
#     text = fastwriteread(text)
#   fldlines = grepl(prime_field_re, text)
#   fldfac = cumsum(fldlines)
#   fldnames = gsub("^([[:upper:]]+).*", "\\1", text[fldlines])[fldfac]
#   spl = split(text, fldnames)
#   resthang = list(LOCUS = readLocus(spl[["LOCUS"]]))
#   resthang[["FEATURES"]] = readFeatures(spl[["FEATURES"]],
#                                         source.only = !ret.anno, partial = partial)
#   seqtype = .seqTypeFromLocus(resthang$LOCUS)
#   resthang$ORIGIN = if (ret.seq)
#     readOrigin(spl[["ORIGIN"]], seqtype = seqtype)
#   else NULL
#   if (ret.anno) {
#     resthang2 = mapply(function(field, lines, verbose) {
#       switch(field, DEFINITION = readDefinition(lines),
#              ACCESSION = readAccession(lines), VERSION = readVersions(lines),
#              KEYWORDS = readKeywords(lines), SOURCE = readSource(lines),
#              NULL)
#     }, lines = spl, field = names(spl), SIMPLIFY = FALSE,
#     verbose = verbose)
#     resthang2$FEATURES = resthang2$FEATURES[sapply(resthang2$FEATURES,
#                                                    function(x) length(x) > 0)]
#     resthang2 = resthang2[!names(resthang2) %in% names(resthang)]
#     resthang = c(resthang, resthang2)
#   }
#   origin = resthang$ORIGIN
#   if (ret.seq && length(origin) > 0) {
#     typs = sapply(resthang$FEATURES, function(x) x$type[1])
#     srcs = fill_stack_df(resthang$FEATURES[typs == "source"])
#     dss = extractAt(origin, ranges(srcs))
#     names(dss) = as.character(seqnames(srcs))
#     if (!ret.anno)
#       resthang = dss
#     else resthang$ORIGIN = dss
#   }
#   else if (!ret.anno) {
#     stop("Asked for only sequence (ret.anno=FALSE) from a file with no sequence information")
#   }
#   af = proc.time()["elapsed"]
#   if (verbose)
#     message("Done Parsing raw GenBank file text. [ ", af -
#               bf, " seconds ]")
#   resthang
# }
#
#
#
# parseGenBank = function(file, text = readLines(file),  partial = NA,
#                         verbose = FALSE,
#                         ret.anno = TRUE,
#                         ret.seq = TRUE) {
#   if(!ret.anno && !ret.seq)
#     stop("Must return at least one of annotations or sequence.")
#   bf = proc.time()["elapsed"]
#   if(missing(text) && !file.exists(file))
#     stop("No text provided and file does not exist or was not specified. Either an existing file or text to parse must be provided.")
#   if(length(text) == 1)
#     text = fastwriteread(text)
#
#   fldlines = grepl(prime_field_re, text)
#   fldfac = cumsum(fldlines)
#   fldnames = gsub("^([[:upper:]]+).*", "\\1", text[fldlines])[fldfac]
#
#   spl = split(text, fldnames)
#
#   resthang = list(LOCUS = readLocus(spl[["LOCUS"]]))
#   resthang[["FEATURES"]] = readFeatures(spl[["FEATURES"]],
#                                         source.only=!ret.anno,
#                                         partial = partial)
#   seqtype = .seqTypeFromLocus(resthang$LOCUS)
#   resthang$ORIGIN = if(ret.seq)
#     readOrigin(spl[["ORIGIN"]],
#                seqtype = seqtype)
#   else NULL
#
#   if(ret.anno) {
#     resthang2 = mapply(function(field, lines, verbose) {
#       switch(field,
#              DEFINITION = readDefinition(lines),
#              ACCESSION = readAccession(lines),
#              VERSION = readVersions(lines),
#              KEYWORDS = readKeywords(lines),
#              SOURCE = readSource(lines),
#              ## don't read FEATURES, ORIGIN, or LOCUS because they are
#              ## already in resthang from above
#              NULL)
#     }, lines = spl, field = names(spl), SIMPLIFY=FALSE, verbose = verbose)
#     resthang2$FEATURES = resthang2$FEATURES[sapply(resthang2$FEATURES,
#                                                    function(x) length(x)>0)]
#     resthang2 = resthang2[!names(resthang2) %in% names(resthang)]
#     resthang = c(resthang, resthang2)
#   }
#   ##DNAString to DNAStringSet
#   origin = resthang$ORIGIN
#   if(ret.seq && length(origin) > 0) {
#     typs = sapply(resthang$FEATURES, function(x) x$type[1])
#     srcs = fill_stack_df(resthang$FEATURES[typs == "source"])
#     ## dss = DNAStringSet(lapply(GRanges(ranges(srcs), function(x) origin[x])))
#     dss = extractAt(origin, ranges(srcs))
#     names(dss) = as.character(seqnames(srcs))
#     if(!ret.anno)
#       resthang = dss
#     else
#       resthang$ORIGIN = dss
#   } else if (!ret.anno) { ##implies ret.seq is TRUE
#     stop("Asked for only sequence (ret.anno=FALSE) from a file with no sequence information")
#   }
#   af = proc.time()["elapsed"]
#   if(verbose)
#     message("Done Parsing raw GenBank file text. [ ", af-bf, " seconds ]")
#   resthang
#
#
# }

