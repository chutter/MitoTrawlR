#' @title setupCheck
#'
#' @description Checks whether all external programs required by MitoTrawlR
#'   can be found on disk. If \code{anaconda.environment} is provided, all
#'   tool paths are automatically set to \code{<anaconda.environment>/bin}.
#'   Individual path overrides can be supplied if tools are installed in
#'   non-standard locations. Prints a found/not-found message for each tool and
#'   returns a single logical indicating overall pass/fail.
#'
#' @param anaconda.environment path to the root of a conda environment (e.g.,
#'   \code{"/opt/miniconda3/envs/MitoTrawlR"}). When provided all tool paths
#'   are derived from \code{<anaconda.environment>/bin}.
#'
#' @param fastp.path path to the directory containing \code{fastp}. Overridden
#'   by \code{anaconda.environment} when that is non-NULL.
#'
#' @param samtools.path path to the directory containing \code{samtools}.
#'
#' @param bwa.path path to the directory containing \code{bwa}.
#'
#' @param spades.path path to the directory containing \code{spades.py}.
#'
#' @param bbmap.path path to the directory containing \code{bbmap.sh}.
#'
#' @param blast.path path to the directory containing \code{blastn} and
#'   \code{makeblastdb}.
#'
#' @param mafft.path path to the directory containing \code{mafft}.
#'
#' @param iqtree.path path to the directory containing \code{iqtree3}.
#'
#' @param trimAl.path path to the directory containing \code{trimal}.
#'
#' @param julia.path path to the directory containing \code{julia} (currently
#'   unused; reserved for TAPER support).
#'
#' @param taper.path path to the directory containing the TAPER script
#'   (currently unused).
#'
#' @param tRNAscan.path path to the directory containing \code{tRNAscan-SE}.
#'
#' @return logical; TRUE if all required tools were found, FALSE otherwise.
#'
#' @examples
#' \dontrun{
#' MitoTrawlR::setupCheck(
#'   anaconda.environment = "/opt/miniconda3/envs/MitoTrawlR"
#' )
#' }
#'
#' @export

setupCheck = function(anaconda.environment = NULL,
                      fastp.path = NULL,
                      samtools.path = NULL,
                      bwa.path = NULL,
                      spades.path = NULL,
                      bbmap.path = NULL,
                      blast.path = NULL,
                      cap3.path = NULL,
                      mafft.path = NULL,
                      iqtree.path = NULL,
                      trimAl.path = NULL,
                      julia.path = NULL,
                      taper.path = NULL,
                      tRNAscan.path = NULL) {

  #anaconda.environment = "/Users/chutter/conda/MitoTrawlR"

  if (is.null(anaconda.environment) == FALSE) {
    fastp.path    = paste0(anaconda.environment, "/bin")
    samtools.path = paste0(anaconda.environment, "/bin")
    bwa.path      = paste0(anaconda.environment, "/bin")
    spades.path   = paste0(anaconda.environment, "/bin")
    bbmap.path    = paste0(anaconda.environment, "/bin")
    blast.path    = paste0(anaconda.environment, "/bin")
    cap3.path     = paste0(anaconda.environment, "/bin")
    mafft.path    = paste0(anaconda.environment, "/bin")
    iqtree.path   = paste0(anaconda.environment, "/bin")
    trimAl.path   = paste0(anaconda.environment, "/bin")
    taper.path    = paste0(anaconda.environment, "/bin")
    julia.path    = paste0(anaconda.environment, "/bin")
    tRNAscan.path = paste0(anaconda.environment, "/bin")
  }

  check.tool = function(path, tool) {
    found = file.exists(paste0(path, "/", tool))
    if (found) {
      message(tool, " was found.")
    } else {
      message(tool, " could not be found.")
    }
    found
  }

  results = c(
    check.tool(fastp.path,     "fastp"),
    check.tool(samtools.path,  "samtools"),
    check.tool(bwa.path,       "bwa"),
    check.tool(spades.path,    "spades.py"),
    check.tool(bbmap.path,     "bbmap.sh"),
    check.tool(blast.path,     "blastn"),
    check.tool(blast.path,     "makeblastdb"),
    check.tool(cap3.path,      "cap3"),
    check.tool(mafft.path,     "mafft"),
    check.tool(iqtree.path,    "iqtree3"),
    check.tool(trimAl.path,    "trimal"),
    check.tool(tRNAscan.path,  "tRNAscan-SE")
  )

  pass = all(results)
  if (pass) {
    message("All tools found. MitoTrawlR environment is ready.")
  } else {
    message("Some tools were not found. Please check the paths above.")
  }

  return(pass)

}#end function
