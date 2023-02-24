#' Get genetic data from the Barcode of Life Database (BOLD) for a set of taxon names.
#'
#' @description `get_bold_data` uses taxon names from a tree topology, a character
#' vector of names or a `datelifeQuery` object, to search for genetic markers
#' in the Barcode of Life Database (BOLD).
#'
#' @inheritParams datelife_search
#' @param marker A character vector indicating the gene from BOLD system to be
#' used for branch length estimation. It searches "COI" marker by default.
#' @inheritDotParams get_otol_synthetic_tree
#' @return A `phylo` object. If there are enough BOLD sequences available for the
#'   `input` taxon names, the function returns a tree with branch lengths proportional
#'   to relative substitution rate. If not enough BOLD sequences are available
#'   for the `input` taxon names, the function returns the topology given as
#'   `input`, or a synthetic Open Tree of Life for the taxon names given in
#'   `input`, obtained with [get_otol_synthetic_tree()].
#' @details
#'   If `input` is a `phylo` object or a newick string, it is used as backbone topology.
#'   If `input` is a character vector of taxon names, an induced synthetic OpenTree
#'   subtree is used as backbone.
#' @importFrom BiocManager install
#' @export
get_bold_data <- function(input = c("Rhea americana", "Struthio camelus", "Gallus gallus"),
                          marker = "COI",
                          ...) {
  ##############################################################################
  # checking arguments and packages
  ##############################################################################
  if (!requireNamespace("msa", quietly=TRUE)) {
    stop("'msa' package is not installed. Please install it from Bioconductor with",
         " BiocManager::install('msa', dependencies = TRUE)")
  }
  if (!requireNamespace("Biostrings", quietly=TRUE)) {
    stop("'Biostring' package is not installed. Please install it from Bioconductor with",
         " BiocManager::install('Biostrings', dependencies = TRUE)")
  }
  # input check (accepts newick strings too)
  datelife_query <- input
  if (suppressMessages(!is_datelife_query(input))) {
    datelife_query <- make_datelife_query(input)
  }
  ##############################################################################
  # searching sequences in BOLD
  ##############################################################################
  message("---> Searching for ", marker,
          " sequences available in the Barcode of Life Database (BOLD) for 'input' taxon names.")
  phy$edge.length <- NULL # making sure there are no branch lengths in phy
  phy$tip.label <- gsub(" ", "_", phy$tip.label) # so phangorn::acctran works
  bold_input <- gsub("_", " ", phy$tip.label) # so bold search works
  sequences <- c()
  progression <- utils::txtProgressBar(min = 0, max = length(bold_input), style = 3)
  for (i in seq(length(bold_input))) {
    ss <- bold::bold_seqspec(taxon = bold_input[i])
    if (inherits(ss, "data.frame")) {
      sequences <- rbind(sequences, ss)
    }
    # allows up to 335 names, then it gives Error: Request-URI Too Long (HTTP 414)
    # even if marker is specified, it will return other markers,
    # so in here we just get all sequences and then filter after
    utils::setTxtProgressBar(progression, i)
  }
  # cat("\n") # just to make the progress bar look better
  sequences <- sequences[grepl(marker, sequences$markercode), ] # filter other markers
  if (length(sequences) == 1) {
    # it is length == 80 when there is at least 1 sequence available;
    # if this is TRUE, it means there are no sequences in BOLD for the set of input taxa.
    # if (!use_tnrs) message("Setting 'use_tnrs = TRUE' might change this, but it can be slow.\n")
    message("* Names in 'input' do not match the Barcode of Life Database (BOLD) specimen records")
    message("* No sequences were found in the Barcode of Life Database (BOLD) for the given 'input' taxon names")
    return(NA)
  }
  message("BOLD sequence search done!")
  return(sequences)
}
