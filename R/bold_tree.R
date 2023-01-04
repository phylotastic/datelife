#' Use genetic data from the Barcode of Life Database (BOLD) to reconstruct branch lengths on a tree.
#'
#' @description `make_bold_otol_tree` takes taxon names from a tree topology or
#' a vector of names to search for genetic markers in the Barcode of Life Database
#' (BOLD), create an alignment, and reconstruct branch lengths on a tree topology
#' with Maximum Likelihood.
#'
#' @inheritParams datelife_search
#' @param marker A character vector indicating the gene from BOLD system to be
#' used for branch length estimation.
#' @param chronogram Default to `TRUE`, branch lengths returned are estimated with
#' [ape::chronoMPL()]. If `FALSE`, branch lengths returned are estimated with
#' [phangorn::acctran()] and represent relative substitution rates.
#' @param doML Default to `FALSE`. If `TRUE`, it does a ML branch length optimization
#' with [phangorn::optim.pml()].
# Only relevant if chronogram = `TRUE`.
#' @param aligner A character vector indicating whether to use MAFFT or MUSCLE
#' to align BOLD sequences. It is not case sensitive. Default to MUSCLE,
#' supported using the [msa](https://bioconductor.org/packages/release/bioc/html/msa.html)
#' package from Bioconductor, which needs to be installed using [BiocManager::install()].
#' @inheritParams get_otol_synthetic_tree
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
make_bold_otol_tree <- function(input = c("Rhea americana", "Struthio camelus", "Gallus gallus"),
                                marker = "COI",
                                otol_version = "v3",
                                chronogram = TRUE,
                                doML = FALSE,
                                aligner = "muscle",
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
  phy <- datelife_query$phy
  if (!inherits(phy, "phylo")) {
    message("A tree topology is needed as 'input' to reconstruct branch lengths with BOLD data.")
    phy <- get_otol_synthetic_tree(input = datelife_query, otol_version = otol_version, ...)
    # otol returns error with missing taxa in v3 of rotl
  }
  attr(phy, "datelife_query") <- datelife_query
  ##############################################################################
  # searching sequences in BOLD
  # is function get_bold_data now
  ##############################################################################
  message("... Searching for ", marker, " sequences from 'input' taxa available in BOLD.")
  aligner <- match.arg(arg = tolower(aligner), choices = c("muscle", "mafft"), several.ok = FALSE)
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
  cat("\n") # just to make the progress bar look better
  sequences <- sequences[grepl(marker, sequences$markercode), ] # filter other markers
  if (length(sequences) == 1) {
    # it is length == 80 when there is at least 1 sequence available;
    # if this is TRUE, it means there are no sequences in BOLD for the set of input taxa.
    # if (!use_tnrs) message("Setting 'use_tnrs = TRUE' might change this, but it can be slow.\n")
    warning("Names in input do not match BOLD specimen records; no sequences
			were found in BOLD for the set of input taxa.\nReturning a tree topology with no branch lengths from Open Tree of Life!")
    phy$tip.label <- gsub(" ", "_", phy$tip.label)
    attr(phy, "datelife_query") <- datelife_query
    return(phy)
  }
  message("BOLD sequence search done!")
  ##############################################################################
  # checking obtained sequences
  ##############################################################################
  sequences$nucleotide_ATGC <- gsub("[^A,T,G,C]", "", sequences$nucleotides) # preserve good nucleotide data, i.e., only A,T,G,C
  sequences$nucleotide_ATGC_length <- unlist(lapply(sequences$nucleotide_ATGC, nchar)) # add a column in data frame, indicating the amount of good information contained in sequences#nucelotides (ATGC)
  col_number <- max(sapply(strsplit(sequences$nucleotides, ""), length))
  final.sequences <- matrix("-",
    nrow = length(bold_input),
    ncol = col_number
  )
  final.sequences.names <- rep(NA, length(bold_input))
  row.index <- 0
  taxa.to.drop <- c()
  for (i in bold_input) {
    row.index <- row.index + 1
    taxon.index <- which(grepl(i, sequences$species_name))
    # if there are no sequences from any taxon, taxon.index is empty
    # but we make sure this is filtered steps before
    if (length(taxon.index) > 0) {
      seq.index <- which.max(sequences$nucleotide_ATGC_length[taxon.index])
      # sequences[taxon.index,][seq.index,]
      seq <- strsplit(sequences$nucleotides[taxon.index][seq.index], split = "")[[1]]
      final.sequences[row.index, sequence(length(seq))] <- seq
    } else {
      taxa.to.drop <- c(taxa.to.drop, i)
    }
    final.sequences.names[row.index] <- i
  }
  rownames(final.sequences) <- gsub(" ", "_", final.sequences.names)
  # final.sequences <- final.sequences[!is.na(final.sequences.names),]
  # taxa.to.drop <- phy$tip.label[which(!phy$tip.label %in% rownames(final.sequences))]
  # if there are only two sequences in the alignment phangorn::acctran will throw an error
  # this usually happens when the input/otol tree has only two tips
  if (length(bold_input) - length(taxa.to.drop) <= 2) {
    message(
      "BOLD sequences found only for: ",
      bold_input[which(!bold_input %in% taxa.to.drop)]
    )
    warning("There are not enough sequences available in BOLD to reconstruct
						branch lengths. \nReturning a tree topology with no branch lengths from Open Tree of Life!")
    phy$tip.label <- gsub(" ", "_", phy$tip.label)
    attr(phy, "datelife_query") <- datelife_query
    return(phy)
  }
  if (length(taxa.to.drop) > 0) {
    taxa.to.drop.print <- paste(taxa.to.drop, collapse = " | ")
    message(marker, " sequences absent for ", taxa.to.drop.print, ".\nDropping taxa from tree.")
    # warning("No ", marker, " sequences found for ", taxa.to.drop.print, "...", "\n", "\t", "Taxa dropped from tree.")
    taxa.to.drop <- gsub(" ", "_", taxa.to.drop)
    phy <- ape::drop.tip(phy, taxa.to.drop)
  }
  ##############################################################################
  # aligning
  ##############################################################################
  if ("mafft" %in% aligner) {
    message("... Aligning BOLD sequences with MAFFT.")
    alignment <- ape::as.DNAbin(final.sequences)
    alignment <- phangorn::as.phyDat(ips::mafft(alignment))
    message("MAFFT alignment done!")
  }
  if ("muscle" %in% aligner) {
    message("... Aligning BOLD sequences with MUSCLE.")
    in_bold <- sapply(1:nrow(final.sequences), function(x) !"-" %in% final.sequences[x, ])
    vector.sequences <- sapply(1:nrow(final.sequences), function(x) paste0(final.sequences[x, ], collapse = ""))

    msa.sequences <- Biostrings::DNAStringSet(vector.sequences)
    names(msa.sequences) <- rownames(final.sequences)
    msa.sequences <- msa.sequences[in_bold]
    alignment <- msa::msaMuscle(msa.sequences, type = "dna")
    alignment <- phangorn::as.phyDat(alignment)
    message("MUSCLE alignment done!")
    # alignment MUST BE OF CLASS phyDat TO BE READ BY ACCTRAN on next step
  }
  ##############################################################################
  # branch length estimation, using topology in phy
  ##############################################################################
  message("... Estimating tree with PML.")
  xx <- phangorn::acctran(ape::multi2di(phy), alignment)
  pml.object <- phangorn::pml(xx, data = alignment)
  phy <- pml.object$tree
  if (!ape::is.binary(pml.object$tree)) {
    message("Resulting PML tree is non-dichotomous.\nResolving with multi2di.")
    pml.object$tree <- ape::multi2di(pml.object$tree)
    phy <- pml.object$tree
  }
  message("PML tree obtained!")
  ##############################################################################
  # dating
  ##############################################################################
  if (chronogram) {
    message("... Dating PML tree with chronoMPL.")
    pml.object$tree <- ape::chronoMPL(pml.object$tree, se = FALSE, test = FALSE)
    phy <- pml.object$tree
    message("chronoPML chronogram obtained!")
    if (any(phy$edge.length < 0)) {
      warning("There are negative branch lengths in chronoMPL chronogram.\n")
      if (doML) {
        warning("Can't do ML branch length optimization.\n")
      }
    } else {
      if (doML) {
        message("... Optimizing chronoPML chronogram.")
        phy <- phangorn::optim.pml(pml.object,
          data = alignment,
          rearrangement = "none",
          optRooted = TRUE,
          optQ = TRUE
        )$tree
        message("chronoPML chronogram optimized!")
      }
    }
  }
  phy$tip.label <- gsub(" ", "_", phy$tip.label)
  attr(phy, "datelife_query") <- datelife_query
  message("Tree branch lengths reconstructed!")
  return(phy)
}
