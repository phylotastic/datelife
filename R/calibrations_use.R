#' Date a given tree topology using a given set of congruified calibrations or ages
#'
#' @description `use_all_calibrations` generates one or multiple chronograms
#'   (i.e., phylogenetic trees with branch lengths proportional to time) by dating
#'   a tree topology given in `phy`, and secondary calibrations given in
#'   `calibrations`, using the algorithm specified in the argument `dating_method`.
#' @inherit use_calibrations details
#' @param phy A `phylo` object to use as tree topology.
#' @param calibrations A `calibrations` object, an output of [get_all_calibrations()].
#' @inheritParams extract_calibrations_phylo
## extract_calibrations_phylo has param each
#' @param dating_method Tree dating algorithm to use. Options are "bladj" or "pathd8"
#'   (Webb et al., 2008, \doi{10.1093/bioinformatics/btn358}; Britton et al., 2007,
#'   \doi{10.1080/10635150701613783}).
#' @inheritDotParams use_calibrations
#' @return A `phylo` or `multiPhylo` object with branch lengths proportional to time.
#' @inheritSection use_calibrations More
#' @inherit use_calibrations details
#' @references
#' Webb, C. O., Ackerly, D. D., & Kembel, S. W. (2008). "Phylocom: software for
#' the analysis of phylogenetic community structure and trait evolution".
#' Bioinformatics, 24(18), \doi{10.1093/bioinformatics/btn358}.
#'
#' Britton, T., Anderson, C. L., Jacquet, D., Lundqvist, S., & Bremer, K. (2007).
#' "Estimating divergence times in large phylogenetic trees". Systematic biology,
#' 56(5), 741-752. \doi{10.1080/10635150701613783}.
#' @export
use_all_calibrations <- function(phy = NULL,
                                 calibrations = NULL,
                                 each = FALSE,
                                 dating_method = "bladj",
                                 ...) {
  #
  # date the tree with bladj, pathd8 if branch lengths:
  # TODO
  # find a way to inherit params for use_calibrations_pathd8 and use_calibrations_bladj
  # phytools method for argument assignation should work, for example:
  # if (methods::hasArg("expand"))
  # if (methods::hasArg("giveup"))
  # if (methods::hasArg("type"))
  # if (methods::hasArg("root_age"))
  # if (methods::hasArg("match_calibrations"))
  if (each) { # if calibrations is a list of data frames:
    res <- use_calibrations_each(
      phy = phy,
      calibrations = calibrations,
      dating_method = dating_method,
      ... # goes to use_calibrations
    )
  } else {
    res <- use_calibrations(phy,
      calibrations,
      dating_method = dating_method,
      ... # goes to use_calibrations, too
    )
  }
  # TODO: make a 'datelife' class for chronograms with calibrations data.frame attached, produced by datelife
  return(res)
}
#' Date a given tree topology by using a given list of calibrations independently,
#' to generate multiple hypothesis of time of divergence
#'
#' @description \code{use_calibrations_each} wraps \code{use_calibrations} to take each set of
#' given calibrations and use it independently as constraints for BLADJ or PATHd8
#' to date a given tree topology.
#' @inheritParams use_all_calibrations
#' @inheritDotParams use_calibrations
# #' @param ... Arguments to pass to \code{use_calibrations}.
#' @return A `multiPhylo` object of trees with branch lengths proportional to time.
#' @inheritSection use_calibrations More
#' @inherit use_calibrations details
#' @export
# You can get the datelifeResult object or the list of chronograms first
# phylo_all1 <- get_datelife_result(input = my_phy)
# phylo_all2 <- summarize_datelife_result(phylo_all1)
# Either will work the same:
# use_calibrations_each(phy = my_phy, phylo_all = phylo_all1)
# use_calibrations_each(phy = my_phy, phylo_all = phylo_all2)
# You can also get the list of calibrations outside
# my_calibrations <- get_all_calibrations(input = phylo_all1, each = TRUE)
# use_calibrations_each(phy = my_phy, calibrations = my_calibrations)
# All this means that you can use your own set of calibrations, not necessarily the ones found only in datelife.
# phy <- tax_otolall[[1]]
# calibrations <- tax_eachcalall[[1]]
use_calibrations_each <- function(phy = NULL,
                                  calibrations = NULL,
                                  ...) {

  # date the tree with bladj, or pathd8 if branch lengths:
  if (inherits(calibrations, "data.frame")) {
    calibrations <- structure(list(calibrations), class = c("list", "calibrations"))
  }
  chronograms <- lapply(calibrations, function(x) {
    use_calibrations(
      phy = phy,
      calibrations = x,
      ...
    )
  })
  class(chronograms) <- c("multiPhylo", "datelife")
  names(chronograms) <- names(calibrations)
  return(chronograms = chronograms)
}
#' Date a given tree topology using a combined set of given calibrations
#'
#' @description `use_calibrations` combines all given calibrations and uses them as
#' constraints to perform a dating analysis on a given tree topology, using BLADJ
#' if it has no branch lengths, or PATHd8 if the given tree topology has initial
#' branch lengths.
#'
#' @inheritParams use_all_calibrations
# use_all_calibrations has params phy, calibrations, dating_method
#' @inheritParams use_calibrations_bladj
# use_calibrations_bladj has param type
#' @inheritDotParams use_calibrations_pathd8
#' @return A `phylo` object with branch lengths proportional to time.
#' @section More: The output object stores the used `calibrations` and `dating_method` as
#'   `attributes(output)$datelife_calibrations` and `attributes(output)$dating_method`.
#' @details
#' If `phy` has no branch lengths, `dating_method` is ignores, and the function applies secondary
#' calibrations to date the tree with the BLADJ algorithm. See [make_bladj_tree()] and [use_calibrations_bladj()].
#' If `phy` has branch lengths, the function can use the PATHd8 algorithm. See [use_calibrations_pathd8()].
#' @export
use_calibrations <- function(phy = NULL,
                             calibrations = NULL,
                             dating_method = "bladj",
                             type = "median",
                             ...) {
  # check that phy names are in calibrations.df: done in match_all_calibrations inside use_calibrations_bladj
  message("... Using calibrations to date a tree topology.")
  exit <- FALSE
  if (inherits(phy, "multiPhylo")) {
    message(message_multiphylo())
    phy <- phy[[1]]
  }
  if (!inherits(phy, "phylo")) {
    exit <- TRUE
    msg1 <- "'phy' is NOT a phylo object. Check this."
  } else {
    msg1 <- "'phy' is a phylo object. You are good."
  }
  # enhance: add a check for calibrations object structure
  if (!inherits(calibrations, "data.frame")) {
    exit <- TRUE
    msg2 <- "'calibrations' is NOT a data frame. Check this."
  } else {
    msg2 <- "'calibrations' is a data frame. You are good."
  }
  if (exit) {
    stop(
      "This function requires valid values for both arguments 'phy' and 'calibrations'.\n",
      msg1, "\n",
      msg2, "\n"
    )
  }
  dating_method <- match.arg(tolower(dating_method), c("bladj", "pathd8"))
  if (dating_method %in% "bladj") {
    phy$edge.length <- NULL
    chronogram <- use_calibrations_bladj(phy, calibrations, type = type)
  }
  if (dating_method %in% "pathd8") {
    if (is.null(phy$edge.length)) {
      message("\nPATHd8 requires initial branch lengths and 'phy' has no branch lengths, using dating_method = bladj instead.\n")
      chronogram <- use_calibrations_bladj(phy, calibrations, type = type)
      dating_method <- "bladj"
    } else {
      chronogram <- use_calibrations_pathd8(phy, calibrations, ...)
    }
  }
  # add dating_method and calibrations attribute to chronogram
  attr(chronogram, "datelife_calibrations") <- calibrations
  attr(chronogram, "dating_method") <- dating_method
  class(chronogram) <- c(class(chronogram), "datelife")
  message("Success!")
  return(chronogram)
}
