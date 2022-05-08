# `Extract calibrations` functions take chronograms and extract all node ages as
# taxon pair ages.
# These functions should be named `extract ages`.


#' Use congruification to extract secondary calibrations from a `phylo` or `multiPhylo`
#' object with branch lengths proportional to time.
#'
#' @description This function extracts node ages for each taxon
#'   pair given in `input$tip.labels`. It applies the congruification method
#'   described in Eastman et al. (2013) \doi{10.1111/2041-210X.12051},
#'   implemented with the function [geiger::congruify.phylo()], to create a
#'   `data.frame` of taxon pair node ages that can be used as secondary calibrations.
#' @param input A `phylo` or `multiPhylo` object with branch lengths
#' proportional to time.
#' @param each Boolean, default to `FALSE`: all calibrations are returned in
#' the same `data.frame`. If `TRUE`, calibrations from each chronogram are returned
#' in separate data frames.
#' @return An object of class `calibrations`, which is a `data.frame` (if
#'   `each = FALSE`) or a list of `data.frames` (if `each = TRUE`) of node
#'   ages for each pair of taxon names. You can access the `input` data from which
#'   the calibrations were extracted with attributes(output)$chronograms.
#' @references
#' Eastman et al. (2013) "Congruification: support for time scaling large
#' phylogenetic trees". Methods in Ecology and Evolution, 4(7), 688-691,
#' \doi{10.1111/2041-210X.12051}.
#' @export
extract_calibrations_phylo <- function(input = NULL,
                                       each = FALSE) {
  ##############################################################################
  ##############################################################################
  # Initial argument check
  ##############################################################################
  ##############################################################################
  chronograms <- NULL
  if (inherits(input, "multiPhylo")) {
    chronograms <- input
    xx <- sapply(chronograms, "[", "edge.length")
    xx <- unname(sapply(xx, is.null))
    if (all(xx)) {
      warning("trees in 'multiPhylo' input have no branch lengths. \n There are no calibrations to return!")
      return(NA)
    }
    if (any(xx)) {
      ii <- which(xx)
      warning("Some trees in 'multiPhylo' input have no branch lengths.")
      message(
        "Dropping tree(s) ",
        paste(ii, collapse = " - "),
        " out of the analysis."
      )
      chronograms <- chronograms[which(!xx)]
    }
  }
  if (inherits(input, "phylo")) {
    if (is.null(input$edge.length)) {
      warning("'input' tree has no branch lengths. \n There are no calibrations to return!")
      return(NA)
    }
    chronograms <- list(input)
  }
  if (is.null(chronograms)) {
    stop("'input' must be a 'phylo' or 'multiPhylo' object with branch lengths proportional to time.")
  }
  ##############################################################################
  ##############################################################################
  # Self congruifiction to get a data.frame of taxon pair ages
  ##############################################################################
  ##############################################################################
  if (each) {
    calibrations <- vector(mode = "list")
  } else {
    # we cannot set an empty data frame because nrow depends on the number of nodes available on each tree
    calibrations <- data.frame()
  }
  for (i in seq(length(chronograms))) {
    chronograms[[i]]$tip.label <- gsub(" ", "_", chronograms[[i]]$tip.label) # the battle will never end!
    class(chronograms[[i]]) <- "phylo"
    local_df <- suppressWarnings(
      geiger::congruify.phylo(
        reference = chronograms[[i]],
        target = chronograms[[i]],
        scale = NA,
        ncores = 1
      )
    )$calibrations
    # suppressedWarnings bc of warning message when running
    # geiger::congruify.phylo(reference = chronograms[[i]], target = chronograms[[i]], scale = NA)
    # 		Warning message:
    # In if (class(stock) == "phylo") { :
    # the condition has length > 1 and only the first element will be used
    local_df$reference <- names(chronograms)[i]
    if (each) {
      calibrations <- c(calibrations, list(local_df))
    } else {
      if (i == 1) {
        calibrations <- local_df
      } else {
        calibrations <- rbind(calibrations, local_df)
      }
    }
  }
  ##############################################################################
  ##############################################################################
  # Output
  ##############################################################################
  ##############################################################################
  if (each) {
    names(calibrations) <- names(chronograms)
  }
  attr(calibrations, "chronograms") <- chronograms
  # TODO check that class data frame is also preserved. Might wanna do:
  class(calibrations) <- c(class(calibrations), "calibrations")
  # instead of using structure()
  return(calibrations)
}

#' Use congruification to extract secondary calibrations from a `datelifeResult` object.
#'
#' @inherit extract_calibrations_phylo description
#' @details The function takes a `datelifeResult` object and calls
#' [summarize_datelife_result()] with `summary_format = "phylo_all". This goes
#' from a `datelifeResult` object to a `phylo` or `multiPhylo` object that is
#' passed to [extract_calibrations_phylo()].
#'
#' @param input A `datelifeResult` object.
#' @inheritParams get_all_calibrations
#' @inherit extract_calibrations_phylo return
#' @export
extract_calibrations_dateliferesult <- function(input = NULL,
                                                each = FALSE) {
  phyloall <- suppressMessages(
    summarize_datelife_result(
      datelife_result = input,
      summary_format = "phylo_all"
    )
  )
  res <- extract_calibrations_phylo(
    input = phyloall,
    each = each
  )
  attr(res, "datelife_result") <- input
  class(res) <- c("data.frame", "calibrations")
  return(res)
}
