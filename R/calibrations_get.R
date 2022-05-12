# `Get calibrations` functions take taxon names and perform a chronogram search
# on a database. Then, they extract all node ages from matching chronograms as
# taxon pair ages using `extract calibrations` functions.
# These functions should be named `get all ages`.

#' Get secondary calibrations from a chronogram database for a set of given taxon names
#'
#' @aliases datelife_calibrations
#' @description \code{get_all_calibrations} performs a [datelife_search()]
#' and gets divergence times (i.e., secondary calibrations) from a chronogram
#' database for each taxon name pair given as \code{input}.
#'
#' @inheritParams datelife_search
#' @inheritParams extract_calibrations_phylo
#' @inherit extract_calibrations_phylo return
#' @export
get_all_calibrations <- function(input = NULL,
                                 each = FALSE) {
  #
  if (inherits(input, "datelifeQuery")) {
    res <- get_calibrations_datelifequery(
      datelife_query = input,
      each = each
    )
    return(res)
  }
  if (inherits(input, "phylo") | inherits(input, "multiPhylo")) {
    if (inherits(input, "multiPhylo")) {
      input <- unname(unlist(lapply(input, "[", "tip.label")))
    } else {
      input <- input$tip.label
    }
  }
  if (inherits(input, "character")) {
    input <- unique(gsub(" ", "_", input))
    res <- get_calibrations_vector(
      input = input,
      each = each
    )
  }
  return(res)
}

#' Search and extract secondary calibrations for a given character
#' vector of taxon names
#'
#' @description The function searches DateLife's local
#' database of phylogenetic trees with branch lengths proportional to time
#' (chronograms) with [datelife_search()], and extracts available node ages
#' for each pair of given taxon names with [extract_calibrations_phylo()].
#'
#' @details The function calls [datelife_search()]
#' with `summary_format = "phylo_all"` to get all chronograms in the database
#' containing at least two taxa in `input`, and generates a `phylo`
#' or `multiPhylo` object object that will be passed to
#' [extract_calibrations_phylo()].
#'
#' @param input A character vector of taxon names.
#' @inheritParams get_all_calibrations
#' @inherit extract_calibrations_phylo return
#' @export
get_calibrations_vector <- function(input = NULL,
                                    each = FALSE) {
  # TODO: is_datelife_search_input function, or any type of input format checking
  # i.e., a function to trap the case were input is a list
  phyloall <- datelife_search(
    input = input,
    summary_format = "phylo_all"
  )

  res <- extract_calibrations_phylo(
    input = phyloall,
    each = each
  )
  attr(res, "datelife_result") <- attributes(phyloall)$datelife_result
  class(res) <- c("data.frame", "calibrations")
  return(res)
}
#' Search and extract available secondary calibrations for taxon names in a given
#' `datelifeQuery` object
#'
#' @param datelife_query A `datelifeQuery` object.
#' @inheritParams get_all_calibrations
#' @inherit get_calibrations_vector description details
#' @inherit extract_calibrations_phylo return
#' @export
get_calibrations_datelifequery <- function(datelife_query = NULL,
                                           each = FALSE) {
  if (suppressMessages(!is_datelife_query(datelife_query))) {
    stop("'datelife_query' is not a 'datelifeQuery' object.")
  }
  phyloall <- datelife_search(
    input = datelife_query,
    summary_format = "phylo_all"
  )
  res <- extract_calibrations_phylo(
    input = phyloall,
    each = each
  )
  attr(res, "datelife_result") <- attributes(phyloall)$datelife_result
  class(res) <- c("data.frame", "calibrations")
  return(extract_calibrations_phylo(
    input = phyloall,
    each = each
  ))
}
