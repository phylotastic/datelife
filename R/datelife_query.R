#' Go from taxon names to a `datelifeQuery` object
#'
## #' @description
#'
#' @param input Taxon names as one of the following:
#' \describe{
#' 	 \item{A character vector of taxon names}{With taxon names as a single comma separated starting or concatenated with [c()].}
#' 	 \item{A phylogenetic tree with taxon names as tip labels}{As a `phylo` or `multiPhylo`
#' 	 			object, OR as a newick character string.}
#' }
#' @param use_tnrs Whether to use Open Tree of Life's Taxonomic Name Resolution Service (TNRS)
#'   to process input taxon names. Default to `TRUE`, it corrects misspellings and
#'   taxonomic name variations with [tnrs_match()], a wrapper of [rotl::tnrs_match_names()].
# #' @param use_tnrs Boolean; default to `FALSE`. If `TRUE`, use OpenTree's services
# #'   to resolve names. This can dramatically improve the chance of matches, but also
# #'   take much longer.
# #' @param approximate_match Boolean; default to `TRUE`: use a slower TNRS to
# #'   correct misspellings, increasing the chance of matches (including false matches).
#' @param get_spp_from_taxon Whether to search ages for all species belonging to a
#'   given taxon or not. Default to `FALSE`. If `TRUE`, it must have same length as input.
#'   If input is a newick string with some clades it will be converted to a `phylo`
#'   object, and the order of `get_spp_from_taxon` will match `phy$tip.label`.
#' @param taxonomic_source Used if `get_spp_from_taxon = TRUE`. A character vector with the desired taxonomic sources.
#'  Options are "ott", "ncbi", "gbif" or "irmng". The function defaults to "ott".
#' @return A `datelifeQuery` object, which is a list of three elements:
#' \describe{
#' 	 \item{$phy}{A `phylo` object or `NA`, if input is not a tree.}
#' 	 \item{$cleaned_names}{A character vector of cleaned taxon names.}
#' 	 \item{$ott_ids}{A numeric vector of OTT ids if `use_tnrs = TRUE`, or `NULL` if `use_tnrs = FALSE`.}
#' }
#' @details It processes `phylo` objects and newick character string inputs
#'   with [input_process()]. If `input` is a `multiPhylo` object, only the first `phylo`
#'   element will be used. Similarly, if an `input` newick character string has multiple trees,
#'   only the first one will be used.
#' @export
make_datelife_query <- function(input = c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"),
                                use_tnrs = TRUE,
                                get_spp_from_taxon = FALSE,
                                taxonomic_source = "ott") {
  # enhance: add mapped (has tnrs been performed?) and matched (was it matched successfully?) element to phylo object
  # add one for each taxonomy queried: ott, catalogue of life (also contains fossils), enciclopedia of life (common names)
  if (suppressMessages(is_datelife_query(input))) {
    message("'input' is already a 'datelifeQuery' object.")
    return(input)
  }
  # input_process determines if input is newick and transforms it to phylo
  # if input is phylo or multiphylo it will also check if it is a good tree
  # if input is anything else, it will return NA:
  phy_new <- input_process(input = input)
  use_tnrs_global <- FALSE
  if (use_tnrs | any(get_spp_from_taxon)) {
    use_tnrs_global <- TRUE
  }
  message("... Making a DateLife query:")
  if (inherits(phy_new, "phylo")) { # if input IS phylo
    cleaned_names <- phy_new$tip.label
    ott_ids <- NULL
    # if we have ott_ids in phy, don't use_tnrs again:
    if (inherits(phy_new$ott_ids, "numeric") | inherits(phy_new$ott_ids, "integer")) {
      # if(!any(is.na(phy_new$ott_ids))){ #if there are no NAs
      use_tnrs_global <- FALSE
      ott_ids <- phy_new$ott_ids
      if (any(get_spp_from_taxon)) { # to use later, when get_spp_from_taxon = TRUE
        cleaned_names_tnrs <- list(ott_id = phy_new$ott_ids, unique_name = phy_new$tip.label)
      }
    }
  } else { # if input is NOT phylo, it can be a list
    input <- unlist(input) # in case input is given as a list
    # split elements by the commas:
    cleaned_names <- unlist(strsplit(input, ","))
    # clean split elements of lingering unneeded white spaces:
    cleaned_names <- stringr::str_trim(cleaned_names, side = "both")
    ott_ids <- NULL
  }
  if (use_tnrs_global) {
    # process names even if it's a "higher" taxon name:
    cleaned_names_tnrs <- clean_tnrs(tnrs_match(input = cleaned_names),
                                     remove_nonmatches = TRUE)
    # recover original names of invalid taxa and unmatched:
    ii <- !tolower(cleaned_names) %in% cleaned_names_tnrs$search_string
    cleaned_names <- c(cleaned_names_tnrs$unique_name, cleaned_names[ii])
    if (inherits(phy_new, "phylo")) {
      if (is.null(phy_new$ott_ids)) {
        cleaned_names <- gsub(" ", "_", cleaned_names)
        ii <- match(cleaned_names_tnrs$search_string, tolower(phy_new$tip.label))
        # after some tests, decided to use rotl's method instead of taxize::gnr_resolve, and just output the original input and the actual query for users to check out.
        # cleaned_names <- taxize::gnr_resolve(names = cleaned_names, data_source_ids=179, fields="all")$matched_name
        # rename the tip labels with tnrs matched names
        phy_new$tip.label[ii] <- cleaned_names[ii]
        ott_ids <- rep(NA, length(cleaned_names))
        ott_ids[ii] <- cleaned_names_tnrs$ott_id
        phy_new$ott_ids <- ott_ids
      }
    }
  }
  if (any(get_spp_from_taxon)) {
    if (length(get_spp_from_taxon) == 1) {
      get_spp_from_taxon <- rep(get_spp_from_taxon, length(cleaned_names))
    }
    if (length(cleaned_names) != length(get_spp_from_taxon)) {
      message("Specify all taxa in input to get species names from -")
      message("'input' and 'get_spp_from_taxon' arguments must have same length.")
      return(NA)
    }
    # rotl::tol_subtree is very fast but returns subspecies too \o/
    # it has no argument to restrict it to species only
    # so we are using our own function that wraps up their services nicely
    if ("ott" %in% taxonomic_source) {
      species_list <- get_opentree_species(taxon_name = cleaned_names_tnrs$unique_name,
                                           ott_id = cleaned_names_tnrs$ott_id,
                                           synth_tree_only = TRUE)
      cleaned_names <- species_list$species_names
      ott_ids <- species_list$ott_ids
    } else {
      # example: df <- get_ott_children(ott_ids = 698424, ott_rank = "species")
      df <- get_ott_children(ott_ids = cleaned_names_tnrs$ott_id, ott_rank = "species")
      # head(rownames(df[[1]])[grepl("species", df[[1]]$rank)])
      # the following does not work; it gives subspecies back
      # fixing it from get_ott_children function and here too
      cleaned_names <- lapply(df, function(x) rownames(x)[grepl("\\bspecies\\b", x$rank)])
      # enhance: create vector original_taxon with original names: rep(cleaned_names[i], length(cleaned_names[i]))
      original_taxa <- lapply(seq(nrow(cleaned_names_tnrs)), function(i) {
        rep(cleaned_names_tnrs$unique_name[i], length(cleaned_names[[i]]))
      })
      ott_ids <- lapply(df, function(x) x$ott_id[grepl("\\bspecies\\b", x$rank)])
      cleaned_names <- unlist(cleaned_names)
      ott_ids <- unlist(ott_ids)
   }
  }
  # Make sure that we are using underscores and not spaces:
  cleaned_names <- gsub("_", " ", cleaned_names)
  cleaned_names_print <- paste(utils::head(cleaned_names, 10), collapse = ", ")
  if (inherits(ott_ids, "numeric") | inherits(ott_ids, "integer")) {
    names(ott_ids) <- cleaned_names
  }
  # enhance: add original_taxa vector (from get_spp_from_taxon) to output here:
  datelife_query_return <- list(
    cleaned_names = cleaned_names, ott_ids = ott_ids,
    phy = phy_new
  )
  #TODO: print working taxa to a file instead of screen:
  message("Working with ",
          length(cleaned_names),
          " taxa: \n",
          cleaned_names_print,
          ifelse(length(cleaned_names) <=10, ".-", ", ..."))
  message("DateLife query done!\n")
  return(structure(datelife_query_return, class = "datelifeQuery"))
}
#' Process a phylo object or a character string to determine if it's correct newick
#'
#' @inheritParams make_datelife_query
#' @return A `phylo` object or `NA` if input is not a tree .
#' @export
input_process <- function(input) {
  message("... Phylo-processing 'input':")
  input_class <- "phylo"
  ott_ids <- NULL
  # TODO remove the multiPhylo if option from here?
  # make a method for input processing on multiPhylo objects
  if (inherits(input, "multiPhylo")) {
    message(message_multiphylo())
    input <- input[[1]]
  }
  if (inherits(input, "phylo")) {
    input_class <- class(input) # in case it has other classes to inherit
    ott_ids <- input$ott_ids # if phy already has ott ids, save them for later
    input <- ape::write.tree(input) # converts to newick
  }
  if (!is.character(input)) {
    message("'input' must be a character vector of taxon names, a newick string, or a 'phylo' or 'multiPhylo' object.")
    return(NA)
  }
  input <- gsub("\\+", " ", input) # clean the string
  input <- stringr::str_trim(input, side = "both") # clean the string
  phy_out <- NA
  if (any(grepl("\\(.*\\).*;", input))) { # our test for newick
    input <- input[grepl("\\(.*\\).*;", input)] # leave only the elements that are newick strings
    if (length(input) > 1) {
      message("'input' has several newick strings. Only the first one will be used...")
    }
    # phytools read.newick is not working
    # phy_out <- ape::collapse.singles(phytools::read.newick(text = gsub(" ", "_", input[1])))
    phy_out <- ape::collapse.singles(ape::read.tree(text = gsub(" ", "_", input[1])))
    phy_out$ott_ids <- ott_ids
    class(phy_out) <- input_class
    message("'input' is a phylogeny and it is correctly formatted...")
    # ape::read.tree creates NaN edge lengths for tree without branch lengths
    # clean it up:
    if (!is.null(phy_out$edge.length)) {
      if (any(is.na(phy_out$edge.length))) {
        warning("'input' has NA or NaN as branch lengths...")
        # phy_out$edge.length <- NULL
      }
    }
  } else {
    # not a requirement for input to be a phylogeny at this point
    message("'input' is not a phylogeny.") # so message instead of warning or stop
    return(NA)
  }
  return(phy_out)
}

#' Check if input is a `datelifeQuery` object
#'
#' @description `is_datelife_query` checks for two things to be `TRUE` or `FALSE`.
#' First, that `input` is of class {datelifeQuery}.
#' Second, that `input` is a list that contains at least two elements of a `datelifeQuery` object:
#' \describe{
#' 	 \item{cleaned_names}{A character vector of taxon names.}
#' 	 \item{phy}{Either NA or a `phylo` object.}
#' }
#' @param input An object to be checked as an object with essential properties of a 'datelifeQuery' object.
#' @return Is determined by the second condition.
#' @details If the object has the correct format but it has a class different than
#'  `datelifeQuery`, the class is not modified.
#' @export
is_datelife_query <- function(input) {
  if (is.list(input) & "phy" %in% names(input) & "cleaned_names" %in% names(input)) {
    if (inherits(input, "datelifeQuery")) {
      message("'input' is a 'datelifeQuery' object.")
    } else {
      message(
        "'input' has the elements of a 'datelifeQuery' object but is of class '",
        class(input),
        "'."
      )
      # class(input) <- "datelifeQuery"
    }
    return(TRUE)
  } else {
    message("'input' is not a 'datelifeQuery' object.")
    return(FALSE)
  }
}
