#' Get a taxon summary of a `datelifeResult` object.
#'
# #' To be renamed `summary_taxon`.
#'
#' @param datelife_result A `datelifeResult` object, usually an output of [get_datelife_result()].
#' @param datelife_query A `datelifeQuery` object, usually an output of [make_datelife_query()].
#' @return A `datelifeTaxonSummary` object, which is a list of 4 elements:
#' \describe{
#'    \item{$matrix}{Data as a presence/absence matrix of taxon names across chronograms.}
#'    \item{$summary}{A `data.frame` with taxon names as [row.names()] and two
#'     columns, one with the number of chronograms that contain a taxon name and
#'     the other one with the total number of chronograms that have at least 2
#'     taxon names.}
#'    \item{$summary2}{A `data.frame` with chronogram citations as [row.names()]
#'     and two columns, one with the number of taxon names found in each chronogram
#'     and the other one with the total number of taxon names.}
#'    \item{$absent_taxa}{A character vector of taxon names that are not found
#'     in the chronogram database.}
#'}
#' @export
get_taxon_summary <- function(datelife_result = NULL,
                              datelife_query = NULL) {

  # datelife_result <- check_datelife_result(datelife_result)
  if (is.null(datelife_result) | !inherits(datelife_result, "datelifeResult")) {
    warning("'datelife_result' argument must be a list of patristic matrices (you can get one with get_datelife_result()).")
    return(NA)
  }

  if (suppressMessages(is_datelife_query(input = datelife_query))) {
    if (is.null(attributes(datelife_result)$datelife_query)) {
      cleaned_names <- datelife_query$cleaned_names
    } else {
      input <- attributes(datelife_result)$datelife_query
      cleaned_names <- attributes(datelife_result)$datelife_query$cleaned_names
    }
  } else {
    message("'datelife_query' argument was not provided.")
    message("Taxa absent in all chronograms are not reported.")
    cleaned_names <- unique(rapply(datelife_result, rownames))
  }
  # results.index <- datelife_result_study_index(datelife_result, cache)
  return.object <- NA
  input.match <- unique(rapply(datelife_result, rownames))
  # if(any(!input.match %in% input)) warning("input does not contain all or any taxa from filteredresults object")
  absent.input <- cleaned_names[!cleaned_names %in% input.match]
  if (length(absent.input) <= 0) {
    if (is.null(datelife_query)) {
      absent.input <- NA # we cannot know if there are any taxa that are completely absent, because original query was not provided
    } else {
      absent.input <- "None" # we know there are no absent taxa
    }
  }
  taxon_list <- vector(mode = "list")
  # tax <- unique(rapply(datelife_result, rownames)) #rownames(datelife_result[[1]])
  for (result.index in sequence(length(datelife_result))) {
    n <- rownames(datelife_result[[result.index]])
    m <- match(input.match, n)
    taxon_list[[result.index]] <- n[m]
  }
  taxon_matrix <- do.call(rbind, taxon_list) # transforms a list of names into a matrix of names
  taxon_matrix <- !is.na(taxon_matrix) # makes a boolean matrix
  colnames(taxon_matrix) <- input.match
  rownames(taxon_matrix) <- paste0("Chronogram", sequence(nrow(taxon_matrix)))
  chronogram_names <- names(datelife_result)
  names(chronogram_names) <- rownames(taxon_matrix)
  # from here, replace by print: take it to a print function, so we only store the matrix, absent taxa and chronogram names
  # tax <- unique(rapply(datelife_result, rownames)) #rownames(datelife_result[[1]])
  x <- rapply(datelife_result, rownames)
  # get the proportion of chronograms where a taxon is found
  prop <- c()
  for (taxon in input.match) {
    prop <- c(prop, paste0(length(which(taxon == x)), "/", length(datelife_result)))
  }
  taxon_summary <- data.frame(taxon = input.match, chronograms = prop)
  taxon_number <- sapply(seq(nrow(taxon_matrix)), function(x) sum(taxon_matrix[x, ]))
  taxon_summary2 <- data.frame(
    chronogram = names(datelife_result),
    taxon_number = taxon_number, total_taxa = rep(ncol(taxon_matrix), length(datelife_result))
  )
  res <- list(
    matrix = taxon_matrix, summary = taxon_summary, summary2 = taxon_summary2,
    absent_taxa = unique(absent.input)
  )
  # end replace by print
  # res <- list(matrix = taxon_matrix, chronogram_names = chronogram_names, absent_taxa = unique(absent.input))
  class(res) <- "datelifeTaxonSummary"
  return(res)
}
# print.datelifeTaxonSummary <- function(taxon_summary){
#
# }
#' Summarize a `datelifeResult` object.
#'
#' @description Get different types of summaries from a `datelifeResult`
#' object, an output from [get_datelife_result()].
#' This allows rapid processing of data.
#' If you need a list of chronograms from your `datelifeResult` object, this
#' is the function you are looking for.
#'
#' @inheritParams get_taxon_summary
#' @param summary_format A character vector of length one, indicating the output
#' format for results of the DateLife search. Available output formats are:
#' \describe{
#' 	 \item{"citations"}{A character vector of references where chronograms with
#' 	 				some or all of the target taxa are published (source chronograms).}
#' 	 \item{"mrca"}{A named numeric vector of most recent common ancestor (mrca)
#' 	 				ages of target taxa defined in input, obtained from the source chronograms.
#' 	 				Names of mrca vector are equal to citations.}
#' 	 \item{"newick_all"}{A named character vector of newick strings corresponding
#' 	 				to target chronograms derived from source chronograms. Names of newick_all
#' 	 				vector are equal to citations.}
#' 	 \item{"newick_sdm"}{Only if multiple source chronograms are available. A
#' 	 				character vector with a single newick string corresponding to a target
#' 	 				chronogram obtained with SDM supertree method (Criscuolo et al. 2006).}
#' 	 \item{"newick_median"}{Only if multiple source chronograms are available.
#' 	 				A character vector with a single newick string corresponding to a target
#' 	 				chronogram from the median of all source chronograms.}
#' 	 \item{"phylo_sdm"}{Only if multiple source chronograms are available. A
#' 	 				phylo object with a single target chronogram obtained with SDM supertree
#' 	 				method (Criscuolo et al. 2006).}
#' 	 \item{"phylo_median"}{Only if multiple source chronograms are available. A
#' 	 				phylo object with a single target chronogram obtained from source
#' 	 				chronograms with median method.}
#' 	 \item{"phylo_all"}{A named list of phylo objects corresponding to each target
#' 	 				chronogram obtained from available source chronograms. Names of
#' 	 				phylo_all list correspond to citations.}
#' 	 \item{"phylo_biggest"}{The chronogram with the most taxa. In the case of a
#' 	 				tie, the chronogram with clade age closest to the median age of the
#' 	 				equally large trees is returned.}
#' 	 \item{"html"}{A character vector with an html string that can be saved and
#' 	 				then opened in any web browser. It contains a 4 column table with data on
#' 	 				target taxa: mrca, number of taxa, citations of source chronogram and
#' 	 				newick target chronogram.}
#' 	 \item{"data_frame"}{A 4 column `data.frame` with data on target taxa: mrca, number of
#' 	 				taxa, citations of source chronograms and newick string.}
#' }
#' @inheritParams datelife_result_MRCA
# datelife_result_MRCA has param na_rm
#' @param summary_print A character vector specifying the type of summary information
#'   to be printed to screen. Options are:
#'   \describe{
#'   	 \item{"citations"}{Prints references of chronograms where target taxa are found.}
#'   	 \item{"taxa"}{Prints a summary of the number of chronograms where each target
#'   	 				taxon is found.}
#'   	 \item{"none"}{Nothing is printed to screen.}
#'   }
#'   Defaults to `c("citations", "taxa")`, which displays both.
#' @param taxon_summary A character vector specifying if data on target taxa missing
#'   in source chronograms should be added to the output as a `"summary"` or as a
#'   presence/absence `"matrix"`. Default to `"none"`, no information on taxon_summary
#'   added to the output.
#' @param criterion Defaults to `criterion = "taxa"`. Used for chronogram summarizing, i.e., obtaining a single
#'   summary chronogram from a group of input chronograms.
#'   For summarizing approaches that return a single summary tree from a group of
#'   phylogenetic trees, it is necessary that the latter form a grove, roughly,
#'   a sufficiently overlapping set of taxa between trees, see Ané et al. (2009) \doi{10.1007/s00026-009-0017-x}.
#'   In rare cases, a group of trees can have multiple groves. This argument indicates
#'   whether to get the grove with the most trees (`criterion = "trees"`) or the
#'   most taxa (`criterion = "taxa"`).
#' @return The output is determined by the argument `summary_format`:
#' \describe{
#'   \item{If `summary_format = "citations"`}{The function returns a character
#'     vector of references.}
#'   \item{If `summary_format = "mrca"`}{The function returns a named numeric
#'     vector of most recent common ancestor (mrca) ages.}
#'   \item{If `summary_format = "newick_[all, sdm, or median]"`}{The function
#'     returns output chronograms as newick strings.}
#'   \item{If `summary_format = "phylo_[all, sdm, median, or biggest]"`}{The
#'     function returns output chronograms as `phylo` or `multiPhylo` objects.}
#'   \item{If `summary_format = "html" or "data_frame"`}{The function returns a
#'     4 column table with data on mrca ages, number of taxa, references, and output chronograms as newick strings.}
#' }
#' @references
#' Ané, C., Eulenstein, O., Piaggio-Talice, R., & Sanderson, M. J. (2009).
#' "Groves of phylogenetic trees". Annals of Combinatorics, 13(2), 139-167,
#' \doi{10.1007/s00026-009-0017-x}.
#' @export
summarize_datelife_result <- function(datelife_result = NULL,
                                      datelife_query = NULL,
                                      summary_format = "phylo_all",
                                      na_rm = TRUE,
                                      summary_print = c("citations", "taxa"),
                                      taxon_summary = c("none", "summary", "matrix"),
                                      criterion = "taxa") {
  taxon_summ <- get_taxon_summary(
    datelife_result = datelife_result,
    datelife_query = datelife_query
  )
  if (length(taxon_summ) == 1) {
    message("get_taxon_summary failed.")
    return(NA)
  }
  summary_format.in <- match.arg(summary_format, choices = c(
    "citations",
    "mrca",
    "newick_all",
    "newick_sdm",
    "newick_median",
    "phylo_sdm",
    "phylo_median",
    "phylo_biggest",
    "phylo_all",
    "html",
    "data_frame"
  ))
  taxon_summary.in <- match.arg(taxon_summary, choices = c("none", "summary", "matrix"))
  summary_print.in <- match.arg(summary_print, c("citations", "taxa", "none"), several.ok = TRUE)
  if (summary_format.in == "citations") {
    return.object <- names(datelife_result)
  }
  mrcas <- datelife_result_MRCA(datelife_result, na_rm = na_rm) # this is later used for median and sdm
  if (summary_format.in == "mrca") {
    return.object <- mrcas
  }
  # TODO replace by datelife_result_newick_all
  if (summary_format.in == "newick_all") {
    trees <- sapply(datelife_result, patristic_matrix_to_newick)
    return.object <- trees[which(!is.na(trees))]
    newick_all <- return.object
  }
  # TODO replace by datelife_result_phylo_all
  if (summary_format.in == "phylo_all") {
    trees <- suppressWarnings(lapply(datelife_result, patristic_matrix_to_phylo)) # suppress warning "Converting from patristic distance matrix to a tree resulted in some negative branch lengths"
    return.object <- trees[which(!is.na(trees))]
    class(return.object) <- "multiPhylo"
    phylo_all <- return.object
    phylo_biggest <- get_biggest_multiphylo(trees) # NAs in trees are removed in get_biggest_multiphylo
  }
  if (summary_format.in == "phylo_biggest") {
    trees <- lapply(datelife_result, patristic_matrix_to_phylo)
    return.object <- get_biggest_multiphylo(trees) # NAs in trees are removed in get_biggest_multiphylo
    # TODO: use get_biggest_datelife_result instead?
  }
  # test if n_overlap = 2 is enough to summarize results with sdm and median:
  if (summary_format.in %in% c("newick_sdm", "phylo_sdm", "newick_median", "phylo_median")) {
    best_grove <- get_best_grove(datelife_result, criterion = criterion, n = 2)$best_grove
  }
  if (inherits(datelife_query, "datelifeQuery")) {
    target_tree <- datelife_query$phy
  } else {
    target_tree <- NULL
  }
  if (grepl("median", summary_format.in)) {
    return.object <- datelife_result_median(best_grove, target_tree = target_tree)
  }
  if (grepl("sdm", summary_format.in)) {
    return.object <- datelife_result_sdm_phylo(best_grove, target_tree = target_tree)
  }
  if (summary_format.in %in% c("newick_sdm", "newick_median")) {
    return.object <- ape::write.tree(return.object)
  }
  # if(summary_format.in %in% c("phylo_sdm", "phylo_median")) {
  # 	return.object <- tree
  # }
  if (summary_format.in == "html") {
    out.vector1 <- "<table border='1'><tr><th>MRCA Age (MY)</th><th>Ntax</th><th>Citation</th><th>Newick"
    if (taxon_summary.in == "matrix") {
      out.vector1 <- paste(out.vector1, paste("</th><th>", colnames(taxon_summ$matrix), sep = "", collapse = ""), sep = "")
    }
    out.vector1 <- paste(out.vector1, "</th></tr>", sep = "")
    ages <- datelife_result_MRCA(datelife_result, na_rm = na_rm)
    trees <- sapply(datelife_result, patristic_matrix_to_newick)
    out.vector2 <- c()
    for (result.index in sequence(length(datelife_result))) {
      out.vector2 <- paste(out.vector2, "<tr><td>", ages[result.index], "</td><td>", sum(!is.na(diag(datelife_result[[result.index]]))), "</td><td>", names(datelife_result)[result.index], "</td><td>", trees[result.index], sep = "")
      if (taxon_summary.in == "matrix") {
        out.vector2 <- paste(out.vector2, paste("</td><td>", taxon_summ$matrix[result.index, ], sep = "", collapse = ""), sep = "")
      }
      out.vector2 <- paste(out.vector2, "</td></tr>", sep = "")
    }
    out.vector <- paste(out.vector1, out.vector2, "</table>", sep = "")
    if (taxon_summary.in == "summary") {
      taxon_summ$summary.html <- as.matrix(taxon_summ$summary)
      out.vector3 <- "<p></p><table border='1'><tr><th>taxon</th><th>chronograms</th><tr>"
      for (summary.index in sequence(nrow(taxon_summ$summary.html))) {
        out.vector3 <- paste(out.vector3, paste("</td><td>", taxon_summ$summary.html[summary.index, ], sep = "", collapse = ""), "</td></tr>", sep = "")
      }
      out.vector <- paste(out.vector, out.vector3, "</table>", sep = "")
    }
    if (taxon_summary.in != "none") {
      out.vector4 <- "<p></p><table border='1'><tr><th> </th><th>Absent Taxa</th><tr>"
      for (i in 1:length(taxon_summ$absent_taxa)) {
        out.vector4 <- paste(out.vector4, "<tr><td>", i, "</td><td>", taxon_summ$absent_taxa[i], "</td><tr>", sep = "")
      }
      out.vector4 <- paste(out.vector4, "</table>", sep = "")
      out.vector <- paste(out.vector, out.vector4, sep = "")
    }
    return.object <- out.vector
  }
  if (summary_format.in == "data_frame") {
    out.df <- data.frame()
    ages <- datelife_result_MRCA(datelife_result, na_rm = na_rm)
    trees <- sapply(datelife_result, patristic_matrix_to_newick)
    for (result.index in sequence(length(datelife_result))) {
      out.line <- data.frame(Age = ages[result.index], Ntax = sum(!is.na(diag(datelife_result[[result.index]]))), Citation = names(datelife_result)[result.index], Newick = trees[result.index])
      if (result.index == 1) {
        out.df <- out.line
      } else {
        out.df <- rbind(out.df, out.line)
      }
    }
    if (taxon_summary.in == "matrix") {
      out.df <- cbind(out.df, taxon_summ$matrix)
    }
    rownames(out.df) <- NULL
    return.object <- out.df
  }
  if (taxon_summary.in != "none" & summary_format.in != "html") {
    return.object <- list(return.object)
    if (taxon_summary.in == "matrix") {
      if (summary_format.in != "data_frame") {
        return.object <- c(return.object, list(taxon_distribution = taxon_summ$matrix))
      }
    }
    if (taxon_summary.in == "summary") {
      return.object <- c(return.object, list(taxon_distribution = taxon_summ$summary))
    }
    return.object <- c(return.object, list(absent_taxa = data.frame(taxon = taxon_summ$absent_taxa)))
    names(return.object)[1] <- summary_format.in
    if (summary_format.in == "data_frame") {
      names(return.object)[1] <- "results"
      if (taxon_summary.in == "matrix") {
        names(return.object)[1] <- "results_and_missing_taxa"
      }
    }
  }
  if (any("citations" %in% summary_print.in) & !any(summary_format.in %in% c("citations", "html", "data_frame"))) {
    if (summary_format.in == "citations") {
      message("Input taxa found in trees from:")
    } else {
      message("Source chronograms from:", "\n")
    }
    for (i in 1:length(datelife_result)) {
      message(i, ": ", names(datelife_result)[i], "\n")
    }
  }
  if (any(grepl("taxa", summary_print.in)) & taxon_summary.in != "summary") {
    message("Input taxa presence across source chronograms:")
    message(paste0(utils::capture.output(taxon_summ$summary), collapse = "\n"), "\n")
    message("Input taxa completely absent from source chronograms:")
    message(paste0(utils::capture.output(data.frame(taxon = taxon_summ$absent_taxa)), collapse = "\n"), "\n")
  }
  return(return.object)
}


#' Get the tree with the most tips from a multiPhylo object: the biggest tree.
#'
#' @param trees A list of trees as `multiPhylo` or as a generic `list` object.
#' @return The largest tree from those given in `trees`, as a `phylo` object with an additional `$citation` element containing the reference of the original publication.
#' @export
get_biggest_multiphylo <- function(trees) {
  trees <- trees[which(!is.na(trees))] # removes NAs, which will return an error later on next logical:
  tree_citation <- names(trees)
  return.object <- trees[which(sapply(trees, ape::Ntip) == max(sapply(trees, ape::Ntip)))]
  tree_citation <- tree_citation[which(sapply(trees, ape::Ntip) == max(sapply(trees, ape::Ntip)))]
  if (length(return.object) > 1) { # there are more than one tree with same number of taxa. Rather than take the first by default, take the one with the most intermediate depth (this assumes that the root node is the same for all trees). An example is the Bininda-Emonds et al. mammal tree: there are three trees with min, max, and best guess calibrations. So, take the one in the middle.
    max.branching.time <- function(x) {
      return(max(ape::branching.times(x)))
    }
    tree.depths <- sapply(return.object, max.branching.time)
    return.object <- return.object[which.min(abs(tree.depths - stats::median(tree.depths)))]
    tree_citation <- tree_citation[which.min(abs(tree.depths - stats::median(tree.depths)))]
  }
  if (!inherits(return.object, "phylo")) {
    return.object <- return.object[[1]]
  }
  return.object$citation <- tree_citation
  return.object
}

get_biggest_phylo <- get_biggest_multiphylo

# TODO: choose from patristic matrix, not phylo objects, it will be faster.
# get_biggest_datelife_result
