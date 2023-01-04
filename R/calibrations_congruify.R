# #' Congruify nodes of a tree topology to nodes from a source chronogram, and find the mrca nodes
# #' @inheritParams phylo_check
# #' @return An object of class data frame.
# #' @export
# congruify_and_mrca <- function(phy) {
#   UseMethod("congruify_and_mrca", phy)
# }

#' Congruify nodes of a tree topology to nodes from a source chronogram, and find the mrca nodes
#'
#' @description \code{congruify_and_mrca} congruifies a target tree against a single
#'   source chronogram, and gets nodes of target tree that correspond to the most
#'   recent common ancestor (mrca) of taxon pairs from the congruified calibrations.
#'    It uses [phytools::findMRCA()] to get mrca nodes.
#' @inheritParams phylo_check
# #' or a vector of taxon names (see details).
#' @param source_chronogram A `phylo` object, output of [datelife_search()].
#' @param reference A character string indicating the study reference that the `source_chronogram` comes from.
#' @return a `data.frame` of node ages from `source_chronograms` and corresponding
#' mrca nodes in target tree `phy`.
#' @export
congruify_and_mrca_phylo <- function(phy,
                                     source_chronogram,
                                     reference) {
    #
    if (!inherits(phy, "phylo")) {
      return(NA)
    }
    if (!inherits(source_chronogram, "phylo")) {
      return(NA)
    }
    if (missing(reference)) {
      reference <- "source_chronogram"
    }
    ############################################################################
    # homogenize tip labels:
    ############################################################################
    phy$tip.label <- sub(" ", "_", phy$tip.label)
    source_chronogram$tip.label <- sub(" ", "_", source_chronogram$tip.label)
    class(source_chronogram) <- "phylo"
    ############################################################################
    # congruify source chronogram to target topology:
    ############################################################################
    congruified <- suppressWarnings(geiger::congruify.phylo(
                                    reference = source_chronogram,
                                    target = phy,
                                    scale = NA,
                                    ncores = 1))
    if (!inherits(congruified, "list")) {
      warning("Congruification failed for study: ", reference)
      return(NA)
    }
    ############################################################################
    # get mrca nodes for available ages:
    ############################################################################
    mrcas <- mrca_calibrations(phy = congruified$target,
                               calibrations = congruified$calibrations)
    calibs_matched <- mrcas$matched_calibrations
    ############################################################################
    # add column of study reference:
    calibs_matched$reference <- rep(reference, nrow(calibs_matched))
    ############################################################################
    # reorder columns:
    # colnames are "MRCA" "MaxAge" "MinAge" "taxonA" "taxonB" "mrca_node_number"
    # "mrca_node_name" "reference"
    calibs_matched <- calibs_matched[ , c("mrca_node_name",
                                          "taxonA",
                                          "taxonB",
                                          "MinAge",
                                          "MaxAge",
                                          "reference",
                                          "mrca_node_number",
                                          "MRCA")]
    ############################################################################
    # order rows:
    row_order <- order(calibs_matched$mrca_node_number,
                    calibs_matched$MaxAge,
                    calibs_matched$reference,
                    calibs_matched$taxonA,
                    calibs_matched$taxonB)
    calibs_matched <- calibs_matched[row_order,]
    ############################################################################
    # add topology as attribute
    attributes(calibs_matched)$phy <- mrcas$matched_phy
    return(structure(calibs_matched,
                     class = c("congruifiedCalibrations", "data.frame")))
}

#' Congruify nodes of a tree topology to nodes from a source chronogram, and find the mrca nodes
#'
#' @description \code{congruify_and_mrca_multiPhylo} congruifies a target tree against all
#'   source chronograms in a `multiPhylo` object, and gets nodes of target tree
#'   that correspond to the most recent common ancestor (mrca) of taxon pairs
#'   in the congruified calibrations.
#'   It calls [congruify_and_mrca_phylo()], and [phytools::findMRCA()] to get mrca nodes.
#' @inheritParams phylo_check
# #' or a vector of taxon names (see details).
#' @param source_chronograms A `multiPhylo` object, output of [datelife_search()].
#' @return a `data.frame` of node ages from `source_chronograms` and corresponding
#' mrca nodes in target tree `phy`. `attributes(return)$phy` stores the congruified and mrca matched phylogeny.
#' @export
congruify_and_mrca_multiPhylo <- function(phy,
                                          source_chronograms) {
    ############################################################################
    # congruify and mrca each source chronogram:
    res <- lapply(seq(length(source_chronograms)),
                 function(i) {
                   xx <- congruify_and_mrca_phylo(phy = phy,
                                            source_chronogram = source_chronograms[[i]],
                                            reference = names(source_chronograms)[i])
                   return(xx)
                 })
    ############################################################################
    # get congruified and matched phylogenies, and remove null ones:
    matched_multiPhylo <- lapply(res, function(x) attributes(x)$phy)
    is_null <- sapply(matched_multiPhylo, is.null)
    matched_multiPhylo <- matched_multiPhylo[!is_null]
    ############################################################################
    # get one congruified and matched phy to assign as attribute later:
    phy <- matched_multiPhylo[[1]]
    ############################################################################
    # identify elements of res that are not tables:
    is_data_frame <- sapply(res, inherits, "data.frame")
    ############################################################################
    # merge the list of tables and return a single table (data.table object):
    res <- data.table::rbindlist(res[is_data_frame])
    ############################################################################
    # Convert to simple data.frame object:
    res <- as.data.frame.data.frame(res)
    ############################################################################
    # reorder columns:
    # colnames are "MRCA" "MaxAge" "MinAge" "taxonA" "taxonB" "mrca_node_number"
    # "mrca_node_name" "reference"
    res <- res[ , c("mrca_node_name",
                    "taxonA",
                    "taxonB",
                    "MinAge",
                    "MaxAge",
                    "reference",
                    "mrca_node_number",
                    "MRCA")]
    ############################################################################
    # order rows:
    row_order <- order(res$mrca_node_number,
                       res$MaxAge,
                       res$reference,
                       res$taxonA,
                       res$taxonB)
    res <- res[row_order,]
    ############################################################################
    # add topology as attribute
    attributes(res)$phy <- phy
    return(structure(res,
                     class = c("congruifiedCalibrations", "data.frame")))
}

#' Get summary statistics of ages in a `congruifiedCalibrations` object.
#'
#' @description Function `summarize_congruifiedCalibrations` returns a table of
#'   summary statistics for each node in `congruified_calibrations` argument.
#' @param congruified_calibrations A `congruifiedCalibrations` object, output of [congruify_and_mrca_multiPhylo()].
#' @param age_column A character string indicating the name of the column to be summarized.
#' @return A `data.frame` of summarized ages.
#' @export
summarize_congruifiedCalibrations <- function(congruified_calibrations,
                                              age_column) {
    ############################################################################
    if (missing(age_column)) {
      age_column = "MaxAge"
    }
    age_column <- age_column[1]
    if (!age_column %in% colnames(congruified_calibrations)) {
      stop("'age_column' provided is not a column of 'congruified_calibrations'.")
    }
    message("... Using ", age_column, " column to summarize ages.")
    if (all(congruified_calibrations$MinAge == congruified_calibrations$MaxAge)) {
      message("... Minimum and maximum ages are equal in 'congruified_calibrations' data.frame.")
    }
    ############################################################################
    # create vectors of summary statistics per node in congruified_calibrations
    min_ages <- q1 <- median_ages <- mean_ages <- q3 <- max_ages <- sd_ages <- var_ages <- c()
    for (node in unique(congruified_calibrations$mrca_node_name)) {
      rowsies <- congruified_calibrations$mrca_node_name %in% node
      min_ages <- c(min_ages, min(congruified_calibrations[rowsies, age_column]))
      q1 <- c(q1, stats::quantile(congruified_calibrations[rowsies, age_column], 0.25))
      median_ages <- c(median_ages, stats::median(congruified_calibrations[rowsies, age_column]))
      mean_ages <- c(mean_ages, mean(congruified_calibrations[rowsies, age_column]))
      q3 <- c(q3, stats::quantile(congruified_calibrations[rowsies, age_column], 0.75))
      max_ages <- c(max_ages, max(congruified_calibrations[rowsies, age_column]))
      sd_ages <- c(sd_ages, stats::sd(congruified_calibrations[rowsies, age_column]))
      var_ages <- c(var_ages, stats::var(congruified_calibrations[rowsies, age_column]))
    }
    ############################################################################
    # assemble summary table as data.frame
    summary_table <- data.frame(unique(congruified_calibrations$mrca_node_name),
                                min_ages,
                                q1,
                                median_ages,
                                mean_ages,
                                q3,
                                max_ages,
                                var_ages,
                                sd_ages)
    #
    ############################################################################
    # Give meaningful names to columns
    colnames(summary_table) <- c("Node Name",
                                 "Min Age",
                                 "Q1",
                                 "Median Age",
                                 "Mean Age",
                                 "Q3",
                                 "Max Age",
                                 "Variance",
                                 "SD")
    ############################################################################
    # TODO: add calibration_distribution as attribute
    # attributes(summary_table)$calibration_distribution <- NA
    return(summary_table)
}
#
# #' Congruify nodes of a tree topology to nodes from taxon pair calibrations
# #'
# #' @description \code{congruify_calibrations} get nodes of a tree topology given in
# #'   `phy` that correspond to the most recent common ancestor (mrca) of taxon
# #'   pairs given in `calibrations`. It uses [phytools::findMRCA()] to get mrca nodes.
# #'
# #' @inheritParams phylo_check
# # #' or a vector of taxon names (see details).
# #' @param chronograms A `phylo` or `multiPhylo` object, output of [datelife_search()].
# #' @param calibrations A `calibrations` object, an output of
# #'   [extract_calibrations_phylo()].
# #' @return A list of two elements:
# #' \describe{
# #' 	\item{matched_phy}{A `phylo` object with nodes renamed to match results of
# #'   the mrca search. Nodes are renamed using [tree_add_nodelabels()].}
# #' 	\item{matched_calibrations}{A `matchedCalibrations` object, which is the input `calibrations`
# #'    object with two additional columns storing results from the mrca search with
# #'    [phytools::findMRCA()]: `$mrca_node_number` and `$mrca_node_name`.}
# #' 	}
# #' @details The function takes pairs of taxon names in a calibrations data frame,
# #' and looks for them in the vector of tip labels of the tree. If both are present,
# #' then it gets the node that represents the most recent
# #' common ancestor (mrca) for that pair of taxa in the tree.
# #' Nodes of input `phy` can be named or not. They will be renamed.
# congruify_calibrations <- function(phy, chronograms, calibrations) {
#   ##############################################################################
#   ##############################################################################
#   # Extract all calibrations as a data.frame if none are provided
#   ##############################################################################
#   ##############################################################################
#   if (missing(calibrations)) {
#     extracted_calibrations <- extract_calibrations_phylo(input = chronograms,
#                                                          each = TRUE)
#   }
#   ##############################################################################
#   ##############################################################################
#   # Self congruify tree topology
#   # This Matches target tree nodelabels to calibration table
#   # so we have all taxon pairs and the corresponding real node labels
#   ##############################################################################
#   ##############################################################################
#   if (is.null(phy$edge.length)) {
#     message("... Adding temporary branch lengths to 'phy' for congruification.")
#     phy <- ape::compute.brlen(phy)
#   }
#   phy_calibrations <- suppressWarnings(
#                         geiger::congruify.phylo(reference = phy,
#                                                 target = phy,
#                                                 scale = NA,
#                                                 ncores = 1)$calibrations)
#   phy_matched <- mrca_calibrations(phy = phy,
#                                    calibrations = phy_calibrations)
#   #
#   ##############################################################################
#   ##############################################################################
#   # Congruify source chronograms to phy_matched$matched_phy and themselves
#   ##############################################################################
#   ##############################################################################
#   # Function
#   # Congruify a source chronogram with the target and itself
#   # return the calibrations data.frame only
#   run_congruification_self <- function(phy,
#                                   chronograms,
#                                   index) {
#     self <- chronograms[[index]]
#     phy$tip.label <- sub(" ", "_", phy$tip.label)
#     self$tip.label <- sub(" ", "_", self$tip.label)
#
#     cc_phy <- suppressWarnings(geiger::congruify.phylo(reference = phy,
#                                target = self,
#                                scale = NA,
#                                ncores = 1))
#     cc_self <- suppressWarnings(geiger::congruify.phylo(reference = self,
#                                target = self,
#                                scale = NA,
#                                ncores = 1))
#     list(phy = cc_phy, self = cc_self)
#   }
#   congruified_self <- lapply(seq(chronograms),
#                             function(i) {
#                               run_congruification_self(phy = phy_matched$matched_phy,
#                                                   chronograms,
#                                                   index = i)
#                             })
#   names(congruified_self) <- names(chronograms)
#
#   ##############################################################################
#   ##############################################################################
#   # Cross congruify source chronograms and target tree
#   ##############################################################################
#   ##############################################################################
#   fix_congruification <- function(congruified_self,
#                                   index = 1,
#                                   extracted_calibrations) {
#     source_chronogram <- congruified_self[[index]]$self$reference
#     target_tree <- congruified_self[[index]]$phy$reference
#     xx <- target_tree$node.label
#      # if number of congruified nodes in target tree is the same as node number from source chronogram
#     if (length(xx[ xx!= ""]) == source_chronogram$Nnode) {
#         # take extracted calibrations data frame as is
#         calibs <- extracted_calibrations[[index]]
#         calibs$congruent <- rep(TRUE, nrow(calibs))
#         calibs_matched <- mrca_calibrations(phy = source_chronogram,
#                                             calibrations = calibs)
#         # TODO: check that the following works:
#         # col_names <- colnames(calibs_matched$matched_calibrations)
#         # ii <- col_names %in% c("mrca_node_number", "mrca_node_name")
#         # col_names[ii] <- c("source_mrca_node_number", "source_mrca_node_name")
#         # colnames(calibs_matched$matched_calibrations) <- col_names
#       } else {
#         # case for index = 11, 10, 8, 7
#         # take congruified calibrations data frame
#         # TODO: needs checking, as the number of congruent nodes in target tree is higher than nodes available in output$calibrations, meaning, that some valid node data is being dropped from final results for some reason
#         # first, get most complete calibrations data.frame
#         # it is usually the one extracted from target tree
#         calibs <- congruified_self[[index]]$phy$calibrations
#         calibs_matched <- mrca_calibrations(phy = source_chronogram,
#                                             calibrations = calibs)
#         # use actual ages from source chronogram
#         mrca_node_name <- calibs_matched$matched_calibrations$mrca_node_name
#         source_bt <- ape::branching.times(calibs_matched$matched_phy)
#         source_ages <- source_bt[mrca_node_name]
#         calibs$MaxAge <- calibs$MinAge <- unname(source_ages)
#         # add congruent column
#         calibs$congruent <- rep(TRUE, nrow(calibs))
#         # add source node numbers and names
#         calibs$source_mrca_node_name <- mrca_node_name
#         calibs$source_mrca_node_number <- as.numeric(sub("n", "", mrca_node_name))
#         # add nodes that were not congruent
#         nn <- is.na(match(names(source_bt), mrca_node_name))
#         nn_ages <- ape::branching.times(calibs_matched$matched_phy)[nn]
#         nn_node_name <- names(nn_ages)
#         nn_node_number <- as.numeric(sub("n", "", nn_node_name))
#
#         non_congruent <- data.frame(MinAge = unname(nn_ages),
#                                     MaxAge = unname(nn_ages),
#                                     congruent = rep(FALSE, length(nn_ages)),
#                                     source_mrca_node_number = nn_node_number,
#                                     source_mrca_node_name = nn_node_name)
#         calibs <- data.table::rbindlist(list(calibs, non_congruent),
#                                         fill = TRUE)
#       }
#     calibs$reference <- rep(names(congruified_self)[index], nrow(calibs))
#     return(calibs)
#   }
#   # TODO: fix error in cross validation data set when running the following:
#   fixed_congruification <- lapply(seq(congruified_self),
#                                   function(i) {
#                                     tryCatch(fix_congruification(congruified_self,
#                                                         index = i,
#                                                         extracted_calibrations),
#                                               error = function(e) NA) } )
#   errored <- is.na(fixed_congruification)
#   if (any(errored)) {
#     message("fixing congruification errored for:\n",
#               paste(which(errored), "--", names(congruified_self)[errored], "\n"))
#   }
#   ##############################################################################
#   ##############################################################################
#   # Merge calibrations data frames into a single one
#   # and match to phy nodes
#   ##############################################################################
#   ##############################################################################
#   calibrations <-   data.table::rbindlist(fixed_congruification[!errored], fill = TRUE)
#   calibrations_results <- mrca_calibrations(phy = phy,
#                                             calibrations = calibrations)
#   return(structure(calibrations_results,
#                    class = "congruifiedCalibrations"))
# }
