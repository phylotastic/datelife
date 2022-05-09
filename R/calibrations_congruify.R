#' Congruify nodes of a tree topology to nodes from a source chronogram
#'
#' @description \code{congruify_and_mrca} congruifies a target tree against a source chronogram gets nodes of a tree topology given in
#'   `phy` that correspond to the most recent common ancestor (mrca) of taxon
#'   pairs from the congruified calibrations. It uses [phytools::findMRCA()] to get mrca nodes.
#' @inheritParams phylo_check
# #' or a vector of taxon names (see details).
#' @param source_chronogram A `phylo` object, output of [datelife_search()].
#' @param study A character string indicating the name of the study the `source_chronogram` comes from.
congruify_and_mrca <- function(phy,
                               source_chronogram,
                               study) {
    #
    if (!inherits(phy, "phylo")) {
      return(NA)
    }
    if (!inherits(source_chronogram, "phylo")) {
      return(NA)
    }
    if (missing(study)) {
      study <- "source_chronogram"
    }
    # homogenize tip labels:
    phy$tip.label <- sub(" ", "_", phy$tip.label)
    source_chronogram$tip.label <- sub(" ", "_", source_chronogram$tip.label)
    class(source_chronogram) <- "phylo"
    congruified <- suppressWarnings(geiger::congruify.phylo(reference = source_chronogram,
                               target = phy,
                               scale = NA,
                               ncores = 1))

    mrcas <- mrca_calibrations(phy = congruified$target,
                               calibrations = congruified$calibrations)
    calibs_matched <- mrcas$matched_calibrations
    calibs_matched$study <- rep(study, nrow(calibs_matched))
    return(structure(calibs_matched,
                     class = "congruifiedCalibrations"))
}

#' Congruify nodes of a tree topology to nodes from taxon pair calibrations
#'
#' @description \code{congruify_calibrations} get nodes of a tree topology given in
#'   `phy` that correspond to the most recent common ancestor (mrca) of taxon
#'   pairs given in `calibrations`. It uses [phytools::findMRCA()] to get mrca nodes.
#'
#' @inheritParams phylo_check
# #' or a vector of taxon names (see details).
#' @param chronograms A `phylo` or `multiPhylo` object, output of [datelife_search()].
#' @param calibrations A `calibrations` object, an output of
#'   [extract_calibrations_phylo()].
#' @return A list of two elements:
#' \describe{
#' 	\item{matched_phy}{A `phylo` object with nodes renamed to match results of
#'   the mrca search. Nodes are renamed using [tree_add_nodelabels()].}
#' 	\item{matched_calibrations}{A `matchedCalibrations` object, which is the input `calibrations`
#'    object with two additional columns storing results from the mrca search with
#'    [phytools::findMRCA()]: `$mrca_node_number` and `$mrca_node_name`.}
#' 	}
#' @details The function takes pairs of taxon names in a calibrations data frame,
#' and looks for them in the vector of tip labels of the tree. If both are present,
#' then it gets the node that represents the most recent
#' common ancestor (mrca) for that pair of taxa in the tree.
#' Nodes of input `phy` can be named or not. They will be renamed.
congruify_calibrations <- function(phy, chronograms, calibrations) {
  ##############################################################################
  ##############################################################################
  # Extract all calibrations as a data.frame if none are provided
  ##############################################################################
  ##############################################################################
  if (missing(calibrations)) {
    extracted_calibrations <- extract_calibrations_phylo(input = chronograms,
                                                         each = TRUE)
  }
  ##############################################################################
  ##############################################################################
  # Match target tree nodelabels to calibration table
  # so we have all taxon pairs and the corresponding real node labels
  ##############################################################################
  ##############################################################################
  if (is.null(phy$edge.length)) {
    message("... Adding temporary branch lengths to 'phy' for congruification.")
    phy <- ape::compute.brlen(phy)
  }
  phy_calibrations <- suppressWarnings(
                        geiger::congruify.phylo(reference = phy,
                                                target = phy,
                                                scale = NA,
                                                ncores = 1)$calibrations)
  phy_matched <- mrca_calibrations(phy = phy,
                                   calibrations = phy_calibrations)
  #
  ##############################################################################
  ##############################################################################
  # Congruify source chronograms to phy_matched$matched_phy and themselves
  ##############################################################################
  ##############################################################################
  # Function
  # Congruify a source chronogram with the target and itself
  # return the calibrations data.frame only
  run_congruification_self <- function(phy,
                                  chronograms,
                                  index) {
    self <- chronograms[[index]]
    phy$tip.label <- sub(" ", "_", phy$tip.label)
    self$tip.label <- sub(" ", "_", self$tip.label)

    cc_phy <- suppressWarnings(geiger::congruify.phylo(reference = phy,
                               target = self,
                               scale = NA,
                               ncores = 1))
    cc_self <- suppressWarnings(geiger::congruify.phylo(reference = self,
                               target = self,
                               scale = NA,
                               ncores = 1))
    list(phy = cc_phy, self = cc_self)
  }
  congruified_self <- lapply(seq(chronograms),
                            function(i) {
                              run_congruification_self(phy = phy_matched$matched_phy,
                                                  chronograms,
                                                  index = i)
                            })
  names(congruified_self) <- names(chronograms)

  ##############################################################################
  ##############################################################################
  # Cross congruify source chronograms and target tree
  ##############################################################################
  ##############################################################################
  fix_congruification <- function(congruified_self,
                                  index = 1,
                                  extracted_calibrations) {
    source_chronogram <- congruified_self[[index]]$self$reference
    target_tree <- congruified_self[[index]]$phy$reference
    xx <- target_tree$node.label
     # if number of congruified nodes in target tree is the same as node number from source chronogram
    if (length(xx[ xx!= ""]) == source_chronogram$Nnode) {
        # take extracted calibrations data frame as is
        calibs <- extracted_calibrations[[index]]
        calibs$congruent <- rep(TRUE, nrow(calibs))
        calibs_matched <- mrca_calibrations(phy = source_chronogram,
                                            calibrations = calibs)
        # TODO: check that the following works:
        # col_names <- colnames(calibs_matched$matched_calibrations)
        # ii <- col_names %in% c("mrca_node_number", "mrca_node_name")
        # col_names[ii] <- c("source_mrca_node_number", "source_mrca_node_name")
        # colnames(calibs_matched$matched_calibrations) <- col_names
      } else {
        # case for index = 11, 10, 8, 7
        # take congruified calibrations data frame
        # TODO: needs checking, as the number of congruent nodes in target tree is higher than nodes available in output$calibrations, meaning, that some valid node data is being dropped from final results for some reason
        # first, get most complete calibrations data.frame
        # it is usually the one extracted from target tree
        calibs <- congruified_self[[index]]$phy$calibrations
        calibs_matched <- mrca_calibrations(phy = source_chronogram,
                                            calibrations = calibs)
        # use actual ages from source chronogram
        mrca_node_name <- calibs_matched$matched_calibrations$mrca_node_name
        source_bt <- ape::branching.times(calibs_matched$matched_phy)
        source_ages <- source_bt[mrca_node_name]
        calibs$MaxAge <- calibs$MinAge <- unname(source_ages)
        # add congruent column
        calibs$congruent <- rep(TRUE, nrow(calibs))
        # add source node numbers and names
        calibs$source_mrca_node_name <- mrca_node_name
        calibs$source_mrca_node_number <- as.numeric(sub("n", "", mrca_node_name))
        # add nodes that were not congruent
        nn <- is.na(match(names(source_bt), mrca_node_name))
        nn_ages <- ape::branching.times(calibs_matched$matched_phy)[nn]
        nn_node_name <- names(nn_ages)
        nn_node_number <- as.numeric(sub("n", "", nn_node_name))

        non_congruent <- data.frame(MinAge = unname(nn_ages),
                                    MaxAge = unname(nn_ages),
                                    congruent = rep(FALSE, length(nn_ages)),
                                    source_mrca_node_number = nn_node_number,
                                    source_mrca_node_name = nn_node_name)
        calibs <- data.table::rbindlist(list(calibs, non_congruent),
                                        fill = TRUE)
      }
    calibs$reference <- rep(names(congruified_self)[index], nrow(calibs))
    return(calibs)
  }
  # TODO: fix error in cross validation data set when running the following:
  fixed_congruification <- lapply(seq(congruified_self),
                                  function(i) {
                                    tryCatch(fix_congruification(congruified_self,
                                                        index = i,
                                                        extracted_calibrations),
                                              error = function(e) NA) } )
  errored <- is.na(fixed_congruification)
  if (any(errored)) {
    message("fixing congruification errored for:\n",
              paste(which(errored), "--", names(congruified_self)[errored], "\n"))
  }
  ##############################################################################
  ##############################################################################
  # Merge calibrations data frames into a single one
  # and match to phy nodes
  ##############################################################################
  ##############################################################################
  calibrations <-   data.table::rbindlist(fixed_congruification[!errored], fill = TRUE)
  calibrations_results <- mrca_calibrations(phy = phy,
                                            calibrations = calibrations)
  return(structure(calibrations_results,
                   class = "congruifiedCalibrations"))
}
