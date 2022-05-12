#' Use calibrations to date a topology with the BLADJ algorithm.
#'
#' @description The function `use_calibrations_bladj` prepares the input for BLADJ
#' and calls [make_bladj_tree()].
#' @param phy A `phylo` object with or without branch lengths.
#' @param calibrations A `data.frame` of secondary calibrations for any pair of taxon
#' names in `phy`, usually obtained with [get_all_calibrations()].
#' @param type The type of age to use as calibration. Options are "median", "mean", "min", or "max".
#' @param root_age Not implemented yet. Numeric specifying the age of the root.
#' If there are no calibrations for it. If NULL or not numeric, the maximum
#' calibration plus a unit of the mean differences will be used as root calibration.
#' If there is only one internal calibration, the root age will be set to 10% more
#' than the age of the calibration.
#' @param match_calibrations Boolean, default to TRUE. It will run
#' [match_all_calibrations()] function. Set to FALSE if your ages have
#' already been matched.
#' @return A `phylo` object with branch lengths proportional to time.
#' @details
#' The BLADJ algorithm is part of the Phylocom software, presented in Webb et al.
#' (2008) \doi{10.1093/bioinformatics/btn358}.
#' @references
#' Webb, C. O., Ackerly, D. D., & Kembel, S. W. (2008). "Phylocom: software for
#' the analysis of phylogenetic community structure and trait evolution".
#' Bioinformatics, 24(18), \doi{10.1093/bioinformatics/btn358}.

#' @export
use_calibrations_bladj <- function(phy = NULL,
                                   calibrations,
                                   type = "median",
                                   root_age = NULL,
                                   match_calibrations = TRUE,
                                   congruified = FALSE) {
  #
  message("... Using ages as secondary calibrations with BLADJ.")
  if (inherits(calibrations, "matchedCalibrations")) {
      phy <- calibrations$phy
      calibs <- calibrations$matched_calibrations
    stop("Not implemented yet, apologies.")
  }
  if (inherits(calibrations, "congruifiedCalibrations")) {
    phy <- attributes(calibrations)$phy

  }
  type <- match.arg(tolower(type), c("mean", "min", "max", "median"))

  if (nrow(calibs$matched_calibrations) < 1) {
    stop("Nodes in 'calibrations' (determined by taxon pairs) do not match any nodes
            in 'phy'.\n\t Dating analysis is not possible with this set of calibrations.")
  }
  age_distributions <- calibs$phy$calibration_distribution
  node_names <- calibs$matched_calibrations$NodeNames


  if ("mean" %in% type) {
    node_ages <- sapply(age_distributions, mean)
    # length(node_ages)
  }
  if ("min" %in% type) {
    node_ages <- sapply(age_distributions, min)
  }
  if ("max" %in% type) {
    node_ages <- sapply(age_distributions, max)
  }
  if ("median" %in% type) {
    node_ages <- sapply(age_distributions, stats::median)
  }
  # check that the root node has a calibration
  # otherwise, add one
  # bladj stops if there is no calibration for the root
  if (!"n1" %in% node_names) {
    if (length(node_ages) > 1) {
      # if there is only one calibration the line below gives give NaN
      root_age <- max(node_ages) + mean(abs(diff(sort(node_ages))))
    } else {
      # case when there is only one calibration:
      root_age <- 1.1 * max(node_ages)
    }
    node_ages <- c(root_age, node_ages)
    node_names <- c("n1", node_names)
  }
  # bladj freezes if calibrations are in conflict
  chronogram <- make_bladj_tree(
    tree = calibs$phy, nodeages = node_ages,
    nodenames = node_names
  )
  # TODO: something about these extra list elements, set up as attributes??
  chronogram$dating_method <- "bladj"
  chronogram$calibration_distribution <- ages_distribution
  chronogram$calibrations <- calibs$matched_calibrations
  chronogram$used_calibrations <- stats::setNames(node_ages, node_names)
  return(chronogram)
}


#' Check for conflicting calibrations.
#'
#' `check_conflicting_calibrations` checks if calibrations are younger or older
#'  relative to descendants and ancestors, respectively.
#'
#' @details It removes conflicting calibrations if needed, but BLADJ works as long as it has an age for the root.
#' @inheritParams phylo_check
#' @param calibration_distribution A list of node age distributions, named with `phy`'s node numbers.
check_conflicting_calibrations <- function(phy,
                                           calibration_distribution) {
  #
  message("... Checking for conflicting calibrations.")
  if (!inherits(phy, "phylo")) {
    message("'phy' is not a phylo object")
    return(NA)
  }
  des <- lapply(as.numeric(names(calibration_distribution)), function(i) phy$edge[phy$edge[, 1] == i, 2])
  anc <- lapply(as.numeric(names(calibration_distribution)), function(i) phy$edge[phy$edge[, 2] == i, 1])
  foo1 <- function(node, des) {
    des <- des[!is.na(des)]
    any(node < des)
  }
  foo2 <- function(node, anc) {
    anc <- anc[!is.na(anc)]
    any(node > anc)
  }
  calibration_distribution_des <- sapply(des, function(x) unlist(calibration_distribution[as.character(x)[as.character(x) %in% names(calibration_distribution)]]))
  calibration_distribution_anc <- lapply(anc, function(x) unlist(calibration_distribution[as.character(x)[as.character(x) %in% names(calibration_distribution)]]))
  node_younger <- mapply(foo1, calibration_distribution, calibration_distribution_des) # which calibration_distribution are younger than their descendants
  node_older <- mapply(foo2, calibration_distribution, calibration_distribution_anc) # which calibration_distribution are older than their ancestor
  list(des, anc, calibration_distribution_des, calibration_distribution_anc, node_younger, node_older)
}

# fix_nodeages_for_bladj <- function(phy, calibration_distribution, remove = TRUE, remove_which = "younger", expand = FALSE){
#   remove <- match.arg(tolower(remove), c("older", "younger"))
#   calibration_distribution_original <- calibration_distribution
#   cc <- check_conflicting_calibrations(phy, calibration_distribution)
#   if("older" %in% remove){
#     calibration_distribution <- calibration_distribution[!cc$node_older]
#   }
#   if("younger" %in% remove){
#     calibration_distribution <- calibration_distribution[!cc$node_younger]
#   }
#   # the root needs to be fixed in bladj always,
#   # so if it is the root that is younger than any of its descendants, we need to fix that
#   # we will just add the difference:
#   if(cc$node_younger[1]){
#     calibration_distribution[1] <- max(cc$calibration_distribution_des[[1]], na.rm = TRUE) +
#       abs(max(cc$calibration_distribution_des[[1]], na.rm = TRUE) - calibration_distribution[1])
#     cc$node_younger[1] <- FALSE
#   }
#   calibration_distribution
# }
