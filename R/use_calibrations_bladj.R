#' Use calibrations to date a topology with the BLADJ algorithm.
#'
#' @description The function `use_calibrations_bladj` prepares the input for BLADJ
#' and calls [make_bladj_tree()].
#' @param phy A `phylo` object with or without branch lengths.
#' @param calibrations A `data.frame` of secondary calibrations for any pair of taxon
#' names in `phy`, usually obtained with [get_all_calibrations()].
#' @param type The type of age to use as calibration. Options are "median", "mean", "min", or "max".
#' @param root_age Numeric specifying the age of the root. Only used if there are
#' no ages for the root node in  `calibrations` argument.
#' If missing, NULL, or not numeric, the value of the oldest calibration plus one
#' unit of the mean differences across calibrations, will be used as root calibration.
#' If there is one single age point provided as `calibrations`, the root age will
#' be set to 10% more than the age of the single calibration.
#' @return A chronogram: a `phylo` object with branch lengths proportional to time.
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
                                   root_age) {
  ############################################################################
  # initial checks
  type <- match.arg(tolower(type), c("mean", "min", "max", "median"))
  message("... Using ", type, " ages as secondary calibrations with BLADJ.")
  ############################################################################
  if (inherits(calibrations, "matchedCalibrations")) {
      phy <- calibrations$phy
      age_distributions <- calibrations$phy$calibration_distribution
      node_names <- calibrations$matched_calibrations$NodeNames
      calibrations <- calibrations$matched_calibrations
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
    # stop("Not finished yet, apologies.")
  }
  ############################################################################
  if (inherits(calibrations, "congruifiedCalibrations")) {
    message("...'calibrations' is a 'congruifiedCalibrations' object.")
    phy <- attributes(calibrations)$phy
    summary_ages <- summarize_congruifiedCalibrations(calibrations)
    node_names <- summary_ages$`Node Name`
    if ("mean" %in% type) {
      node_ages <- summary_ages$`Mean Age`
      # length(node_ages)
    }
    if ("min" %in% type) {
      node_ages <- summary_ages$`Min Age`
    }
    if ("max" %in% type) {
      node_ages <- summary_ages$`Max Age`
    }
    if ("median" %in% type) {
      node_ages <- summary_ages$`Median Age`
    }
  }
  ############################################################################
  if (nrow(calibrations) < 1) {
    stop("Nodes in 'calibrations' (determined by taxon pairs) do not match any nodes
            in 'phy'.\n\t Dating analysis is not possible with this set of calibrations.")
  }
  ############################################################################
  # check that the root node has a calibration
  # otherwise, add one because
  # bladj stops if there is no calibration for the root
  # first, get the node label of the root:
  root_node_name <- min(phy$node.label)
  if (!root_node_name %in% node_names) {
    if (missing(root_age)) {
      if (length(node_ages) > 1) {
        # if there is only one calibration the line below gives give NaN
        root_age <- max(node_ages) + mean(abs(diff(sort(node_ages))))
      } else {
        # case when there is only one calibration:
        root_age <- 1.1 * max(node_ages)
      }
    }
    node_ages <- c(root_age, node_ages)
    node_names <- c(root_node_name, node_names)
  }
  ############################################################################
  # bladj does not seem to freeze no more if calibrations are in conflict
  # it does ignore ages of descendant nodes that are older than ancestral nodes
  chronogram <- make_bladj_tree(tree = phy,
                                nodeages = node_ages,
                                nodenames = node_names)
  ############################################################################
  attributes(chronogram)$dating_method <- "bladj"
  # attributes(chronogram)$calibration_distribution <- ages_distribution
  attributes(chronogram)$used_calibrations <- stats::setNames(node_ages, node_names)
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
