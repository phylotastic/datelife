#' Use calibrations to date a topology with the BLADJ algorithm.
#'
#' @description The function prepares the input for BLADJ and calls [make_bladj_tree()]
#' @param calibrations A `data.frame` of secondary calibrations for any pair of taxon
#' names in `phy`, usually obtained with [get_all_calibrations()].
#' @param type The type of age to use as calibration. Options are "median", "mean", "min", or "max".
#' @param root_age Not implemented yet. Numeric specifying the age of the root. If there are no calibrations for it. If NULL or not numeric, the maximum calibration plus a unit of the mean differences will be used as root calibration. If there is only one internal calibration, the root age will be set to 10% more than the age of the calibration.
#' @return A `phylo` object with branch lengths proportional to time.
#' @details
#' The BLADJ algorithm is part of the Phylocom software, presented in Webb et al.
#' (2008) \doi{10.1093/bioinformatics/btn358}.
#' @references
#' Webb, C. O., Ackerly, D. D., & Kembel, S. W. (2008). "Phylocom: software for
#' the analysis of phylogenetic community structure and trait evolution".
#' Bioinformatics, 24(18), \doi{10.1093/bioinformatics/btn358}.

#' @export
use_calibrations_bladj.matchedCalibrations <- function(calibrations,
                                                       type = "mean",
                                                       root_age = NULL) {
  #
  message("... Using secondary calibrations with BLADJ.")
  is_good_calibrations <- inherits(calibrations, "congruifiedCalibrations") | !inherits(calibrations, "matchedCalibrations")
  if (!is_good_calibrations) {
    stop("'calibrations' must be a `matchedCalibrations` or `congruifiedCalibrations` object.\n")
  }
  type <- match.arg(tolower(type), c("mean", "min", "max", "median"))
  if (!all(calibrations$MinAge == calibrations$MaxAge)) {
    stop("Min and max calibrations provided are different, please choose one.")
  }
  if (nrow(calibrations) < 1) {
    stop("Nodes in 'calibrations' (determined by taxon pairs) do not match any nodes
						in 'phy'.\n\t Dating analysis is not possible with this set of calibrations.")
  }
  if ("mean" %in% type) {
    #node_ages <- sapply(calibrations$phy$calibration_distribution, mean)
    node_ages <- sapply(unique(calibrations$mrca_node_name), function(x) {
      node <- calibrations$mrca_node_name == x
      mean(calibrations$MinAge[node])})
  }
  if ("min" %in% type) {
    node_ages <- sapply(attributes(calibrations)$phy$calibration_distribution, min)
  }
  if ("max" %in% type) {
    node_ages <- sapply(attributes(calibrations)$phy$calibration_distribution, max)
  }
  if ("median" %in% type) {
    node_ages <- sapply(attributes(calibrations)$phy$calibration_distribution, stats::median)
  }
  # check that the root node has a calibration
  # otherwise, add one
  # bladj does not run if there is no calibration for the root
  # bladj will run if calibrations are in conflict
  node_names <- names(node_ages)
  # length(node_names)
  root_node_name <- paste0("n", ape::Ntip(attributes(calibrations)$phy) + 1)
  if (!root_node_name %in% node_names) {
    if (length(node_ages) > 1) {
      root_age <- max(node_ages) + mean(abs(diff(sort(node_ages))))
    } else {
      # if there is only one calibration the line above will give NaN
      root_age <- 1.1 * max(node_ages)
    }
    # add root node name and age:
    node_ages <- c(root_age, node_ages)
    node_names <- c(root_node_name, node_names)
  }
  chronogram <- make_bladj_tree(
    tree = attributes(calibrations)$phy, nodeages = node_ages,
    nodenames = node_names
  )
  # TODO: something about these extra list elements, set up as attributes??
  chronogram$dating_method <- "bladj"
  # chronogram$calibration_distribution <- calibrations$phy$calibration_distribution
  chronogram$calibrations <- calibrations
  chronogram$used_calibrations <- stats::setNames(node_ages, node_names)
  return(chronogram)
}
