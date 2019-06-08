# calibrations <- get_all_calibrations(phyloall)
# phy <- use_calibrations_bladj(phy, calibrations, type = "median")

# calibrations <- get_all_calibrations(cetaceae_phyloall)
# phy <- cetaceae_phyloall[[2]]
# plot(use_all_calibrations_bladj(phy,calibrations, use = "Mean"))

#' Use calibrations to date a topology with bladj.
#' @param phy A phylo object
#' @param calibrations A data frame of calibrations from get_all_calibrations function, or a subset of it.
#' @param type The type of age to use as calibration.
#' @return A phylo object
#' @export
use_calibrations_bladj <- function(phy, calibrations, type = "median"){
	type <- match.arg(tolower(type), c("mean", "min", "max", "median"))
	calibs <- match_all_calibrations(phy, calibrations)
    if(nrow(calibs$matched_calibrations) == 0){
      message("Nodes in calibrations (determined by taxon pairs) do not match any nodes in phy; phy cannot be dated")
      return(NA)
    }
	if("mean" %in% type){
	  node_ages <- sapply(calibs$phy$calibration_distribution, mean)
    }
    if("min" %in% type){
	  node_ages <- sapply(calibs$phy$calibration_distribution, min)
    }
    if("max" %in% type){
	  node_ages <- sapply(calibs$phy$calibration_distribution, max)
    }
	if("median" %in% type){
	  node_ages <- sapply(calibs$phy$calibration_distribution, stats::median)
    }
  # check that the root node has a calibration
  # otherwise, add one
  # bladj will run if calibrations are in conflict
  # it will not run if there is no calibration for the root
  node_names <- calibs$matched_calibrations$NodeNames
  if(!"n1" %in% calibs$matched_calibrations$NodeNames){
    node_ages <- c(max(node_ages) + mean(abs(diff(sort(node_ages)))), node_ages)
    node_names <- c("n1", node_names)
  }
	new_phy <- make_bladj_tree(tree = calibs$phy, nodeages = node_ages,
	  nodenames = node_names)
	new_phy$dating_method <- "bladj"
	new_phy$calibration_distribution <- calibs$phy$calibration_distribution
    new_phy$calibrations <- calibs$matched_calibrations
	return(new_phy)
}
#' function to check for conflicting calibrations
#' if calibrations are younger or older relative to descendants and ancestors, respectively
#' @details it removes them if needed, but bladj works as long as it has an age for the root
#' @param calibration_distribution is a list of node ages distributions, named with the node number from phy
# calibration_distribution <- calibs$phy$calibration_distribution
check_conflicting_calibrations <- function(phy, calibration_distribution){
  if(!inherits(phy, "phylo")){
    message("phy is not a phylo object")
    return(NA)
  }
  des <- lapply(as.numeric(names(calibration_distribution)), function(i) phy$edge[phy$edge[,1] == i,2])
  anc <- lapply(as.numeric(names(calibration_distribution)), function(i) phy$edge[phy$edge[,2] == i,1])
  foo1 <- function(node, des){
    des <- des[!is.na(des)]
    any(node < des)
  }
  foo2 <- function(node, anc){
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
