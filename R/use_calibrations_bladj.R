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
    if(nrow(calibs$calibrations) == 0){
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
  node_names <- calibs$calibrations$NodeNames
  if(!"n1" %in% calibs$calibrations$NodeNames){
    node_ages <- c(max(node_ages) + mean(abs(diff(sort(node_ages)))), node_ages)
    node_names <- c("n1", node_names)
  }
	new_phy <- make_bladj_tree(tree = calibs$phy, nodeages = node_ages,
	  nodenames = node_names)
	new_phy$dating_method <- "bladj"
	new_phy$calibration_distribution <- calibs$phy$calibration_distribution
    new_phy$calibrations <- calibs$calibrations
	return(new_phy)
}
# #' functions to check if calibrations are younger or older relative to descendants and ancestors, respectively
# #' @details it removes them if needed, but bladj works as long as it has an age for the root
# check_nodeages_bladj <- function(phy, node_ages){
#   if(!inherits(phy, "phylo")){
#     message("phy is not a phylo object")
#     return(NA)
#   }
#   des <- lapply(as.numeric(names(node_ages)), function(i) phy$edge[phy$edge[,1] == i,2])
#   anc <- lapply(as.numeric(names(node_ages)), function(i) phy$edge[phy$edge[,2] == i,1])
#   foo1 <- function(node, des){
#     des <- des[!is.na(des)]
#     any(node < des)
#   }
#   foo2 <- function(node, anc){
#     anc <- anc[!is.na(anc)]
#     any(node > anc)
#   }
#   node_ages_des <- lapply(des, function(x) node_ages[as.character(x)])
#   node_ages_anc <- lapply(anc, function(x) node_ages[as.character(x)])
#   node_younger <- mapply(foo1, node_ages, node_ages_des) # which node_ages are younger than their descendants
#   node_older <- mapply(foo2, node_ages, node_ages_anc) # which node_ages are older than their ancestor
#   list(des, anc, node_ages_des, node_ages_anc, node_younger, node_older)
# }
# fix_nodeages_for_bladj <- function(phy, node_ages, remove = TRUE, remove_which = "younger", expand = FALSE){
#   remove <- match.arg(tolower(remove), c("older", "younger"))
#   node_ages_original <- node_ages
#   cc <- check_nodeages_bladj(phy, node_ages)
#   if("older" %in% remove){
#     node_ages <- node_ages[!cc$node_older]
#   }
#   if("younger" %in% remove){
#     node_ages <- node_ages[!cc$node_younger]
#   }
#   # the root needs to be fixed in bladj always,
#   # so if it is the root that is younger than any of its descendants, we need to fix that
#   # we will just add the difference:
#   if(cc$node_younger[1]){
#     node_ages[1] <- max(cc$node_ages_des[[1]], na.rm = TRUE) +
#       abs(max(cc$node_ages_des[[1]], na.rm = TRUE) - node_ages[1])
#     cc$node_younger[1] <- FALSE
#   }
#   node_ages
# }
