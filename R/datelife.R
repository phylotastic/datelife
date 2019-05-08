# calibrations <- get_all_calibrations(phyloall)
# phy <- use_calibrations_bladj(phy, calibrations, type = "median")

# calibrations <- get_all_calibrations(cetaceae_phyloall)
# phy <- cetaceae_phyloall[[2]]
# plot(use_all_calibrations_bladj(phy,calibrations, use = "Mean"))
#' Use calibrations to date a topology with bladj.
#' @param phy A phylo object
#' @param calibrations A data frame of calibrations from get_all_calibrations function
#' @param type A character vector
#' @return A phylo object
#' @export
use_calibrations_bladj <- function(phy, calibrations, type = "median"){
	type <- match.arg(tolower(type), c("mean", "min", "max", "median"))
	calibs <- map_all_calibrations(phy, calibrations)
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
	new_phy <- make_bladj_tree(tree = calibs$phy, nodeages = node_ages,
	  nodenames = as.character(calibs$calibrations$MRCA))
	new_phy$dating_method <- "bladj"
	new_phy$calibration_distribution <- stats::setNames(calibs$phy$calibration_distribution)
	return(new_phy)
}
