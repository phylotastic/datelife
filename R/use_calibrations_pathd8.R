#' Use calibrations to date a tree with branch lengths with PATHd8.
#' @param phy A phylo object with branch lengths
#' @inheritParams use_calibrations_bladj
#' @param expand How much to expand by each step to get consistent calibrations. Should be between 0 and 1.
#' @param giveup How many expansions to try before giving up
#' @return A phylo object with branch length proportional to time.
#' @export
#' @details
#' This function will try to use the calibrations as fixed ages.
#' If that fails (often due to conflict between calibrations), it will expand the
#' range of the minimum age and maximum age and try again. And repeat.
#' If expand = 0, it uses the summarized calibrations.
#' In some cases, it returns edge lengths in relative time (with maximum tree depth = 1)
#' instead of absolute time, as given by calibrations. In this case, the function returns NA.
#' This is an issue from PATHd8.
# i=1
# phy <- tax_phyloall_bold[[3]][[1]]
# calibrations <- tax_othercalall[[3]][[1]]
# phy$edge.length
# calibrations$MaxAge
# calibs <- match_all_calibrations(tax_phyloall_bold[[1]][[i]], tax_othercalall[[1]][[i]])
# calibs$phy$edge.length
# used_calibrations <- calibs$matched_calibrations
# used_calibrations$MaxAge
# chronogram <- geiger::PATHd8.phylo(calibs$phy, used_calibrations)
# chronogram$edge.length
# plot(chronogram, main = i)
# ape::axisPhylo()
use_calibrations_pathd8 <- function(phy, calibrations, expand = 0.1, giveup = 100){
    phy <- input_process(phy, verbose = FALSE)
    if (!inherits(phy, "phylo")){
		    message("'phy' is not a phylo object")
        return(NA)
	  }
    if(is.null(phy$edge.length)){
        message("'phy' does not have branch lengths, consider using a dating method that does not require data, such as BLADJ or MrBayes.")
        return(NA)
    }
    # fix any negative branch lengths, otherwise pathd8 will silently not work:
    phy <- tree_fix_brlen(phy, fixing_criterion = "negative", fixing_method = 0, ultrametric = FALSE)
    # following line checks wether all calibrations are present in phy, and returns that in calibs$present_calibrations
    calibs <- match_all_calibrations(phy, calibrations)
    if(nrow(calibs$present_calibrations) < 1){
      message("\nDating analysis is not possible with this set of calibrations.")
      return(NA)
    }
    chronogram <- NA
    if(expand != 0) {
        attempts = 0
        used_calibrations <- calibs$present_calibrations[,c("MRCA", "MaxAge", "MinAge", "taxonA", "taxonB")]
        try(chronogram <- suppressWarnings(geiger::PATHd8.phylo(phy, used_calibrations)), silent = TRUE)
        # original_calibrations <- get_all_calibrations(phy2)
        while(!inherits(chronogram, "phylo") & attempts < giveup) {
            # print(rep(attempts, 10))
            used_calibrations <- calibs$present_calibrations
            used_calibrations$MaxAge <- as.numeric(used_calibrations$MaxAge) * ((1+expand)^attempts)
            used_calibrations$MinAge <- used_calibrations$MinAge * ((1-expand)^attempts)
            # We will have no fixed ages. Pathd8 just quietly gives up. So instead, we add a tiny branch with a zero calibration
            # between it and its sister.
            made.up.edgelength <- min(1e-9, .001*min(phy$edge.length))
            phy2 <- phytools::bind.tip(ape::reorder.phylo(phy), "tinytip", edge.length = made.up.edgelength, where = 1, position = made.up.edgelength) #bind tip has weird behavior for non-reordered trees
            used_calibrations[dim(used_calibrations)[1]+1,]<- c("fixed", 0, 0, phy$tip.label[1], "tinytip", "none")
            chronogram <- tryCatch(geiger::PATHd8.phylo(phy2, used_calibrations), error = function(e) NULL)
            # data that errors in geiger::PATHd8.phylo(phy2, used_calibrations) (and hence the tryCatch):
            # phy2 <- cetaceae_phyloall[2]# cetaceae_phyloall[3] and cetaceae_phyloall[4] also error
            # used_calibrations <- get_all_calibrations(cetaceae_phyloall[2])
            # Error in write.pathd8(phy, used_calibrations, base): Some calibrations not encountered in tree
            if(!is.null(chronogram)) {
                # chronogram$edge.length[which(chronogram$edge.length < 0)] <- 0
                #sometimes pathd8 returns tiny negative branch lengths. https://github.com/phylotastic/datelife/issues/11
                # fix this with tree_fix_brlen
                chronogram <- tree_fix_brlen(chronogram, fixing_criterion = "negative", fixing_method = 0, ultrametric = FALSE)
                chronogram <- ape::drop.tip(chronogram, "tinytip")
            }
            attempts <- attempts+1
        }
        if(attempts > 0 & inherits(chronogram, "phylo")) {
            message("Dates are even more approximate than usual: had to expand constraints to have them agree.")
        }
    } else { # alternative to expansion: summarize calibrations
        used_calibrations <- calibs$matched_calibrations
        try(chronogram <- geiger::PATHd8.phylo(calibs$phy, used_calibrations), silent = TRUE)
    }
    if(inherits(chronogram, "phylo")) {
        # sometimes pathd8 returns tiny negative branch lengths.
        # https://github.com/phylotastic/datelife/issues/11
        problem <- NULL
        if(is.null(chronogram$edge.length ) | all(is.na(chronogram$edge.length))){
            chronogram$edge.length <- NULL
            problem <- "PATHd8 returned a tree with no branch lengths."
        } else {
          if(all(chronogram$edge.length == 0)){
              chronogram$edge.length <- NULL
              problem <- "PATHd8 returned a tree with branch lengths equal to 0."
          }
          if(is.null(problem)){ # then fix negative branch lengths again
              # chronogram$edge.length[which(chronogram$edge.length<0)] <- 0
              chronogram <- tree_fix_brlen(tree = chronogram, fixing_criterion =
                  "negative", fixing_method = 0, ultrametric = TRUE)
          }
          if(round(max(ape::node.depth.edgelength(chronogram)), digits = 3) == 1){
              problem <- "Edge lengths seem to be relative to maximum age = 1 (and not absolute to time given by calibrations)."
          }
        }
        if(!is.null(problem)){
            message(paste(problem, "\nThis is an issue from PATHd8; returning tree with a $problem."))
        }
        chronogram$problem <- problem
        chronogram$dating_method <- "pathd8"
    	  chronogram$calibration_distribution <- calibs$phy$calibration_distribution
        chronogram$used_calibrations <- used_calibrations
        chronogram$present_calibrations <- calibs$present_calibrations
    } else {
      warning("Dating analysis with PATHd8 failed with this tree and set of calibrations.")
      return(NA)
    }
	return(chronogram)
}
