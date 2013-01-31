#' Test timing of name cleaning on datelife trees
#' 
#' @import plyr ggplot2 taxize stringr
#' @param listoftrees List of trees, each tree in phylo format, not multiPhylo
#' 		format. 
#' @param plot Plot results or not.
#' @examples \dontrun{
#' treefilenames <- dir("/Users/scottmac2/phyloorchard/pkg/data")
#' l_ply(treefilenames, function(x) load(paste("/Users/scottmac2/phyloorchard/pkg/data/",x,sep=""), .GlobalEnv))
#' trees <- sapply(treefilenames, function(x) str_replace(x, ".rda", ""), USE.NAMES=F)
#' binemonds <- BinindaEmondsEtAl2007[[1]]
#' makesmallsubtrees <- function(phylo){
#' 		drops <- c(2000,2500,3000,3500,4000,4100,4200,4300,4350,4400,4420,4450,4470,4500,4505)
#' 		llply(rev(drops), function(x) drop.tip(phylo, tip=1:x))
#' }
#' mytreelist <- makesmallsubtrees(phylo = binemonds)
#' library(ggplot2); library(plyr); library(lubridate)
#' out <- time_cleaning(listoftrees = mytreelist, plotresults=TRUE)
#' out[[1]] # data.frame
#' out[[2]] # plot
#' }
time_cleaning <- function(listoftrees, plotresults=FALSE){
	checknames_safe <- plyr::failwith(NULL, checknames)
	temp <- llply(listoftrees, function(tree){
			start <- Sys.time()
			result <- suppressMessages(checknames_safe(phylo=tree, source_="MSW3", splitby=100, byfilename=FALSE))
			stop_ <- Sys.time()
			elapsed <- as.duration(stop_-start)
			as.numeric(elapsed)
		}, .progress="text"
	)
	
	names(temp) <- sapply(listoftrees, Ntip, USE.NAMES=F)
	resultsout <- ldply(temp)
	names(resultsout) <- c("num_tips","elapsed")
	resultsout$num_tips <- as.numeric(resultsout$num_tips)
	
	if(plotresults){
		p <- ggplot(resultsout, aes(num_tips, elapsed)) + 
			theme_bw(base_size=18) +
			geom_point(size = 4) + 
			labs(y = "Elapsed time", x = "Number of tips")
		list(resultsout, p)
	} else
	{ resultsout }
}