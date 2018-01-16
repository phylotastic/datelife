#skeleton of code to make open access info on dating life
# library(ape)
# library(abind)
# library(phangorn)
# library(compare)
# library(geiger) #note this uses the next version of geiger, which has congruifier code. The relevant code is in auteur-congruify.R
# library(datelife2) #Klaus Schliep's code for pruning efficiently. Eventually, this will be moved into phangorn
# library(parallel)
# library(doMC)
# library(stringr)
# library(ggplot2)
# library(taxize)
# library(plyr)
# library(rotl)
# source("/Library/WebServer/Sites/datelife.org/datelife/R/cleaning.r")

#' Core function to input a vector of species, newick string, or phylo object to get a chronogram or dates back.
#' @aliases datelife
#' @param input Target taxa names in the form of a vector of characters, a newick character string, or a phylo object.
#' @param output.format The desired output format for target chronograms (chronograms of target taxa). See details.
#' @param showSummary A character vector specifying type of summary information to be printed: "citations" for the references of chronograms from cache where target taxa are found, "taxa" for a summary of the number of chronograms where each target taxon is found, or "none" if nothing should be printed. Default to display both c("citations", "taxa").
#' @param missing.taxa A character vector specifying if data on target taxa missing in source chronograms should be added to the output as a "summary" or as a presence/absence "matrix". Default to "none", no information on missing.taxa added to the output.
#' @param partial If TRUE, use source chronograms even if they only match some of the desired taxa
#' @param usetnrs If TRUE, use OpenTree's services to resolve names. This can dramatically improve the chance of matches, but also take much longer.
#' @param approximatematch If TRUE, use a slower TNRS to correct mispellings, increasing the chance of matches (including false matches).
#' @param update_cache default to FALSE
#' @param cache The cached set of chronograms and other info from data(opentree_chronograms).
#' @param method The method used for congruification. PATHd8 only right now, r8s and treePL later.
#' @param bold Logical. If TRUE, use Barcode of Life Data Systems (BOLD)  and Open Tree of Life (OToL) backbone to estimate branch lengths of target taxa using GetBoldOToLTree function.
#' @param marker A character vector with the name of the gene from Barcode of Life Data Systems (BOLD) to be used for branch length estimation.
#' @param sppfromtaxon boolean vector, default to FALSE. If TRUE, will get all species names from taxon names given in input. Must have same length as input. If input is a newick string , with some clades it will be converted to phylo object phy, and the order of sppfromtaxon will match phy$tip.label.
#' @param verbose Boolean. If TRUE, it gives printed updates to the user.
#' @inheritDotParams GetBoldOToLTree otol_version chronogram doML
#' @export
#' @details
#' Available output formats are:
#'
#' citations: A character vector of references where chronograms with some or all of the target taxa are published (source chronograms).
#'
#' mrca: A named numeric vector of most recent common ancestor (mrca) ages of target taxa defined in input, obtained from the source chronograms. Names of mrca vector are equal to citations.
#'
#' newick.all: A named character vector of newick strings corresponding to target chronograms derived from source chronograms. Names of newick.all vector are equal to citations.
#'
#' newick.sdm: Only if multiple source chronograms are available. A character vector with a single newick string corresponding to a target chronogram obtained with SDM supertree method (Criscuolo et al. 2006).
#'
#' newick.median: Only if multiple source chronograms are available. A character vector with a single newick string corresponding to a target chronogram from the median of all source chronograms.
#'
#' phylo.sdm: Only if multiple source chronograms are available. A phylo object with a single target chronogram obtained with SDM supertree method (Criscuolo et al. 2006).
#'
#' phylo.median: Only if multiple source chronograms are available. A phylo object with a single target chronogram obtained from source chronograms with median method.
#'
#' phylo.all: A named list of phylo objects corresponding to each target chronogram obtained from available source chronograms. Names of phylo.all list correspond to citations.
#'
#' html: A character vector with an html string that can be saved and then opened in any web browser. It contains a 4 column table with data on target taxa: mrca, number of taxa, citations of source chronogram and newick target chronogram.
#'
#' data.frame A data.frame with data on target taxa: mrca, number of taxa, citations of source chronograms and newick string.
#' @examples
#' # obtain median ages from a set of source chronograms in newick format:
#' ages <- EstimateDates(c("Rhea americana", "Pterocnemia pennata", "Struthio camelus",
#' 		"Mus musculus"), output.format="newick.median")
#' # save the tree in newick format
#' write(ages, file="some.bird.ages.txt")
#'
#' # obtain median ages from a set of source chronograms in phylo format
#' # will produce same tree as above but in r phylo format:
#' ages.again <- EstimateDates(c("Rhea americana", "Pterocnemia pennata", "Struthio camelus",
#' 		"Mus musculus"), output.format="phylo.median")
#' plot(ages.again)
#' library(ape)
#' ape::axisPhylo()
#' mtext("Time (million years ago)", side=1, line=2, at = (max(get("last_plot.phylo",
#' 		envir = .PlotPhyloEnv)$xx) * 0.5))
#' write.tree(ages.again, file="some.bird.tree.again.txt") # saves phylo object in newick format
#'
#' # obtain mrca ages and target chronograms from all source chronograms
#' # generate an html  output readable in any web browser:
#' ages.html <- EstimateDates(c("Rhea americana", "Pterocnemia pennata", "Struthio camelus",
#' 		"Mus musculus"), output.format="html")
#' write(ages.html, file="some.bird.trees.html")
#' system("open some.bird.trees.html")

EstimateDates <- function(input=c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"),
		output.format = "phylo.all", partial = TRUE, usetnrs = FALSE, approximatematch = TRUE, update_cache = FALSE, cache = get("opentree_chronograms"), method="PATHd8", bold=FALSE, showSummary= c("citations", "taxa"), missing.taxa = c("none", "summary", "matrix"),  marker = "COI", sppfromtaxon=FALSE, verbose=FALSE, ...) {
			#... only defines arguments to be passed to GetBoldOToLTree for now
			# find a way not to repeat partial and cache arguments, which are used in both GetFilteredResults and SummarizeResults
			if(update_cache){
				cache <- UpdateCache(save = TRUE, verbose=verbose)
			}
			input.here <- ProcessInput(input=input, usetnrs=usetnrs, approximatematch=approximatematch, sppfromtaxon=sppfromtaxon, verbose=verbose)
			filtered.results.here <- GetFilteredResults(input = input.here, partial = partial, usetnrs = usetnrs, approximatematch = approximatematch, update_cache = FALSE, cache = cache, method = method, bold = bold, marker = marker, verbose=verbose, ...)
			return(SummarizeResults(input = input.here$cleaned.names, filtered.results = filtered.results.here, output.format = output.format, partial = partial, update_cache = FALSE, cache = cache, showSummary = showSummary, missing.taxa = missing.taxa, verbose=verbose))
}

#' Go from a vector of species, newick string, or phylo object to a list of patristic matrices
#' @inheritParams EstimateDates
#' @inheritParams ProcessInput
#' @inheritParams GetBoldOToLTree
#' @inheritDotParams GetBoldOToLTree
#' @return List of patristic matrices
#' @export
GetFilteredResults <- function(input=c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"), partial=TRUE, usetnrs=FALSE, approximatematch=TRUE, update_cache = FALSE, cache=get("opentree_chronograms"), method="PATHd8", bold=FALSE, marker = "COI", sppfromtaxon=FALSE, verbose=FALSE, ...) {
	if(update_cache){
		cache <- UpdateCache(save = TRUE, verbose=verbose)
	}
	input <- CheckInput(input=input, usetnrs = usetnrs, approximatematch = approximatematch, sppfromtaxon = sppfromtaxon, verbose=verbose)
	tree <- input$phy
	cleaned.names <- input$cleaned.names
	CheckCleanedNames(cleaned.names=cleaned.names, sppfromtaxon=sppfromtaxon, verbose=verbose)
  results.list <- lapply(cache$trees, GetSubsetArrayDispatch, taxa=cleaned.names, phy=tree, method=method)
  filtered.results <- ProcessResultsList(results.list, cleaned.names, partial)
	CheckFilteredResults(filtered.results, usetnrs)
	if(bold){
		 bold.OToLTree <- GetBoldOToLTree(input = input, usetnrs = FALSE, approximatematch = FALSE, marker = marker, verbose=verbose,  ...)
		 bold.data <- GetSubsetArrayBothFromPhylo(reference.tree.in = bold.OToLTree, taxa.in = cleaned.names, phy.in = NULL, phy4.in = NULL, method.in = method)
		 bold.data.processed <- ProcessResultsList(results.list=list(bold.data), taxa=cleaned.names, partial)
	 	 names(bold.data.processed) <-  paste("BoldOToL tree (using ", marker, " as marker)", sep="")
	   filtered.results <- c(filtered.results, bold.data.processed)
	}
#	cat("\n")
	return(filtered.results)
}

#' checks if input has already been processed, otherwise it uses ProcessInput
#' @inheritParams EstimateDates
#' @inheritDotParams ProcessInput
CheckInput <- function(input, ...){
	badformat <- TRUE
	if(is.list(input) & "phy" %in% names(input) & "cleaned.names" %in% names(input)) badformat <- FALSE
	if(badformat){
		input <- ProcessInput(input = input, ...)
	}
	return(input)
}
#' checks that we have at least two taxon names to perform a search
#' @inheritParams EstimateDates
#' @param cleaned.names A character vector; usually an output from ProcessInput function
CheckCleanedNames <- function(cleaned.names, sppfromtaxon, verbose=FALSE){
	if(length(cleaned.names)==1){
		if(verbose) cat("Cannot perform a search of divergence times with just one taxon.", "\n")
		if(sppfromtaxon) {
			if(verbose) cat("Clade contains only one lineage.", "\n")
		} else {
			if(verbose) cat("Performing a clade search? set sppfromtaxon=TRUE.", "\n")
		}
		stop("input is length 1")
	}
}
#' checks if we obtained an empty search with the set of input taxon names
#' @inheritParams EstimateDates
#' @param filtered.results An object output from GetFilteredResults function
CheckFilteredResults <- function(filtered.results, usetnrs, verbose=FALSE){
	if(length(filtered.results) < 1) {
		warning("Output is empty.", call. = FALSE)
		if(verbose) cat("Input species were not found in any chronograms available in cache.", "\n")
		if(!usetnrs & verbose) cat("Setting usetnrs = TRUE might change this, but it is time consuming.", "\n")
	}
}
#' Take input phylo object or character string and figure out if it's correct newick format or a list of species
#' @inheritParams EstimateDates
#' @inheritParams ProcessInput
#' @return A phylo object or NA if no tree
#' @export
ProcessPhy <- function(input, verbose=FALSE){
  	if(class(input) == "multiPhylo") stop("Only one phylogeny can be processed at a time.")
	if(class(input) == "phylo") {
		input <- ape::write.tree(input)
	}
 	input <- gsub("\\+"," ",input)
  	input <- stringr::str_trim(input, side = "both")
  	phy.new.in <- NA
   	# if(length(input) == 1) {
    	# if(showSummary)cat("\t", "Input is length 1.", "\n")
	  	if(any(grepl("\\(.*\\).*;", input))) { #our test for newick
	  		if(length(input)>1) stop("Only one phylogeny can be processed at a time.")
	    	phy.new.in <- ape::collapse.singles(phytools::read.newick(text=gsub(" ", "_", input)))
	    	if(verbose) {cat("\t", "Input is a phylogeny and it is correcly formatted.", "\n")}
	  	} else {
	  		if(verbose) {cat("Input is not a phylogeny.")} #not a warning nor stop, 'cause it is not a requirement for input to be a phylogeny at this step
	  	}
  	# }
	return(phy.new.in)
}



#' Cleans taxon names from input character vector, phylo object or newick character string. Process the two latter with ProcessPhy first.
#' @inheritParams EstimateDates
#' @inheritDotParams rphylotastic::GetSpeciesFromTaxon filters
#' @return A list with the phy (or NA, if no tree) and cleaned vector of taxa
#' @export
ProcessInput <- function(input=c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"), usetnrs=FALSE, approximatematch=TRUE, sppfromtaxon=FALSE, verbose=FALSE, ...) {
	if(verbose) cat("Processing input...", "\n")
	phy.new <- ProcessPhy(input = input, verbose=verbose)
	# cleaned.names <- ""
	if(!is.na(phy.new[1])) {
	    # if(usetnrs) {
	      	# phy.new$tip.label <- gsub("_", " ", cleaned.names)
	    # }
	  	# cleaned.names <- gsub("_", " ", phy.new$tip.label)
	  	input <- phy.new$tip.label
	}
	if(length(input)==1) {
		input <- strsplit(input, ',')[[1]]
		if(!sppfromtaxon[1]) {
				if(verbose) {
					cat("Datelife needs at least two input taxon names to perform a search.", "\n")
					cat("Setting sppfromtaxon = TRUE gets all species from a clade and accepts only one taxon name as input.", "\n")
				}
				stop("Input is length 1 and not in a good newick format.")
		}
	}
	cleaned.input <- stringr::str_trim(input, side = "both")
    if (usetnrs) {
		cleaned.input <- rotl::tnrs_match_names(cleaned.input)$unique_name # process even if it's a "higher" taxon name
    }
	cleaned.names <- gsub("_", " ", cleaned.input)
    if(any(sppfromtaxon)){
    	if(length(sppfromtaxon)==1) sppfromtaxon <- rep(sppfromtaxon,length(cleaned.input))
    	if(length(cleaned.input)!=length(sppfromtaxon)){
    		if(verbose) cat("Specify all taxa in input to get species names from.", "\n")
    		stop("input and sppfromtaxon arguments must have same length.")
    	}
    	species.names <- vector()
    	index <- 1
	    for (i in sppfromtaxon){
	    	if (i) {
	    		spp <- rphylotastic::GetSpeciesFromTaxon(taxon = cleaned.names[index], ...)
	    		if(length(spp)==0) {
	    			if(verbose) cat("\t", " No species names found for taxon ", cleaned.names[index], ".", "\n", sep="")
	    			if (!usetnrs & verbose) cat("\t", "Setting usetnrs = TRUE might change this, but it can be slow.", "\n")
	    			warning(paste("No species names available for input taxon '", cleaned.names[index], "'", sep=""))
	    		}
	    		species.names <- c(species.names, spp)
	    	} else {
	    		species.names <- c(species.names, cleaned.names[index])
	    	}
	    	index <- index + 1
	    }
		cleaned.names <- gsub("_", " ", species.names)
	}
	cleaned.names <- unique(cleaned.names)
    if(verbose) cat("OK.", "\n")
  	cleaned.names.print <- paste(cleaned.names, collapse = " | ")
  	if(verbose) cat("Working with the following taxa:", "\n", "\t", cleaned.names.print, "\n")
   	return(list(cleaned.names=cleaned.names, phy=phy.new))
}

#' Takes a tree and fixes negative branch lengths in several ways
#' @param phy A tree either as a newick character string or as a phylo object
#' @param method A character vector specifying the method to fix negative branch lengths: "zero", "bladj" or "mrbayes"
#' @return A phylo object with no negative branch lengths
#' @export
FixNegBrLen <- function(phy=NULL, method = "zero"){
	phy <- ProcessPhy(input=phy, verbose=FALSE)
	if (!inherits(phy, "phylo"))
		stop("phy must be a newick character string or in phylo format")
	if(is.null(phy$edge.length))
		stop("phy must have branch lengths")
	if(!ape::is.ultrametric(phy))
		stop("branch lengths must be relative to time")
	method <- match.arg(method, c("zero", "bladj", "mrbayes"))

	index <- which(phy$edge.length<0)  # identifies edge numbers with negative edge lengths value

	if(method=="zero") {# chunk for neg br len to zero
		for (i in index){
			# snode <- pos.phy$edge[i,1]
			# pool  <- pos.phy$edge[seq(nrow(pos.phy$edge))[-i], 1]
			# sisedge <- which(pool==snode) # determines position of sister edge
			# pos.phy$edge.length[sisedge] <- pos.phy$edge.length[sisedge] - pos.phy$edge.length[i]
			# adds neg branch length to sister branch, should add error to both sides???? or only to the daughter branches??
			cnode <- phy$edge[i,2]
			dauedge <- which(phy$edge[,1]==cnode)
			phy$edge.length[dauedge] <- phy$edge.length[dauedge] + phy$edge.length[i]
			phy$edge.length[i] <- 0
			fixed.phy <- phy
		}
	}

	if(any(method==c("bladj", "mrbayes"))) { #chunk for bladj and mrbayes
	# if(method=="bladj") { #chunk for bladj
		phy <- add_node_labels(phy=phy)  # all nodes need to be named
		cnode <- phy$edge[index,2]  # we assume that the negative edge length is the one that needs to be changed (but it could be the sister edge that should be shorter)
		tofix <- cnode-length(phy$tip.label)  # so, we take the crown node number of the negative branch lengths
		if(method=="bladj")
			fixed.phy <- GetBladjTree(nodenames = phy$node.label[-tofix], nodeages = get_node_data(phy=phy, node_data="node_age")$node_age[-tofix], phy = phy, phyformat = "phylo")
		if(method=="mrbayes") {
			mrbayes.file <- paste0("datelife_query", "_negBrLen_fixed.nexus") # make this an argument
			ncalibration <- get_node_data(phy, node_data=c("descendant_tips_label", "node_age"))
			ncalibration <- lapply(ncalibration, "[", seq(phy$Nnode)[-tofix])
			phy1 <- AddOutgroup(phy=phy, outgroup="fake_outgroup", processphy=FALSE)
			fixed.phy <- GetMrBayesTree(phy=phy1, ncalibration=ncalibration, file=mrbayes.file)
			fixed.phy <- ape::drop.tip(fixed.phy, "fake_outgroup")
		}
	}
	return(fixed.phy)
}

#' Gets ages, node numbers, node names and descendant tips number and label of all nodes from a dated tree
#' @inheritParams FixNegBrLen
#' @param node_data A character vector containing one or all from: "node_number", "node_label", "node_age", "descendant_tips_number", "descendant_tips_label"
#' @return A list
get_node_data <- function(phy=NULL, node_data=c("node_number", "node_label", "node_age", "descendant_tips_number", "descendant_tips_label")){
	res <- vector(mode="list")
	phy <- ProcessPhy(input=phy)
	phy <- CheckPhylo(phylo = phy, dated=TRUE)
	nn <- get_node_numbers(phy)
	if("node_number" %in% node_data){
		res <- c(res, list(node_number = nn))
	}
	if("node_label" %in% node_data){
		res <- c(res, list(node_label = phy$node.label))
	}
	if("node_age" %in% node_data){
		res <- c(res, list(node_age = ape::branching.times(phy)))
	}
	if(any(c("descendant_tips_number", "descendant_tips_label") %in% node_data)){
		dt_num <- lapply(nn, function(x) GetTips(tree=phy, node=x))# tip numbers stemming from node
		names(dt_num) <- nn
		dt_lab <- lapply(dt_num, function(x) phy$tip.label[x]) # tip labels corresponding to those tips numbers
		names(dt_lab) <- nn
		if("descendant_tips_number" %in% node_data) {
			res <- c(res, list(descendant_tips_number=dt_num))
		}
		if("descendant_tips_label" %in% node_data) {
			res <- c(res, list(descendant_tips_label=dt_lab))
		}
	}
	return(res)
}

#' Adds labels to nodes with no asigned label
#' @inheritParams FixNegBrLen
#' @param node_prefix Character vector. If length 1, it will be used to name all nodes with no labels, followed by a number which can be the node_number or consecutive, as specified in node_number
#' @param node_index Character vector choosing one of "consecutive" or "node_number" as node label index. It will use consecutive numbers from 1 to total node number in the first case and node numbers in the second case.
#' @return A phylo object
add_node_labels <- function(phy=NULL, node_prefix="n", node_index="node_number"){
	phy <- ProcessPhy(input = phy)
	phy <- CheckPhylo(phylo=phy, dated=FALSE)
	node_index <- match.arg(arg=node_index, choices=c("consecutive","node_number"), several.ok = FALSE)
	if("node_number" %in% node_index){
		node_number <- get_node_numbers(phylo = phy)
	}
	if("consecutive" %in% node_index){
		node_number <- seq(phy$Nnode)
	}
	if(is.null(phy$node.label)){
		phy$node.label <- paste0(node_prefix, node_index)
	} else {
		en <- which(phy$node.label=="")
		phy$node.label[en] <- paste0(node_prefix, en)
	}
	return(phy)
}

#' Gets node numbers from any phylogeny
#' @inheritParams CheckPhylo
#' @return A numeric vector with node numbers
get_node_numbers <- function(phylo){
	node_numbers <- (length(phylo$tip.label)+1):(length(phylo$tip.label)+phylo$Nnode)
	return(node_numbers)
}

#' Function to add an outgroup to any phylogeny, in phylo or newick format
#' @inheritParams GetMrBayesTree
#' @param outgroup A character vector with the name of the outgroup. If it has length>1, only first element will be used.
#' @param processphy Boolean. If true, phy will be processed with ProcessPhy function.
#' @return A phylo object.
#' @export
AddOutgroup <- function(phy=NULL, outgroup="outgroup", processphy=TRUE){
    if (processphy) {
			phy <- ProcessPhy(input=phy, verbose=FALSE)
			phy <- CheckPhylo(phylo=phy, dated=FALSE)
    }
    # if(phy$edge.length) # comment cause it doesn't matter if phy has branch lengths
    outgroup_edge <- c()
    ingroup_edge <- c()
    if(!is.null(phy$edge.length)){
    	mbt <- max(ape::branching.times(phy))
    	outgroup_edge <- mbt + mbt*0.10
    	ingroup_edge <- paste0(":", outgroup_edge-mbt)
    	outgroup_edge <- paste0(":", outgroup_edge)
    }
	phy <- gsub(";", "", ape::write.tree(phy))
	phy <- paste("(", phy, ingroup_edge, ",", outgroup[1], outgroup_edge, ");", sep="")
	phy <-  phytools::read.newick(text=phy)
	return(phy)
}

#' Checks if phy is a phylo class object and a chronogram
#' @param phylo A phylo object
#' @param dated Boolean. If TRUE it checks if phylo object has branch lengths and is ultrametric.
#' @return A phylo object
CheckPhylo <- function(phylo, dated=FALSE){
	if (!inherits(phylo, "phylo")){
		stop("phy is not a phylo object")
	}
	if(dated){
		if(is.null(phylo$edge.length))
			stop("phylo object must have branch lengths")
		if(!ape::is.ultrametric(phylo)){
			warning("branch lengths in phylo object should be proportional to time")
			stop("phylo object must be ultrametric")  # is this true?  # Think how to incorporate trees with extinct taxa
		}
	}
	return(phylo)
}

#' To get tip numbers descending from any given node of a tree
#' @inheritParams phytools::getDescendants
#' @return A numeric vector with tip numbers descending from a node
#' @export
GetTips <- function(tree=NULL, node=NULL, curr=NULL){
	des <- phytools::getDescendants(tree=tree, node=node, curr=NULL)
	tips <- des[which(des<=length(tree$tip.label))]
	return(tips)
}

#' Takes a tree and uses bladj to estimate node ages and branch lengths given a set of fixed node ages and respective node names
#' @param nodenames A character vector with node names from tree with fixed ages
#' @param nodeages A numeric vector with known or fixed node ages from tree
#' @param phy A tree either as a newick character string or phylo format
#' @param phyformat A character vector specifying tree output format, either "newick" (default) or "phylo"
#' @return A newick or phylo tree with non negative branch lengths
#' @export
GetBladjTree <- function(nodenames=NULL, nodeages=NULL, phy=NULL, phyformat="newick"){
	phyformat <- match.arg(phyformat, choices = c("newick", "phylo"))
	if(is.null(phy$node.label)) stop("phy must have node labels")
	if(!is.null(phy$edge.length)) phy$edge.length <- NULL
	m <- match(nodenames, phy$node.label)
	if(any(is.na(m))) stop("all nodenames must be in phy$node.label") # add a printed line saying which nodenames are not in phy$node.label
	if(length(nodenames)!=length(nodeages)) stop("nodenames and nodeages must have the same length")
	if(!is.character(nodenames)) stop("nodenames must be a character vector")
	if(!is.numeric(nodeages)) stop("nodeages must be a numeric vector")
	ages_df <- data.frame(
		a=nodenames,
		b=nodeages
	)
	new.phy <- phylocomr::ph_bladj(ages = ages_df, phylo = phy)
	attributes(new.phy) <- NULL
	if(phyformat == "phylo") new.phy <- ape::read.tree(text = new.phy)
	# plot(new.phy)
	return(new.phy)
}

#' Takes a constraint tree and uses mrBayes to get node ages and branch lengths given a set of node calibrations
#' @param phy The constraint tree: a phylo object or a newick character string, with or without branch lengths.
#' @param ncalibration The node calibrations: a phylo object with branch lengths proportional to time; in this case all nodes from ncalibration will be used as calibration points. Alternatively, a list with two elements: the first is a character vector with node names from phy to calibrate; the second is a numeric vector with the corresponding ages to use as calibrations.
#' @inheritParams CheckMissingTaxa
#' @param file A character vector specifying the name of mrBayes run file and outputs (can specify directory too).
#' @return A phylo tree with branch lengths proportional to time. It will save all mrBayes outputs in the working directory.
#' @export
GetMrBayesTree <- function(phy=NULL, ncalibration=NULL, missingTaxa = NULL, file="mrbayes_run.nexus"){
	MakeMrBayesRunFile(phy= phy, ncalibration= ncalibration, missingTaxa = missingTaxa, file=file)
	new.tree <- MrBayesRun(file=file)
	return(new.tree)
}

#' Makes a mrBayes run block file with a constraint topology and a set of node calibrations
#' @inheritParams GetMrBayesTree
#' @return A MrBayes block run file in nexus format.
#' @export
MakeMrBayesRunFile <- function(phy= NULL, ncalibration= NULL, missingTaxa = NULL, file="mrbayes_run.nexus"){
  phy <- ProcessPhy(input=phy, verbose=FALSE) # add outgroup=TRUE argument
	if (!inherits(phy, "phylo")) {
		stop("phy must be a newick character string or in phylo format")
    }
  phy <- ConvertSpacesToUnderscores(phy)
	taxa <- phy$tip.label
	constraints <- paleotree::createMrBayesConstraints(tree= phy, partial=FALSE) # this works perfectly
	calibrations <- GetMrBayesNodeCalibrations(phy=phy, ncalibration=ncalibration, ncalibrationType = "fixed")
	og <- IdentifySingletonOutgroup(phy) #if(outgroup)
	if(!is.na(og)) ogroup <- paste0("outgroup ", og, ";")
	missingTaxa <- CheckMissingTaxa(missingTaxa)
	if(!is.null(missingTaxa)){
		if(is.vector(missingTaxa)){
			taxa <- c(taxa, missingTaxa)
		}
	}
	bayes_data <- c(paste("   Begin DATA; \nDimensions ntax=", length(taxa), "nchar=1;"),
	"Format datatype=DNA gap=- missing=?;",
	"Matrix\n",
	paste(taxa, "?"),
	";")

    bayes_set <- c("   Begin MRBAYES;",
    	"unlink shape=(all) tratio=(all) statefreq=(all) revmat=(all) pinvar=(all);\n",
    	constraints[-length(constraints)],
    	ogroup, "",
    	constraints[length(constraints)], "",
    	"prset nodeagepr=calibrated;", "",
    	calibrations, "\n",
    	"   set usebeagle=no Beaglesse=no;", "",
    	paste("prset ", c("brlenspr=clock:birthdeath", "Extinctionpr = Fixed(0)",
    	"Speciationpr=exponential(1)", "clockvarpr=ibr", "ibrvarpr=exponential(10)"), ";", sep=""),
    	"mcmcp nruns=1 nchains=1 ngen=50000000 samplefreq=1000;",
    	"mcmc;", "",
    	paste0("sumt filename=", file, " burnin=5000000 contype=halfcompat;\n"),
    	"end;"
    	)

	all <- c(bayes_data, "\n", bayes_set)
	write(all, file)
}

#' Makes up dates on an already dated tree for missing taxa.
#' @param  missingTaxa either a tree (as phylo object or as a newick character string, with or without branch lengths) with all the taxa you want at the end, or a data frame assigning missing taxa to nodes in dated tree.
#' @param datedTree a phylo object with branch lengths proportional to absolute time
#' @param method character string; one of "bladj" or "mrbayes"
#' @inheritParams GetMrBayesTree
#' @return A phylo object
MakeUpDates <- function(datedTree = NULL, missingTaxa = NULL, method = "mrbayes", file="mrbayes_makeupdates.nexus"){
	datedTree <- ProcessPhy(input = datedTree)
	datedTree <- CheckPhylo(phylo=datedTree, dated=TRUE)
	method <- match.arg(method, c("bladj", "mrbayes"))
	mrbayes.file <- file
	datedTree <- AddOutgroup(phy=datedTree, outgroup="fake_outgroup", processphy=FALSE)
	ncalibration <- get_node_data(phy=datedTree, node_data=c("node_age", "descendant_tips_label"))
	new.phy <- GetMrBayesTree(phy=datedTree, ncalibration= ncalibration, missingTaxa=missingTaxa, file=mrbayes.file)
	new.phy <- ape::drop.tip(new.phy, "fake_outgroup")
	return(new.phy)
}

#' Checks that missingTaxa argument is ok to be used by MakeMrBayesRunFile function
#' @param missingTaxa A phylo object, a newick character string or a dataframe with taxonomic assignations
#' @return A phylo object, a newick character string or a dataframe with taxonomic assignations
CheckMissingTaxa <- function(missingTaxa){
	if(is.data.frame(missingTaxa)){ # or is.matrix??
		stop("not implemented yet")
		# checkPastisData
		return(missingTaxa)
	}
	if(is.vector(missingTaxa)){
		missingTaxa <- as.character(missingTaxa)
		return(missingTaxa)
	}
	if(is.null(missingTaxa)) {
		return(NULL)
	}
	missingTaxa <- ProcessPhy(missingTaxa)
	if(inherits(missingTaxa, "phylo")){
		stop("not implemented yet")
	}
	tryCatch(CheckPhylo(phylo = missingTaxa), error=function(e) stop("missing taxa must be one of: NULL; a vector with species names; a dataframe with taxonomic assignations; a newick character string; a phylo object"))
	# IMPORTANT: Add a check that taxa in dated.trees is in reference.tree and viceversa
	return(missingTaxa)
	# if (is.null(reference.tree)){
		# if (is.null(missing.taxa)) {
			# cat("specify a reference.tree or missing.taxa to be added to the dated.trees")
			# stop("")
		# } else {
			# # construct a tree with phylotastic or take that from otol?
			# cat("Constructing a reference.tree with taxa from dated.tree and missing.taxa")
		# }
	# } else {
		# if(!is.null(missing.taxa)){
			# cat("A reference.tree was given, missing.taxa argument is ignored")
		# }
}

#' Runs MrBayes from R.
#' @inheritParams GetMrBayesTree
#' @return MrBayes outputs.
MrBayesRun <- function(file=NULL){
	# code borrowed from phyloch::mrbayes()
	if(is.null(file)) stop("You must provide a block file for MrBayes run")
	if (.Platform$OS.type == "unix"){
		system(paste("mb > execute", file))
	} else {
		system(paste("mrbayes ", file, ".bayes", sep = ""))
	}
	tr <- ape::read.nexus(paste(file, ".con.tre", sep = ""))
	return(tr)
}

#' Writes a node calibrations block for a MrBayes run file.
#' @inheritParams GetMrBayesTree
#' @param ncalibrationType A character string specifying the type of calibration. Only "fixed" is implemented for now.
#' @return A set of MrBayes calibration commands printed in console as character strings or as a text file with name specified in file.
#' @export
# This function is set to match node names with constraints obtained from paleotree::GetMrBayesConstraints
GetMrBayesNodeCalibrations <- function(phy=NULL, ncalibration=NULL, ncalibrationType = "fixed", file = NULL){
	phy <- ProcessPhy(input=phy, verbose=FALSE)
    if (!inherits(phy, "phylo")) {
		stop("phy must be a newick character string or in phylo format")
    }
	if(length(ncalibration)==2){  # if it is a list of descendant tips labels and node ages, from get_node_data function
		if(!is.list(ncalibration)) {
			stop("ncalibration must be a newick character string, in phylo format or a list with taxon names and dates")
		}
		includes.ncalibration <- lapply(ncalibration$descendant_tips_label, function(x) gsub(" ", "_", x))
		nages <- ncalibration$node_age

	} else {  # if it is a tree
			ncalibration <- ProcessPhy(input=ncalibration, verbose=FALSE)
			if (!inherits(ncalibration, "phylo")) {
					stop("ncalibration must be a newick character string, a phylo object or a list with taxon names and dates")
	    }
			if(is.null(ncalibration$edge.length) | !ape::is.ultrametric(ncalibration)) {
	    	stop("ncalibration tree must have branch lengths and be ultrametric")
			}
	    phy <- ConvertSpacesToUnderscores(phy)
	    ncalibration <- ConvertSpacesToUnderscores(ncalibration)
	    splits.ncalibration <- ape::prop.part(ncalibration)
	    includes.ncalibration <- lapply(splits.ncalibration, function(x) ncalibration$tip.label[x])
			nages <- ape::branching.times(ncalibration)
	}
	nodes <- sapply(includes.ncalibration, function(tax)
				phytools::findMRCA(phy, tax, type="node")) - length(phy$tip.label)
	calibrations <- paste0("calibrate node", nodes-1, " = ", ncalibrationType, "(", nages, ");")
	root <- which(nodes==1)  # tests for the presence of a root calibration, which should be implemented with treeagepr and not with calibrate
	if(length(root)!=0){
		nodes <- nodes[-root]
		nages <- nages[-root]
		calibrations <- c(calibrations, paste0("prset treeagepr = ", ncalibrationType, "(",
		nages[root], ");"))
	}
    if (!is.null(file)) {
        write(calibrations, file)
    }
    else {
        return(calibrations)
    }
}

#' Identifies the presence of a single lineage outgroup in a phylogeny.
#' @inheritParams GetMrBayesTree
#' @return A character vector with the name of the single lineage outgroup. Returns NA if there is none.
#' @export
IdentifySingletonOutgroup <- function(phy=NULL){
    if (!inherits(phy, "phylo")) {
        phy <- ProcessPhy(input=phy, verbose=FALSE)
    }
	phy <- ConvertSpacesToUnderscores(phy)
    outgroup <- NA
    splits <- ape::prop.part(phy)
    if(length(splits)>1){
    	index <- which.max(sapply(splits, length))
    	s1 <- splits[[index]]
    	index2 <- which.max(sapply(splits[-index], length))
    	s2 <- splits[-index][[index2]]
    	if(length(s1)-length(s2) == 1){
    		outgroup <- phy$tip.label[s1[!(s1 %in% s2)]]
    	}
    }
    return(outgroup)
}

#' Are all desired taxa in the patristic.matrix?
#' @param patristic.matrix A patristic matrix, rownames and colnames must be taxa
#' @param taxa A vector of taxon names
#' @return A Boolean
#' @export
AllMatching <- function(patristic.matrix, taxa) {
	return(sum(!(taxa %in% rownames(patristic.matrix) ))==0)
}


#' Find the index of relevant studies in a opentree_chronograms object
#' @param filtered.results The patristic.matrices that will be used
#' @param cache The cache of studies
#' @return A vector with the indices of studies that have relevant info
#' @export
FindMatchingStudyIndex <- function(filtered.results, cache=get("opentree_chronograms")) {
    return(which(names(cache$trees) %in% names(filtered.results)))
}

#' Return the relevant authors for a set of studies
#' @param results.index A vector from FindMatchingStudyIndex() with the indices of the relevant studies
#' @param cache The cache
#' @return A vector with counts of each author, with names equal to author names
#' @export
TabulateRelevantAuthors <- function(results.index, cache=get("opentree_chronograms")) {
	authors <- cache$authors[results.index]
	return(table(unlist(authors)))
}

#' Return the relevant curators for a set of studies
#' @param results.index A vector from FindMatchingStudyIndex() with the indices of the relevant studies
#' @param cache The cache
#' @return A vector with counts of each curator, with names equal to curator names
#' @export
TabulateRelevantCurators <- function(results.index, cache=get("opentree_chronograms")) {
	curators <- cache$curators[results.index]
	return(table(unlist(curators)))
}

#' Take results.list and process it
#' @param results.list A list returned from using GetSubsetArrayDispatch on opentree_chronograms$trees
#' @param taxa A vector of taxa to match
#' @param partial If TRUE, return matrices that have only partial matches
#' @return A list with the patristic.matrices that are not NA
#' @export
ProcessResultsList <- function(results.list, taxa=NULL, partial=FALSE) {
	if(is.null(taxa)) {
		taxa <- unique(unname(unlist(lapply(final.matrices, rownames))))
	}
	patristic.matrices <- lapply(results.list, "[[", "patristic.matrix.array")

	final.matrices <- patristic.matrices[!is.na(patristic.matrices)]

	if(!partial) {
		final.matrices <- final.matrices[sapply(final.matrices, AllMatching, taxa=taxa)]
	}
	if(length(final.matrices)>0) {
		to.delete <- c()
		for (i in sequence(length(final.matrices))) {
			if(all(is.na(final.matrices[[i]]))) {
				to.delete <- c(to.delete, i)
			}
		}
		if(length(to.delete)>0) {
			final.matrices <- final.matrices[-to.delete]
		}
	}
	return(final.matrices)
}

#Note that originally trees were stored as patristic matrices. This was intended
#to make subsetting fast. The downside is large memory usage. Klaus Schliep wrote
#fast tree subsetting for phylo and multiphylo objects, so now trees are stored
#internally as objects of this type, but with the final output after pruning
#going through patristic matrices.

#Some trees are so large that they can't be stored as patristic distance matrices. For all others,
#patristic matrices are better. For example, for the 20,000 HeathEtAl2012 trees of 35 taxa,
#getting a subset down to two taxa takes 0.0475 seconds just for the pruning, 0.0504 seconds
#for pruning and getting a subset, for a single tree (run times go up linearly with number of trees:
#pruning and converting 1000 trees takes 3 seconds). Subsetting 1000 trees from the patristic
#distance matrix takes just 0.0013 seconds.


#in case we want to cache. Not clear we do
ComputePatristicDistance <- function(phy, test=TRUE,tol=0.01) {
	# stores the distance between taxa
	patristic.matrix <- NA
	if(class(phy)=="phylo") {
		if (test) {
			if (!ape::is.ultrametric(phy,tol)) {
				stop("currently we require that chronograms be ultrametric") # can pad them so that terminals all reach to present
			}
		}
		patristic.matrix<-stats::cophenetic(phy)
	}
	return(patristic.matrix)
}

GetSubsetMatrix <- function(patristic.matrix, taxa, phy4=NULL) {
  #gets a subset of the patristic.matrix. If you give it a phylo4 object, it can check to see if taxa are a clade
  patristic.matrix.new <- patristic.matrix[ rownames(patristic.matrix) %in% taxa,colnames(patristic.matrix) %in% taxa ]
  problem.new <- "none"
  final.size <- sum(rownames(patristic.matrix.new) %in% taxa) # returns number of matches
  if (final.size < length(taxa)) {
    problem.new <- "some of the queried taxa are not on this chronogram, so this is probably an underestimate" # fewer taxa on final matrix than we asked for
    if (final.size < 2 ) {
      problem.new <- "insufficient coverage" # we either have one species or zero. Not enough for an MRCA
      patristic.matrix.new <- NA # to make sure no one uses the zero by mistake
    }
  }
  if(!is.null(phy4)) {
    if (length(phylobase::descendants(phy4, phylobase::MRCA(phy4, taxa), type="tips")) > taxa) {
       problem <- "set of taxa not a clade, so this is probably an overestimate"
    }
  }
  return(list(patristic.matrix= patristic.matrix.new,problem= problem.new))
}

#' Summarize a filtered results list from GetFilteredResults function in various ways
#' @param filtered.results A list of patristic matrices; labels correspond to citations
#' @inheritParams EstimateDates
#' @inherit EstimateDates return details
#' @export
SummarizeResults <- function(filtered.results = NULL, output.format = "phylo.all", input = NULL, partial=TRUE, update_cache = FALSE, cache = get("opentree_chronograms"), showSummary = c("citations", "taxa"), missing.taxa = c("none", "summary", "matrix"), verbose=FALSE) {
		# if(!partial) {
		# 	filtered.results <- filtered.results[which(!sapply(filtered.results, anyNA))]
		# } # not necessary cause already filtered in GetFilteredResults
	if(update_cache){
		cache <- UpdateCache(save = TRUE, verbose=verbose)
	}
	if(is.null(filtered.results) | !is.list(filtered.results)){
		stop("filtered.results argument must be a list from GetFilteredResults function.")
	}
	output.format.in <- match.arg(output.format, choices=c("citations", "mrca", "newick.all", "newick.sdm", "newick.median", "phylo.sdm", "phylo.median", "phylo.median", "phylo.all", "html", "data.frame"))
	missing.taxa.in <- match.arg(missing.taxa, choices=c("none", "summary", "matrix"))
	showSummary.in <- match.arg(showSummary, c("citations", "taxa", "none"), several.ok=TRUE)
	if(is.null(input)){
		input.in <- unique(rapply(filtered.results, rownames))
		if(missing.taxa.in !="none") warning("showing taxon summary from taxa found in at least one chronogram, this excludes input taxa not found in any chronogram")
	} else {
		if(!is.character(input)) stop("input must be a character vector")
		input.in <- input
	}
	results.index <- FindMatchingStudyIndex(filtered.results, cache)
	return.object <- NA
	input.match <- unique(rapply(filtered.results, rownames))
	# if(any(!input.match %in% input)) warning("input does not contain all or any taxa from filteredresults object")
	absent.input <- input.in[!input.in %in% input.match]
	if(length(absent.input)<=0) {
		if(is.null(input)){
			absent.input <- "NULL"
		} else {
			absent.input <- "None"
		}
	}
	if(missing.taxa.in == "matrix"){
		missing.taxa.list <- vector(mode="list")
		# tax <- unique(rapply(filtered.results, rownames)) #rownames(filtered.results[[1]])
		for(result.index in sequence(length(filtered.results))){
			n <- rownames(filtered.results[[result.index]])
			m <- match(input.match,n)
			missing.taxa.list[[result.index]] <- n[m]
		}
		missing.taxa.matrix <- do.call(rbind, missing.taxa.list) #transforms a list of names into a matrix of names
		missing.taxa.matrix <- !is.na(missing.taxa.matrix) # makes a boolean matrix
		colnames(missing.taxa.matrix) <- input.match
		rownames(missing.taxa.matrix) <- sequence(nrow(missing.taxa.matrix))
	}
	if(missing.taxa.in == "summary" | any(grepl("taxa", showSummary.in))){ # may add here another condition: | makeup_brlen ==TRUE
		# tax <- unique(rapply(filtered.results, rownames)) #rownames(filtered.results[[1]])
		x <- rapply(filtered.results, rownames)
		prop <- c()
		for (taxon in input.match){
			prop <- c(prop, paste(length(which(taxon==x)), "/", length(filtered.results), sep=""))
		}
		missing.taxa.summary <- data.frame(Taxon=input.match, Chronograms=prop)
	}
	if(output.format.in == "citations") {
		return.object <- names(filtered.results)
	}
	if(output.format.in == "mrca") {
		return.object <- GetAges(filtered.results, partial=partial)
	}
	if(output.format.in == "newick.all") {
		trees <- sapply(filtered.results, PatristicMatrixToNewick)
		return.object <- trees[which(!is.na(trees))]
	}
	if(output.format.in == "newick.sdm") {
		local.results <- RunSDM(filtered.results)
		filtered.results <- local.results$filtered.results
		tree <- local.results$phy
		return.object <- ape::write.tree(tree)
	}
	if(output.format.in == "phylo.sdm") {
		local.results <- RunSDM(filtered.results)
		filtered.results <- local.results$filtered.results
		tree <- local.results$phy
		return.object <- tree
	}
	if(output.format.in == "newick.median") {
		patristic.array <- BindMatrices(filtered.results)
		median.matrix <- SummaryPatristicMatrixArray(patristic.array)
		tree <- PatristicMatrixToNewick(median.matrix)
		return.object <- tree
	}
	if(output.format.in == "phylo.median") {
		patristic.array <- BindMatrices(filtered.results)
		median.matrix <- SummaryPatristicMatrixArray(patristic.array)
		tree <- PatristicMatrixToTree(median.matrix)
		return.object <- tree
	}
	if(output.format.in == "phylo.all") {
		trees <- lapply(filtered.results, PatristicMatrixToTree)
		return.object <- trees[which(!is.na(trees))]
	}
	if(output.format.in == "html") {
		out.vector1 <- "<table border='1'><tr><th>MRCA Age (MY)</th><th>Ntax</th><th>Citation</th><th>Newick"
		if(missing.taxa.in == "matrix"){
			out.vector1 <- paste(out.vector1, paste("</th><th>", colnames(missing.taxa.matrix), sep="", collapse=""), sep="")
		}
		out.vector1 <- paste(out.vector1, "</th></tr>", sep="")
		ages <- GetAges(filtered.results, partial=partial)
		trees <- sapply(filtered.results, PatristicMatrixToNewick)
		out.vector2 <- c()
		for(result.index in sequence(length(filtered.results))) {
			out.vector2 <- paste(out.vector2, "<tr><td>",ages[result.index],"</td><td>",sum(!is.na(diag(filtered.results[[result.index]]))), "</td><td>", names(filtered.results)[result.index], "</td><td>", trees[result.index],  sep="")
			if(missing.taxa.in == "matrix"){
				out.vector2 <- paste(out.vector2, paste("</td><td>", missing.taxa.matrix[result.index,], sep="", collapse=""), sep="")
			}
			out.vector2 <- paste(out.vector2, "</td></tr>", sep="")
		}
		out.vector <- paste(out.vector1, out.vector2, "</table>", sep="")
		if(missing.taxa.in == "summary"){
			missing.taxa.summ <- as.matrix(missing.taxa.summary)
			out.vector3 <- "<p></p><table border='1'><tr><th>Taxon</th><th>Chronograms</th><tr>"
			for (summary.index in sequence(nrow(missing.taxa.summ))){
				out.vector3 <- paste(out.vector3, paste("</td><td>", missing.taxa.summ[summary.index,], sep="", collapse=""), "</td></tr>", sep="")
			}
			out.vector <- paste(out.vector, out.vector3, "</table>", sep="")
		}
		# out.vector4 <- c()
		if(missing.taxa.in != "none") {
			out.vector4 <- "<p></p><table border='1'><tr><th> </th><th>Absent Taxa</th><tr>"
			for (i in 1:length(absent.input)){
				out.vector4 <- paste(out.vector4, "<tr><td>", i, "</td><td>", absent.input[i], "</td><tr>", sep="")
			}
			out.vector4 <- paste(out.vector4, "</table>", sep="")
			out.vector <- paste(out.vector, out.vector4, sep="")
		}
		return.object <- out.vector
	}
	if(output.format.in == "data.frame") {
		out.df <- data.frame()
		ages <- GetAges(filtered.results, partial=partial)
		trees <- sapply(filtered.results, PatristicMatrixToNewick)
		for(result.index in sequence(length(filtered.results))) {
			out.line<- data.frame(Age=ages[result.index],Ntax=sum(!is.na(diag(filtered.results[[result.index]]))), Citation=names(filtered.results)[result.index], Newick= trees[result.index])
			if(result.index == 1) {
				out.df <- out.line
			} else {
				out.df <- rbind(out.df, out.line)
			}
		}
		if(missing.taxa.in == "matrix"){
			out.df <- cbind(out.df, missing.taxa.matrix)
		}
		rownames(out.df) <- NULL
		return.object <- out.df
	}
	if(missing.taxa.in != "none" & output.format.in !="html"){
		return.object <- list(return.object)
		if(missing.taxa.in=="matrix") {
			if(output.format.in !="data.frame") return.object <- c(return.object, list(missing.taxa=missing.taxa.matrix))
		}
		if(missing.taxa.in=="summary") {
			return.object <- c(return.object, list(missing.taxa=missing.taxa.summary))
		}
		return.object <- c(return.object, list(absent.taxa=data.frame(Taxon=absent.input)))
		names(return.object)[1] <- output.format.in
		if(output.format.in =="data.frame"){
			names(return.object)[1] <- "results"
			if(missing.taxa.in == "matrix") names(return.object)[1] <- "results.and.missing.taxa"
		}
	}

	if(any("citations" %in% showSummary.in) & !any(output.format.in %in% c("citations", "html", "data.frame"))) {
		if(output.format.in == "citations"){
			cat("Target taxa found in trees from:", "\n")
			print(names(filtered.results), quote=FALSE)
			cat("\n")
		} else {
			cat("Source chronograms from:", "\n")
			print(names(filtered.results), quote=FALSE)
			cat("\n")
		}
	}
	if(any(grepl("taxa", showSummary.in)) & missing.taxa.in!="summary") {
		cat("Target taxa presence in source chronograms:", "\n")
		print(missing.taxa.summary)
		cat("\n")
		cat("Target taxa completely absent from source chronograms:", "\n")
		print(data.frame(Taxon=absent.input))
		cat("\n")
	}
	return(return.object)
}


#' Figure out which subset function to use
#' @param study.element The thing being passed in: an array or a phylo to serve as reference
#' @param taxa Vector of taxon names to get a subset for
#' @param phy A user tree to congruify in phylo format (ape)
#' @param phy4 A user tree to congruify in phylo4 format (phylobase)
#' @param method Which method to use for congruification
#' @return A patristic matrix with for the taxa.
#' @export
GetSubsetArrayDispatch <- function(study.element, taxa, phy=NULL, phy4=NULL, method="PATHd8") {
  if(class(study.element)=="array") {
    return(GetSubsetArrayBoth(study.element, taxa, phy, phy4, method))
  } else {
    return(GetSubsetArrayBothFromPhylo(reference.tree.in = study.element, taxa.in = taxa, phy.in = phy, phy4.in = phy4, method.in = method))
  }
}

GetSubsetArrayBothFromPhylo <- function(reference.tree.in, taxa.in, phy.in=NULL, phy4.in=NULL, method.in="PATHd8") {
#COMMENTING OUT: OpenTree gives single trees, let's just standardize on those
#  if (class(reference.tree)=="phylo") {
#    reference.tree<-c(reference.tree) #from here in, assumes multiphylo object, even if a single tree
#  }
	congruify=FALSE
	if(!is.null(phy.in[1])) {
		congruify=TRUE
		if(is.na(phy.in[1])) {
			congruify=FALSE
		}
	}
  if (!congruify) {
    return(GetSubsetArrayFromPhylo(reference.tree=reference.tree.in, taxa=taxa.in, phy4=phy4.in, method.in))
  }
  else { #congruify
    return(GetSubsetArrayCongruifyFromPhylo(reference.tree=reference.tree.in, taxa=taxa.in, phy=phy.in, method.in))
  }

}

PruneTree <- function(phy, taxa) {
	return(ape::drop.tip(phy, tip=phy$tip.label[-(which(phy$tip.label %in% taxa))]))
}

GetSubsetArrayFromPhylo <- function(reference.tree, taxa, phy4=NULL, method="PATHd8") {
  final.size<-sum(reference.tree$tip.label %in% taxa) # returns number of matches
  if(final.size>=2) { #it's worth doing the pruning
    reference.tree<-PruneTree(reference.tree, taxa)
    #reference.tree<-pruneTrees(reference.tree, taxa) #pruneTrees is the new, fast fn from Klaus Schliep. Eventually will be in phangorn, currently in datelife2
  }
  problem <- "none"
  patristic.matrix.array <- NA
  if (final.size < length(taxa)) {
    problem <- "missing some taxa on chronogram, so this is probably an underestimate" # fewer taxa on final matrix than we asked for
    if (final.size < 2 ) {
      problem <- "insufficient coverage" # we either have one species or zero. Not enough for an MRCA
      patristic.matrix.array <- NA # to make sure no one uses the zero by mistake
    }
  }
  if (final.size >= 2) {
  	patristic.matrix.array <- ComputePatristicDistance(reference.tree)
  }
  if(!is.null(phy4)) {
    if (length(phylobase::descendants(phy4, phylobase::MRCA(phy4, taxa), type="tips")) > taxa) {
      problem <- "set of taxa not a clade, so this is probably an overestimate"
    }
  }
  return(list(patristic.matrix.array=patristic.matrix.array,problem=problem))
}

GetSubsetArrayCongruifyFromPhylo <- function(reference.tree, taxa, phy=NULL, method="PATHd8") {
  final.size<-sum(reference.tree$tip.label %in% taxa) # returns number of matches
  if(final.size>=2) { #it's worth doing the pruning
   reference.tree<-PruneTree(reference.tree, taxa)
   #reference.tree<-pruneTrees(reference.tree, taxa) #pruneTrees is the new, fast fn from Klaus Schliep. Eventually will be in phangorn, currently in datelife2
  }
  problem.new <- "none"
  patristic.matrix.array.new <- NA
  if (final.size < length(taxa)) {
    problem.new <- "missing some taxa on chronogram, so this is probably an underestimate" # fewer taxa on final matrix than we asked for
    if (final.size < 3 ) {
      problem.new <- "insufficient coverage" # we either have one species or zero. Not enough for an MRCA
      patristic.matrix.array.new <- NA # to make sure no one uses the zero by mistake
      return(list(patristic.matrix.array= patristic.matrix.array.new,problem= problem.new))

    }
  }
  if (final.size >= 3) {
  	patristic.matrix.array.new <-   CongruifyTreeFromPhylo(reference.tree, query.tree=phy, method=method)
  }
  return(list(patristic.matrix.array= patristic.matrix.array.new,problem= problem.new))
}


GetSubsetArrayBoth <- function(patristic.matrix.array, taxa, phy=NULL, phy4=NULL, method="PATHd8") {
  if (is.null(phy)) {
    return(GetSubsetArray(patristic.matrix.array=patristic.matrix.array, taxa=taxa, phy4=phy4))
  }
  else { #congruify
    return(GetSubsetArrayCongruify(patristic.matrix.array=patristic.matrix.array, taxa=taxa, phy=phy, method))
  }
}

CongruifyTree <- function(patristic.matrix, query.tree, method="PATHd8", attempt.fix=TRUE) {
  result.matrix<-matrix(nrow=dim(patristic.matrix)[1], ncol=dim(patristic.matrix)[2])
  if(is.null(query.tree$edge.length)) {
    query.tree$edge.length<-numeric(nrow(query.tree$edge))
  }
#	try(result.matrix<-ComputePatristicDistance(ConvertUnderscoresToSpaces(geiger::congruify.phylo(ConvertSpacesToUnderscores(PatristicMatrixToTree(patristic.matrix)), ConvertSpacesToUnderscores(query.tree), NULL, 0, scale=method)$phy)))
try(result.matrix<-ComputePatristicDistance(CongruifyAndCheck(reference=PatristicMatrixToTree(patristic.matrix), target=query.tree, scale=method, attempt.fix=attempt.fix)))
  return(result.matrix)
}

CongruifyTreeFromPhylo <- function(reference.tree, query.tree, method="PATHd8", attempt.fix=TRUE) {
  result.matrix<-matrix(nrow=ape::Ntip(reference.tree), ncol=ape::Ntip(reference.tree))
  if(is.null(query.tree$edge.length)) {
    query.tree$edge.length<-numeric(nrow(query.tree$edge)) #makes it so that branches that don't match reference tree get zero length
  }
	try(result.matrix<-ComputePatristicDistance(CongruifyAndCheck(reference=reference.tree, target=query.tree, scale=method, attempt.fix=attempt.fix)))
  return(result.matrix)
}

CongruifyAndCheck <- function(reference, target, taxonomy=NULL, tol=0.01, scale="pathd8", attempt.fix=TRUE) {
  if(!ape::is.ultrametric(reference, tol=tol)) {
    return(NA)
  }
	new.tree <- ConvertUnderscoresToSpaces(suppressWarnings(geiger::congruify.phylo(ConvertSpacesToUnderscores(reference), ConvertSpacesToUnderscores(target), taxonomy=taxonomy, tol=tol, scale=scale)$phy)) #suppressing warnings b/c geiger ignores tolerance
	if(anyNA(new.tree$edge.length) & attempt.fix) {
		warning("Congruification resulted in NA edge lengths. Resolving polytomies and making up starting branch lengths")
		new.tree <- ConvertUnderscoresToSpaces(geiger::congruify.phylo(ConvertSpacesToUnderscores(reference), ConvertSpacesToUnderscores(ape::compute.brlen(ape::multi2di(target))), taxonomy, tol, scale)$phy)
		if(anyNA(new.tree$edge.length)) {
			new.tree <- NA
		}
	}
	new.tree$edge.length[which(new.tree$edge.length<0)] <- 0 #sometimes pathd8 returns tiny negative branch lengths. https://github.com/phylotastic/datelife/issues/11
	return(new.tree)
}

#' Convert spaces to underscores in trees
#' @param phy A phylo object
#' @return A phylo object
#' @export
ConvertSpacesToUnderscores <- function(phy) {
	phy$tip.label <- gsub(" ", "_", phy$tip.label)
	return(phy)
}

#' Convert underscores to spaces in trees
#' @param phy A phylo object
#' @return A phylo object
#' @export
ConvertUnderscoresToSpaces <- function(phy) {
	phy$tip.label <- gsub("_", " ", phy$tip.label)
	return(phy)
}

GetSubsetArrayCongruify <- function(patristic.matrix.array, taxa, phy=NULL, method="PATHd8") {
  #gets a subset of the patristic.matrix.array.
  patristic.matrix.array <- patristic.matrix.array[ rownames(patristic.matrix.array) %in% taxa,colnames(patristic.matrix.array) %in% taxa,  ]
  problem <- "none"
  final.size <- sum(rownames(patristic.matrix.array) %in% taxa) # returns number of matches
  if (final.size < length(taxa)) {
    problem <- "missing some taxa on chronogram, so this is probably an underestimate" # fewer taxa on final matrix than we asked for
    if (final.size < 3 ) {
      problem <- "insufficient coverage" # we either have one species or zero. Not enough for an MRCA
      patristic.matrix.array <- NA # to make sure no one uses the zero by mistake
      return(list(patristic.matrix.array=patristic.matrix.array,problem=problem))

    }
  }
  patristic.matrix.list <- SplitArray(patristic.matrix.array)
  patristic.matrix.array<-BindMatrices(lapply(patristic.matrix.list, CongruifyTree, query.tree=phy, scale=method)) #yes, this should be parallel
  return(list(patristic.matrix.array=patristic.matrix.array,problem=problem))
}

GetSubsetArray <- function(patristic.matrix.array, taxa, phy4=NULL) {
  #gets a subset of the patristic.matrix.array. If you give it a phylo4 object, it can check to see if taxa are a clade
  patristic.matrix.array <- patristic.matrix.array[ rownames(patristic.matrix.array) %in% taxa,colnames(patristic.matrix.array) %in% taxa,  ]
  problem <- "none"
  final.size <- sum(rownames(patristic.matrix.array) %in% taxa) # returns number of matches
  if (final.size < length(taxa)) {
    problem <- "missing some taxa on chronogram, so this is probably an underestimate" # fewer taxa on final matrix than we asked for
    if (final.size < 2 ) {
      problem <- "insufficient coverage" # we either have one species or zero. Not enough for an MRCA
      patristic.matrix.array <- NA # to make sure no one uses the zero by mistake
    }
  }
  if(!is.null(phy4)) {
    if (length(phylobase::descendants(phy4, phylobase::MRCA(phy4, taxa), type="tips")) > taxa) {
      problem <- "set of taxa not a clade, so this is probably an overestimate"
    }
  }
  return(list(patristic.matrix.array=patristic.matrix.array,problem=problem))
}

#' Get time of MRCA from patristic matrix
#' @param patristic.matrix A patristic matrix
#' @param partial If TRUE, drop NA from the patristic matrix; if FALSE, will return NA if there are missing entries
#' @return The depth of the MRCA
#' @export
GetAge <- function(patristic.matrix, partial=TRUE) {
  # 0.5 since patristic distance is down to the root and back up
  return(0.5 * max(patristic.matrix, na.rm=partial))
}

#' Get vector of MRCA from list of patristic matrices
#' @param filtered.results List of patristic matrices
#' @param partial If TRUE, drop NA from the patristic matrix; if FALSE, will return NA if there are missing entries
#' @return Vector of MRCA ages with names same as in filtered.results
#' @export
GetAges <- function(filtered.results, partial=TRUE) {
	ages <- sapply(filtered.results, GetAge, partial=partial)
	return(ages)
}

#' Function to reorder a matrix so that row and column labels are in alphabetical order
#' @param patristic.matrix A patristic matrix with row and column names for taxa
#' @return patristic.matrix A patristic matrix with row and column names for taxa in alphabetial order
#' @export
ReorderMatrix <- function(patristic.matrix) {
  return(patristic.matrix[order(rownames(patristic.matrix)),order(colnames(patristic.matrix))])
}

#' Function to fill in empty cells in a patristic matrix for missing taxa
#' @param patristic.matrix A patristic matrix with row and column names for taxa
#' @param all.taxa A vector of the names of all taxa you want, including ones not in the patristic matrix
#' @return Patristic.matrix for all.taxa, with NA for entries between taxa where at least one was not in the original patristic.matrix
#' @export
PadMatrix <- function(patristic.matrix, all.taxa) {
	number.missing <- length(all.taxa) - dim(patristic.matrix)[1]
	final.matrix <- patristic.matrix
	if(number.missing>0) {
		final.matrix <- rbind(patristic.matrix, matrix(nrow=number.missing, ncol=dim(patristic.matrix)[2]))
		final.matrix <- cbind(final.matrix, matrix(ncol=number.missing, nrow=dim(final.matrix)[1]))
		rownames(final.matrix) <- c(rownames(patristic.matrix), all.taxa[-which(all.taxa %in% rownames(patristic.matrix))])
	 	colnames(final.matrix) <- c(colnames(patristic.matrix), all.taxa[-which(all.taxa %in% colnames(patristic.matrix))])
	}
 	return(ReorderMatrix(final.matrix))
}

#' Function to remove missing taxa
#' @param patristic.matrix A patristic matrix with row and column names for taxa
#' @return Patristic.matrix for all.taxa
#' @export
UnpadMatrix <- function(patristic.matrix) {
	bad.ones <- which(apply(is.na(patristic.matrix),2,all))
	if(length(bad.ones)>0) {
		patristic.matrix <- patristic.matrix[-bad.ones, -bad.ones]
	}
	return(patristic.matrix)
}

TestNameOrder <- function(patristic.matrix, standard.rownames, standard.colnames) {
  if (compare::compare(rownames(patristic.matrix),standard.rownames)$result!=TRUE) {
    return(FALSE)
  }
  if (compare::compare(colnames(patristic.matrix),standard.colnames)$result!=TRUE) {
    return(FALSE)
  }
  return(TRUE)
}

#' Convert list of patristic matrices to a 3D array
#' @param patristic.matrix.list List of patristic matrices
#' @param pad If TRUE, pad missing entries
#' @return A 3d array of patristic matrices
#' @export
BindMatrices <- function(patristic.matrix.list, pad=TRUE) {
  all.taxa <- sort(unique(unname(unlist(lapply(patristic.matrix.list, rownames)))))
  if(pad) {
    patristic.matrix.list <- lapply(patristic.matrix.list, PadMatrix, all.taxa=all.taxa)
  }
  original.size<-length(patristic.matrix.list)
  patristic.matrix.list<-lapply(patristic.matrix.list,ReorderMatrix)
  if(length(patristic.matrix.list)<1) {
    stop(paste("The patristic matrices you are trying to bind are too few; input was ", original.size, " and current length is ", length(patristic.matrix.list), sep=""))
  }
  standard.rownames<-rownames(patristic.matrix.list[[1]])
  standard.colnames<-colnames(patristic.matrix.list[[1]])
  matching.names<-sapply(patristic.matrix.list,TestNameOrder,standard.rownames,standard.colnames)
  if (sum(matching.names)!=length(matching.names)) {
    stop("The patristic matrices you are trying to bind do not have the same taxa")
  }
  return(abind::abind(patristic.matrix.list, along=3 ))
}

SplitArray <- function(patristic.matrix.array) {
  return(lapply(sequence(dim(patristic.matrix.array)[3]),AsubForLapply,patristic.matrix.array))
}

AsubForLapply <- function(idx, x, dims=3) {
  return(abind::asub(x, idx, dims))
}

GetQuantiles <- function(ages,probs=c(0.5,0,0.025,0.975,1) ) {
  # just utility wrapper function with different defaults
  return(stats::quantile(ages,probs))
}

VectorToTableRow <- function(x,digits=2) {
  return(paste(paste("<td>",round(x,digits),sep=""),"</td>",sep="",collapse=""))
}

#' Convert patristic matrix to a phylo object
#' @param patristic.matrix A patristic matrix
#' @return A rooted phylo object
#' @export
PatristicMatrixToTree <- function(patristic.matrix) {
  if(anyNA(patristic.matrix)) {
  	patristic.matrix <- patristic.matrix[rowSums(is.na(patristic.matrix)) != ncol(patristic.matrix),colSums(is.na(patristic.matrix)) != nrow(patristic.matrix)]
  }
  if(dim(patristic.matrix)[1] < 2) {
  	return(NA)
  }
	tree <- NA
	if(dim(patristic.matrix)[1] == 2) {
		tree <- ape::rtree(n=2, rooted=TRUE, tip.label=rownames(patristic.matrix), br=patristic.matrix[1,2]/2)
	} else {
  	tree <- ape::nj(patristic.matrix)
	}
  if(ape::Ntip(tree)>2) {
    tree <- 	phangorn::midpoint(tree)
  }
	if(length(which(tree$edge.length<0))>0) {
		warning(paste("Converting from patristic distance matrix to a tree resulted in some negative branch lengths; the largest by magnitude was", min(tree$edge.length)))
		tree$edge.length[which(tree$edge.length<0)] <- 0 #sometimes NJ returns tiny negative branch lengths. https://github.com/phylotastic/datelife/issues/11
	}
  return(tree)
}

#' Convert patristic matrix to a newick string
#' @param patristic.matrix A patristic matrix
#' @return A newick string
#' @export
PatristicMatrixToNewick <- function(patristic.matrix) {
  tree <- PatristicMatrixToTree(patristic.matrix)
  if(class(tree)=="phylo") {
  	return(ape::write.tree(tree))
  }
  return(NA)
}


#' Summarize patristic matrix array (by default, median)
#' @param patristic.matrix.array 3D array of patristic matrices
#' @param fn The function to use ot summarize
#' @return A 2d array with the median (or max, or mean, etc) of the input array
#' @export
SummaryPatristicMatrixArray <- function(patristic.matrix.array,fn=stats::median) {
  return(apply(patristic.matrix.array,MARGIN=c(1,2),fn, na.rm=TRUE))
}

SamplePatristicMatrix <- function(patristic.matrix.array, uncertainty) {
 # if (dim(patristic.matrix.array)[3] == 1) {
 # 	patristic.matrix<-patristic.matrix.array[,,1]
 # 	#need order of node depths, from just the upper triangular and diagonal part of the matrix
 # 	element.order<-order(patristic.matrix[upper.tri(patristic.matrix,diag=FALSE)],decreasing=TRUE)
 # 	new.patristic.matrix<-patristic.matrix*0
 # 	cur.val<-patristic.matrix[upper.tri(patristic.matrix,diag=FALSE)][element.order[1]]
 #   new.patristic.matrix[upper.tri(new.patristic.matrix,diag=FALSE)][element.order[1]] <- cur.val + runif(1, -cur.val*uncertainty/100, cur.val*uncertainty/100)
#	element.order<-element.order[-1]
#  	for (i in sequence(length(element.order))) {
#  		cur.val<-patristic.matrix[upper.tri(patristic.matrix,diag=FALSE)][element.order[i]]
#  		new.patristic.matrix[upper.tri(new.patristic.matrix,diag=FALSE)][element.order[i]] <- cur.val + runif(1, -cur.val*uncertainty/100, min(cur.val*uncertainty/100, min( ))
#  	}
#  }
#  else {
  	return(patristic.matrix <- patristic.matrix.array[,,sample.int(1, size = dim(patristic.matrix.array)[3] )] )
 # }
}


#' Get all calibrations given trees in database
#' @param input vector of names, a newick string, or a phylo object
#' @param partial Boolean; default TRUE: use source trees even if they only match some of the desired taxa
#' @param usetnrs Boolean; default False. If TRUE, use OpenTree's services to resolve names. This can dramatically improve the chance of matches, but also take much longer
#' @param approximatematch Boolean; default TRUE: use a slower TNRS to correct mispellings, increasing the chance of matches (including false matches)
#' @param cache The cached set of chronograms and other info from data(opentree_chronograms)
#' @inheritParams EstimateDates
#' @return data.frame of calibrations
#' @export
GetAllCalibrations <- function(input = c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"), partial = TRUE, usetnrs = FALSE, approximatematch = TRUE, update_cache = FALSE, cache = get("opentree_chronograms"), verbose=FALSE) {
	phylo.results <- EstimateDates(input = input, partial = partial, usetnrs = usetnrs, approximatematch = approximatematch, update_cache = update_cache, cache = cache, output.format = "phylo.all", verbose=verbose)
	constraints.df <- data.frame()
	for (i in sequence(length(phylo.results))) {
		local.df <- geiger::congruify.phylo(reference = phylo.results[[i]], target = phylo.results[[i]], scale = NA)$calibrations
		local.df$reference <- names(phylo.results)[i]
		if(i == 1) {
			constraints.df <- local.df
		} else {
			constraints.df <- rbind(constraints.df, local.df)
		}
	}
	return(constraints.df)
}

#' Use all calibrations given trees in database to date a tree
#' @param phy A phylo object
#' @param partial If TRUE, use source trees even if they only match some of the desired taxa
#' @param usetnrs If TRUE, use OpenTree's services to resolve names. This can dramatically improve the chance of matches, but also take much longer
#' @param approximatematch If TRUE, use a slower TNRS to correct mispellings, increasing the chance of matches (including false matches)
#' @param cache The cached set of chronograms and other info from data(opentree_chronograms)
#' @param expand How much to expand by each step to get consistent calibrations
#' @param giveup How many expansion to try before giving up
#' @inheritParams EstimateDates
#' @return list with chronogram, original calibrations, and expanded calibrations
#' @export
#' @details
#' This will try to use the calibrations as fixed ages.
#' If that fails (often due to conflict between calibrations), it will expand the range of the minage and maxage and try again. And repeat.
#' expand sets the expansion value: should be between 0 and 1
UseAllCalibrations <- function(phy = GetBoldOToLTree(c("Rhea americana",  "Struthio camelus", "Gallus gallus"), chronogram = FALSE, verbose=FALSE), partial = TRUE, usetnrs = FALSE, approximatematch = TRUE, update_cache = FALSE, cache = get("opentree_chronograms"), expand = 0.1, giveup = 100, verbose=FALSE) {
	# if(!geiger::is.phylo(phy)){
	# 	cat("phy argument must be a phylo object.", "\n")
	# 	stop()
	# }
	#	[code above] useful if we decide to allow newick as input, would need to add some newick to phylo code in here... OK! It is now implemented in ProcessPhy function:
	input <- ProcessPhy(input=phy, verbose=verbose)
	calibrations.df <- GetAllCalibrations(input = gsub('_', ' ', phy$tip.label), partial = partial, usetnrs = usetnrs, approximatematch = approximatematch, update_cache = update_cache, cache = cache, verbose=verbose)
	phy$tip.label <- gsub(' ', '_', phy$tip.label) #underscores vs spaces: the battle will never end.
	calibrations.df$taxonA <- gsub(' ', '_', calibrations.df$taxonA)
	calibrations.df$taxonB <- gsub(' ', '_', calibrations.df$taxonB)
	calibrations.df <- calibrations.df[which(calibrations.df$taxonA %in% phy$tip.label),]
	calibrations.df <- calibrations.df[which(calibrations.df$taxonB %in% phy$tip.label),]
	original.calibrations.df <- calibrations.df
	chronogram <- NULL
	try(chronogram <- geiger::PATHd8.phylo(phy, calibrations.df), silent=TRUE)
	if(!is.null(chronogram)) {
		chronogram$edge.length[which(chronogram$edge.length<0)] <- 0 #sometimes pathd8 returns tiny negative branch lengths. https://github.com/phylotastic/datelife/issues/11
	}
	attempts = 0
	if(expand != 0) {
		while(is.null(chronogram) & attempts<giveup) {
			calibrations.df <- original.calibrations.df
			calibrations.df$MaxAge <- calibrations.df$MaxAge * ((1+expand)^attempts)
			calibrations.df$MinAge <- calibrations.df$MinAge * ((1-expand)^attempts)

			# We will have no fixed ages. Pathd8 just quietly gives up. So instead, we add a tiny branch with a zero calibration
			# between it and its sister.
			made.up.edgelength <- min(1e-9, .001*min(phy$edge.length))
			phy2 <- phytools::bind.tip(ape::reorder.phylo(phy), "tinytip", edge.length = made.up.edgelength, where = 1, position = made.up.edgelength) #bind tip has weird behavior for non-reordered trees
			calibrations.df[dim(calibrations.df)[1]+1,]<- c("fixed", 0, 0, phy$tip.label[1], "tinytip", "none")
			try(chronogram <- geiger::PATHd8.phylo(phy2, calibrations.df))
			if(!is.null(chronogram)) {
				chronogram$edge.length[which(chronogram$edge.length < 0)] <- 0 #sometimes pathd8 returns tiny negative branch lengths. https://github.com/phylotastic/datelife/issues/11
				chronogram <- ape::drop.tip(chronogram, "tinytip")
			}
			attempts <- attempts+1
		}
		if(attempts > 0) {
			# print("Dates are even more approximate than usual: had to expand constraints to have them agree")
			cat("Dates are even more approximate than usual: had to expand constraints to have them agree", "\n")
		}
	}
	return(list(phy = chronogram, calibrations.df = calibrations.df, original.calibrations.df = original.calibrations.df))
}

#' Use Barcode of Life data to get branch lengths on the OToL tree of a set of taxa.
#' @inheritParams EstimateDates
#' @param otol_version Version of OToL to use
#' @param chronogram Boolean; default to TRUE:  branch lengths represent time estimated with ape::chronoMPL. If FALSE, branch lengths represent relative substitution rates estimated with phangorn::acctran.
#' @param doML Boolean; if TRUE, does ML branch length optimization with phangorn::optim.pml
#' @inheritParams ProcessInput
#' @return A phylogeny with ML branch lengths
#' @export
GetBoldOToLTree <- function(input = c("Rhea americana",  "Struthio camelus", "Gallus gallus"), usetnrs = FALSE, approximatematch = TRUE, marker = "COI", otol_version = "v2", chronogram = TRUE, doML = FALSE, sppfromtaxon=FALSE, verbose=FALSE) {
	#otol returns error with missing taxa in v3 of rotl
	input <- CheckInput(input=input, usetnrs = usetnrs, approximatematch = approximatematch, sppfromtaxon = sppfromtaxon, verbose=verbose)
	input <- input$cleaned.names
	if (verbose) cat("Searching", marker, "sequences for these taxa in BOLD...", "\n")
	sequences <- bold::bold_seqspec(taxon = input, marker = marker)
	sequences$nucleotide_ATGC <- gsub("[^A,T,G,C]", "", sequences$nucleotides)  # preserve good nucleotide data, i.e., only A,T,G,C
	sequences$nucleotide_ATGC_length <- unlist(lapply(sequences$nucleotide_ATGC, nchar))  # add a column in data.frame, indicating the amount of good information contained in sequences#nucelotides (ATGC)
	if(length(sequences) == 1) {  # it is length == 82 when there is at least 1 sequence available, if this is TRUE, it means there are no sequences in BOLD for the set of input taxa.
		if (verbose) cat("No sequences found in BOLD for input taxa...", "\n")
		# if (!usetnrs) cat("Setting usetnrs=TRUE might change this, but it can be slowish.", "\n")
		warning("Names in input do not match BOLD specimen records. No tree was constructed.")
		return(NA)
	}
	if (verbose) cat("\t", "OK.", "\n")
	rr <- rotl::tnrs_match_names(names = input)
	rr <- rr[!is.na(rr$unique_name),]  # gets rid of names not matched with rotl::tnrs_match_names; otherwise rotl::tol_induced_subtree won't run
	phy <- ape::multi2di(rotl::tol_induced_subtree(ott_ids=rr$ott_id, label_format = "name",  otl_v = otol_version))
	phy$tip.label <- gsub("_ott.*","", phy$tip.label)
	# when there are synonyms among the input names, phy will conserve the accepted name (rr$uniqe_name) instead of the original query name from input (rr$search_string)
	# this produces an error downstream, while using phangorn::pml()
	# to avoid this error, we replace the unique name by the original query name in phy$tip.label:
	mm <- match(phy$tip.label, gsub(" ","_", rr$unique_name))  # this gets the order of tip labels in phy
	phy$tip.label <- gsub(" ","_", input[mm])  # this overlaps the original query over phy$tip.labels in the correct order
	final.sequences <- matrix("-", nrow = length(input), ncol = max(sapply(strsplit(sequences$nucleotides, ""), length)))
	final.sequences.names <- rep(NA, length(input))
	# for (i in sequence(dim(sequences)[1])) {
		# taxon <- sequences$species_name[i]
		# if(!(taxon %in% final.sequences.names)) {
			# seq <- strsplit(sequences$nucleotide[i],"")[[1]]
			# matching.index <- 1+sum(!is.na(final.sequences.names))
			# final.sequences[matching.index, sequence(length(seq))] <- seq
			# final.sequences.names[matching.index] <- taxon
		# }
	# }
	row.index <- 0
	taxa.to.drop <- c()
	for (i in input){
		row.index <- row.index + 1
		taxon.index <- which(grepl(i, sequences$species_name))  # what happens here if there are no sequences from the taxon????
		if (length(taxon.index)>0){
			seq.index <- which.max(sequences$nucleotide_ATGC_length[taxon.index])
			# sequences[taxon.index,][seq.index,]
			seq <- strsplit(sequences$nucleotides[taxon.index][seq.index], split = "")[[1]]
			final.sequences[row.index, sequence(length(seq))] <- seq
		} else {
			taxa.to.drop <- c(taxa.to.drop, i)
		}
		final.sequences.names[row.index] <- i
	}
	rownames(final.sequences) <- gsub(" ", "_", final.sequences.names)
	# final.sequences <- final.sequences[!is.na(final.sequences.names),]
	# taxa.to.drop <- phy$tip.label[which(!phy$tip.label %in% rownames(final.sequences))]
	if(length(input)-length(taxa.to.drop) == 1) {
		if (verbose) cat("BOLD sequences found only for one input name", input[which(!input %in% taxa.to.drop)], "...","\n","\t", "Cannot construct a tree." )
		warning("Not enough sequences available in BOLD. No tree was constructed.")
		# if (usetnrs == FALSE) cat("Setting usetnrs=TRUE might change this, but it is time consuming.", "\n")
		return(NA)
	}
	if(length(taxa.to.drop) > 0) {
		taxa.to.drop.print <- paste(taxa.to.drop, collapse = " | ")
		if (verbose) cat("No", marker, "sequences found for", taxa.to.drop.print, "...", "\n", "\t", "Dropping taxa from tree.", "\n")
		#warning("No ", marker, " sequences found for ", taxa.to.drop.print, "...", "\n", "\t", "Taxa dropped from tree.")
		taxa.to.drop <- gsub(" ", "_", taxa.to.drop)
		phy <- ape::drop.tip(phy, taxa.to.drop)
	}
	if (verbose) cat("Aligning with MAFFT...", "\n")
	alignment <- ape::as.DNAbin(final.sequences)
	alignment <- phangorn::as.phyDat(ips::mafft(alignment))
	if (verbose) cat( "\t", "OK.", "\n", "Estimating BoldOToL tree...", "\n")
	pml.object <- phangorn::pml(phangorn::acctran(phy, alignment), data=alignment)
	phy <- pml.object$tree
	if(!ape::is.binary.tree(pml.object$tree)){
		if (verbose) cat("\t", marker, " sequence data available generates a non-dichotomous tree...", "\n", "\t", "Resolving with multi2di...", "\n")
		pml.object$tree <- ape::multi2di(pml.object$tree)
		phy <- pml.object$tree
	}
	if (verbose) cat("\t", "OK.", "\n")
	if (chronogram) {
		if (verbose) cat("Dating BoldOToL tree with chronoMPL...", "\n")
		pml.object$tree <- ape::chronoMPL(pml.object$tree, se = FALSE, test = FALSE)
		phy <- pml.object$tree
		if (verbose) cat("\t", "OK.", "\n")
	}
	if(any(pml.object$tree$edge.length < 0)) {
		warning("\t", "Negative branch lengths in BOLD chronogram.", "\n")
		if(doML) warning("\t", "\t", "Cannot do ML branch length optimization.", "\n")
	} else {
		if(doML) {
			phy <- phangorn::optim.pml(pml.object, data = alignment, rearrangement = "none", optRooted = TRUE, optQ = TRUE)$tree
		}
	}
	phy$tip.label <- gsub('_', ' ', phy$tip.label)
	if (verbose) cat("Done.", "\n")
	return(phy)
}


#
# ReadDistance <- function(file) {
# 	data <- readLines(file, n=-1)[-1] #read in phylip distance, perhaps from SDM
# 	data <- strsplit(data, " +")
# 	data[sapply(data, length)>0] #trim trailing
# 	final.matrix <- matrix(nrow=length(data), ncol=length(data))
# 	all.names <- rep(NA, length(data))
# 	for (data.index in sequence(length(data))) {
# 		local.line <- data[[data.index]]
# 		all.names[data.index] <- local.line[1]
# 		local.line <- as.numeric(local.line[-1])
#
# 	}
# }

#
# SDM <- function(filtered.results, weights=NULL)) {
#
# 	patristic.array <- BindMatrices(filtered.results)
# 	k <- length(filtered.results)
# 	if(is.null(weights)) {
# 		weights <- rep(1/k, k)
# 	} else {
# 		weights <- weights/sum(weights) #just to make sure total weight is 1
# 	}
#
# 	#Use their appendix and stick it in solve
# 	#a*x=b
# 	#The a matrix has a vector of alpha values, then a_i for matrix 1, a_i for matrix 2..., then
#
# }
#
# #rewrite of code in ape by Andrei Popescu niteloserpopescu@gmail.com to allow for debugging
# apeSDM <- function(...) {
# 	st <- list(...)
# 	k <- length(st)/2
# 	ONEtoK <- seq_len(k)
# 	for (i in ONEtoK) st[[i]] <- as.matrix(st[[i]])
# 	ROWNAMES <- lapply(st[ONEtoK], rownames)
# 	NROWS <- sapply(ROWNAMES, length)
# 	tot <- sum(NROWS)
# 	labels <- unique(unlist(ROWNAMES))
# 	sp <- unlist(st[k + ONEtoK])
# 	astart <- numeric(tot)
# 	astart[1] <- k
# 	for (i in 2:k) astart[i] <- astart[i - 1] + NROWS[i - 1]
# 	n <- length(labels)
# 	miustart <- k + tot
# 	niustart <- miustart + n
# 	lambstart <- niustart + k - 1
# 	X <- matrix(0, n, n, dimnames = list(labels, labels))
# 	V <- w <- X
# 	tmp <- 2 * k + tot + n
# 	col <- numeric(tmp)
# 	for (i in 1:(n - 1)) {
# 		for (j in (i + 1):n) {
# 			for (p in ONEtoK) {
# 				if (is.element(labels[i], ROWNAMES[[p]]) && is.element(labels[j],
# 					ROWNAMES[[p]])) {
# 						w[i, j] <- w[j, i] <- w[i, j] + sp[p]
# 					}
# 				}
# 			}
# 		}
# 		ONEtoN <- seq_len(n)
# 		Q <- matrix(0, tmp, tmp)
# 		for (p in ONEtoK) {
# 			d_p <- st[[p]]
# 			for (l in ONEtoK) {
# 				d <- st[[l]]
# 				sum <- 0
# 				dijp <- -1
# 				if (l == p) {
# 					for (i in ONEtoN) {
# 						for (j in ONEtoN) {
# 							if (i == j)
# 							next
# 							pos <- match(labels[c(i, j)], ROWNAMES[[l]])
# 							if (all(!is.na(pos))) {
# 								ipos <- pos[1L]
# 								jpos <- pos[2L]
# 								dij <- d[ipos, jpos]
# 								sum <- sum + dij * dij - sp[l] * dij *
# 								dij/w[i, j]
# 								tmp2 <- dij - sp[l] * dij/w[i, j]
# 								Q[p, astart[l] + ipos] <- Q[p, astart[l] +
# 								ipos] + tmp2
# 								Q[p, astart[l] + jpos] <- Q[p, astart[l] +
# 								jpos] + tmp2
# 							}
# 						}
# 					}
# 				}
# 				else {
# 					for (i in ONEtoN) {
# 						for (j in ONEtoN) {
# 							if (i == j)
# 							next
# 							pos <- match(labels[c(i, j)], ROWNAMES[[l]])
# 							posp <- match(labels[c(i, j)], ROWNAMES[[p]])
# 							if (all(!is.na(pos)) && all(!is.na(posp))) {
# 								ipos <- pos[1L]
# 								jpos <- pos[2L]
# 								dij <- d[ipos, jpos]
# 								dijp <- d_p[posp[1L], posp[2L]]
# 								sum <- sum - sp[l] * dij * dijp/w[i, j]
# 								tmp2 <- sp[l] * dijp/w[i, j]
# 								Q[p, astart[l] + ipos] <- Q[p, astart[l] +
# 								ipos] - tmp2
# 								Q[p, astart[l] + jpos] <- Q[p, astart[l] +
# 								jpos] - tmp2
# 							}
# 						}
# 					}
# 				}
# 				Q[p, l] <- sum
# 			}
# 			Q[p, lambstart + 1] <- 1
# 		}
# 		r <- k
# 		for (p in ONEtoK) {
# 			dp <- st[[p]]
# 			for (i in ONEtoN) {
# 				if (is.element(labels[i], ROWNAMES[[p]])) {
# 					r <- r + 1
# 					for (l in ONEtoK) {
# 						d <- st[[l]]
# 						if (l == p) {
# 							ipos <- match(labels[i], ROWNAMES[[p]])
# 							for (j in ONEtoN) {
# 								if (i == j)
# 								next
# 								jpos <- match(labels[j], ROWNAMES[[p]])
# 								if (!is.na(jpos)) {
# 									dij <- d[ipos, jpos]
# 									Q[r, l] <- Q[r, l] + dij - sp[l] * dij/w[i,
# 									j]
# 									tmp2 <- 1 - sp[l]/w[i, j]
# 									Q[r, astart[l] + ipos] <- Q[r, astart[l] +
# 									ipos] + tmp2
# 									Q[r, astart[l] + jpos] <- Q[r, astart[l] +
# 									jpos] + tmp2
# 								}
# 							}
# 						}
# 						else {
# 							for (j in ONEtoN) {
# 								if (i == j)
# 								next
# 								if (!is.element(labels[j], rownames(dp)))
# 								next
# 								pos <- match(labels[c(i, j)], ROWNAMES[[l]])
# 								if (all(!is.na(pos))) {
# 									ipos <- pos[1L]
# 									jpos <- pos[2L]
# 									dij <- d[ipos, jpos]
# 									Q[r, l] <- Q[r, l] - sp[l] * dij/w[i,
# 									j]
# 									tmp2 <- sp[l]/w[i, j]
# 									Q[r, astart[l] + ipos] <- Q[r, astart[l] +
# 									ipos] - tmp2
# 									Q[r, astart[l] + jpos] <- Q[r, astart[l] +
# 									jpos] - tmp2
# 								}
# 							}
# 						}
# 					}
# 					if (p < k)
# 					Q[r, ] <- Q[r, ] * sp[p]
# 					Q[r, miustart + i] <- 1
# 					if (p < k)
# 					Q[r, niustart + p] <- 1
# 				}
# 			}
# 		}
# 		r <- r + 1
# 		col[r] <- k
# 		Q[r, ONEtoK] <- 1
# 		for (i in ONEtoN) {
# 			r <- r + 1
# 			for (p in ONEtoK) {
# 				ipos <- match(labels[i], ROWNAMES[[p]])
# 				if (!is.na(ipos))
# 				Q[r, astart[p] + ipos] <- 1
# 			}
# 		}
# 		for (p in 1:(k - 1)) {
# 			r <- r + 1
# 			for (i in ONEtoN) {
# 				ipos <- match(labels[i], ROWNAMES[[p]])
# 				if (!is.na(ipos))
# 				Q[r, astart[p] + ipos] <- 1
# 			}
# 		}
# 		a <- solve(Q, col, 1e-19)
# 		for (i in ONEtoN) {
# 			for (j in ONEtoN) {
# 				if (i == j) {
# 					X[i, j] <- V[i, j] <- 0
# 					next
# 				}
# 				sum <- 0
# 				sumv <- 0
# 				for (p in ONEtoK) {
# 					d <- st[[p]]
# 					pos <- match(labels[c(i, j)], ROWNAMES[[p]])
# 					if (all(!is.na(pos))) {
# 						ipos <- pos[1L]
# 						jpos <- pos[2L]
# 						dij <- d[ipos, jpos]
# 						sum <- sum + sp[p] * (a[p] * dij + a[astart[p] +
# 							ipos] + a[astart[p] + jpos])
# 							sumv <- sumv + sp[p] * (a[p] * dij)^2
# 						}
# 					}
# 					X[i, j] <- sum/w[i, j]
# 					V[i, j] <- sumv/(w[i, j])^2
# 				}
# 	}
# 	list(X, V)
# }
#
#
#
#
#
#' Function to compute the SDM supertree Criscuolo et al. 2006
#' @param filtered.results List of patristic matrices
#' @param weighting flat, taxa, inverse
#' @return A list containing phy (a chronogram), and filtered.results that were actually used
#' @export
#' @details
#' Weighting is how much weight to give each input tree.
#'    flat = all trees have equal weighting
#'    taxa = weight is proportional to number of taxa
#'    inverse = weight is proportional to 1 / number of taxa
#' Criscuolo A, Berry V, Douzery EJ, Gascuel O. SDM: a fast distance-based approach for (super) tree building in phylogenomics. Syst Biol. 2006;55(5):740–55. doi: 10.1080/10635150600969872.
RunSDM <- function(filtered.results, weighting="flat") {
	phy <- NA
	used.studies <- names(filtered.results)
	unpadded.matrices <- lapply(filtered.results, UnpadMatrix)
	good.matrix.indices <- c()
	for(i in sequence(length(unpadded.matrices))) {
		test.result <- NA
		# Rationale here: some chronograms always cause errors with SDM, even when trying to get a consensus of them
		# with themselves. For now, throw out of synthesis.
		try(test.result <- mean(do.call(ape::SDM, c(unpadded.matrices[i], unpadded.matrices[i], rep(1, 2)))[[1]]), silent=TRUE)
		if(is.finite(test.result)) {
			good.matrix.indices <- append(good.matrix.indices,i)
		}
	}
	if(length(good.matrix.indices)>0) {
		unpadded.matrices <- unpadded.matrices[good.matrix.indices]
		used.studies <- used.studies[good.matrix.indices]
		weights = rep(1, length(unpadded.matrices))
		if (weighting=="taxa") {
			weights = unname(sapply(unpadded.matrices, dim)[1,])
		}
		if (weighting=="inverse") {
			weights = 1/unname(sapply(unpadded.matrices, dim)[1,])
		}
		SDM.result <- do.call(ape::SDM, c(unpadded.matrices, weights))[[1]]
		#agnes in package cluster has UPGMA with missing data; might make sense here
		try(phy <- phangorn::upgma(SDM.result))
	}
	return(list(phy=phy, filtered.results=unpadded.matrices))
}
