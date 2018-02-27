# skeleton of code to make open access info on dating life
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
#' @param summary_format The desired output format for target chronograms (chronograms of target taxa). See details.
#' @param summary_print A character vector specifying type of summary information to be printed: "citations" for the references of chronograms from cache where target taxa are found, "taxa" for a summary of the number of chronograms where each target taxon is found, or "none" if nothing should be printed. Default to display both c("citations", "taxa").
#' @param add_taxon_distribution A character vector specifying if data on target taxa missing in source chronograms should be added to the output as a "summary" or as a presence/absence "matrix". Default to "none", no information on add_taxon_distribution added to the output.
#' @param partial If TRUE, use source chronograms even if they only match some of the desired taxa
#' @param use_tnrs If TRUE, use OpenTree's services to resolve names. This can dramatically improve the chance of matches, but also take much longer.
#' @param approximate_match If TRUE, use a slower TNRS to correct mispellings, increasing the chance of matches (including false matches).
#' @param update_cache default to FALSE
#' @param cache The cached set of chronograms and other info from data(opentree_chronograms).
#' @param dating_method The method used for tree dating.
# #' @param bold Logical. If TRUE, use Barcode of Life Data Systems (BOLD)  and Open Tree of Life (OToL) backbone to estimate branch lengths of target taxa using make_bold_otol_tree function.
#' @param get_spp_from_taxon boolean vector, default to FALSE. If TRUE, will get all species names from taxon names given in input. Must have same length as input. If input is a newick string , with some clades it will be converted to phylo object phy, and the order of get_spp_from_taxon will match phy$tip.label.
#' @param verbose Boolean. If TRUE, it gives printed updates to the user.
# #' @inheritDotParams make_bold_otol_tree otol_version chronogram doML
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
#' ages <- datelife_search(c("Rhea americana", "Pterocnemia pennata", "Struthio camelus",
#' 		"Mus musculus"), summary_format="newick.median")
#' # save the tree in newick format
#' write(ages, file="some.bird.ages.txt")
#'
#' # obtain median ages from a set of source chronograms in phylo format
#' # will produce same tree as above but in r phylo format:
#' ages.again <- datelife_search(c("Rhea americana", "Pterocnemia pennata", "Struthio camelus",
#' 		"Mus musculus"), summary_format="phylo.median")
#' plot(ages.again)
#' library(ape)
#' ape::axisPhylo()
#' mtext("Time (million years ago)", side = 1, line = 2, at = (max(get("last_plot.phylo",
#' 		envir = .PlotPhyloEnv)$xx) * 0.5))
#' write.tree(ages.again, file="some.bird.tree.again.txt") # saves phylo object in newick format
#'
#' # obtain mrca ages and target chronograms from all source chronograms
#' # generate an html  output readable in any web browser:
#' ages.html <- datelife_search(c("Rhea americana", "Pterocnemia pennata", "Struthio camelus",
#' 		"Mus musculus"), summary_format="html")
#' write(ages.html, file="some.bird.trees.html")
#' system("open some.bird.trees.html")

datelife_search <- function(input = c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"),
		summary_format = "phylo.all", partial = TRUE, use_tnrs = FALSE, approximate_match = TRUE, update_cache = FALSE, cache = get("opentree_chronograms"), dating_method="PATHd8", summary_print= c("citations", "taxa"), add_taxon_distribution = c("none", "summary", "matrix"),  get_spp_from_taxon = FALSE, verbose = FALSE) {
			# find a way not to repeat partial and cache arguments, which are used in both get_datelife_result and summarize_datelife_result
			if(update_cache){
				cache <- update_datelife_cache(save = TRUE, verbose = verbose)
			}
			datelife_query <- make_datelife_query(input = input, use_tnrs = use_tnrs, approximate_match = approximate_match, get_spp_from_taxon = get_spp_from_taxon, verbose = verbose)
			datelife_result.here <- get_datelife_result(input = datelife_query, partial = partial, use_tnrs = use_tnrs, approximate_match = approximate_match, update_cache = FALSE, cache = cache, dating_method = dating_method, verbose = verbose)
			return(summarize_datelife_result(input = datelife_query$cleaned_names, datelife_result = datelife_result.here, summary_format = summary_format, partial = partial, update_cache = FALSE, cache = cache, summary_print = summary_print, add_taxon_distribution = add_taxon_distribution, verbose = verbose))
}

#' Go from a vector of species, newick string, or phylo object to a list of patristic matrices
#' @inheritParams datelife_search
#' @inheritParams make_datelife_query
# #' @inheritParams make_bold_otol_tree
# #' @inheritDotParams make_bold_otol_tree
#' @return List of patristic matrices
#' @export
get_datelife_result <- function(input = c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"), partial = TRUE, use_tnrs = FALSE, approximate_match = TRUE, update_cache = FALSE, cache = get("opentree_chronograms"), dating_method="PATHd8", get_spp_from_taxon = FALSE, verbose = FALSE) {
	if(update_cache){
		cache <- update_datelife_cache(save = TRUE, verbose = verbose)
	}
	input <- input_check(input = input, use_tnrs = use_tnrs, approximate_match = approximate_match, get_spp_from_taxon = get_spp_from_taxon, verbose = verbose)
	tree <- input$phy
	cleaned_names <- input$cleaned_names
	datelife_query_length_check(cleaned_names = cleaned_names, get_spp_from_taxon = get_spp_from_taxon, verbose = verbose)
  results_list <- lapply(cache$trees, get_subset_array_dispatch, taxa = cleaned_names, phy = tree, dating_method = dating_method)
  datelife_result <- results_list_process(results_list, cleaned_names, partial)
	datelife_result_check(datelife_result, use_tnrs)
	# if(bold){
	# 	 bold.OToLTree <- make_bold_otol_tree(input = input, use_tnrs = FALSE, approximate_match = FALSE, marker = marker, verbose = verbose,  ...)
	# 	 bold.data <- phylo_subset_both(reference_tree.in = bold.OToLTree, taxa.in = cleaned_names, phy.in = NULL, phy4.in = NULL, dating_method.in = dating_method)
	# 	 bold.data.processed <- results_list_process(results_list = list(bold.data), taxa = cleaned_names, partial)
	#  	 names(bold.data.processed) <-  paste("BOLD-OToL tree (using ", marker, " as marker)", sep="")
	#    datelife_result <- c(datelife_result, bold.data.processed)
	# }
#	cat("\n")
	class(datelife_result) <- c("datelifeResult")
	return(datelife_result)
}

#' checks if input is a datelifeQuery object, otherwise it uses make_datelife_query to process it
#' @inheritParams datelife_search
#' @inheritDotParams make_datelife_query
#' @export
input_check <- function(input = NULL, ...){
	if(is.null(input)){
		input <- NA
	}
	if(length(input) == 1 & any(is.na(input))) {
		stop("input argument is NULL or NA")
	}
	badformat <- TRUE
	if(is.list(input) & "phy" %in% names(input) & "cleaned_names" %in% names(input)) {
		badformat <- FALSE
	}
	if(badformat){
		input <- make_datelife_query(input = input, ...)
	}
	return(input)
}
#' checks that we have at least two taxon names to perform a search
#' @inheritParams datelife_search
#' @param cleaned_names A character vector; usually an output from make_datelife_query function
#' @export
datelife_query_length_check <- function(cleaned_names = NULL, get_spp_from_taxon = FALSE, verbose = FALSE){
	if(length(cleaned_names) == 1){
		if(verbose) cat("Cannot perform a search of divergence times with just one taxon.", "\n")
		if(get_spp_from_taxon) {
			if(verbose) cat("Clade contains only one lineage.", "\n")
		} else {
			if(verbose) cat("Performing a clade search? set get_spp_from_taxon = TRUE.", "\n")
		}
		stop("search is length 1")
	}
}
#' checks if we obtained an empty search with the set of input taxon names
#' @inheritParams datelife_search
#' @param datelife_result An object output from get_datelife_result function
#' @export
datelife_result_check <- function(datelife_result, use_tnrs, verbose = FALSE){
	if(length(datelife_result) < 1) {
		warning("Output is empty.", call. = FALSE)
		if(verbose) cat("Input species were not found in any chronograms available in cache.", "\n")
		if(!use_tnrs & verbose) cat("Setting use_tnrs = TRUE might change this, but it is time consuming.", "\n")
	}
}
#' Takes a phylo object or a character string and figure out if it's correct newick format or a list of species
#' @inheritParams datelife_search
#' @inheritParams make_datelife_query
#' @return A phylo object or NA if no tree
#' @export
input_process <- function(input, verbose = FALSE){
  if(class(input) == "multiPhylo") stop("Only one phylogeny can be processed at a time.")
	if(class(input) == "phylo") {
		input <- ape::write.tree(input)
	}
 	input <- gsub("\\+"," ",input)
  	input <- stringr::str_trim(input, side = "both")
  	phy.new.in <- NA
   	# if(length(input) == 1) {
    	# if(summary_print)cat("\t", "Input is length 1.", "\n")
	  	if(any(grepl("\\(.*\\).*;", input))) { #our test for newick
	  		if(length(input)>1) {
					stop("Only one phylogeny can be processed at a time.")
				}
	    	phy.new.in <- ape::collapse.singles(phytools::read.newick(text = gsub(" ", "_", input)))
	    	if(verbose) {
					cat("\t", "Input is a phylogeny and it is correcly formatted.", "\n")
				}
	  	} else {
	  		if(verbose) {
					cat("Input is not a phylogeny.")
				} #not a warning nor stop, 'cause it is not a requirement for input to be a phylogeny at this step
	  	}
  	# }
	return(phy.new.in)
}

#' Cleans taxon names from input character vector, phylo object or newick character string. Process the two latter with input_process first.
#' @inheritParams datelife_search
#' @inheritDotParams rphylotastic::taxon_get_species
#' @return A list with the phy (or NA, if no tree) and cleaned vector of taxa
#' @export
make_datelife_query <- function(input = c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"), use_tnrs = FALSE, approximate_match = TRUE, get_spp_from_taxon = FALSE, verbose = FALSE, ...) {
	if(verbose) cat("Processing input...", "\n")
	phy.new <- input_process(input = input, verbose = verbose)
	# cleaned_names <- ""
	if(!is.na(phy.new[1])) {
	    # if(use_tnrs) {
	      	# phy.new$tip.label <- gsub("_", " ", cleaned_names)
	    # }
	  	# cleaned_names <- gsub("_", " ", phy.new$tip.label)
	  	input <- phy.new$tip.label
	}
	if(length(input) == 1) {
		input <- strsplit(input, ',')[[1]]
		if(!get_spp_from_taxon[1]) {
				if(verbose) {
					cat("Datelife needs at least two input taxon names to perform a search.", "\n")
					cat("Setting get_spp_from_taxon = TRUE gets all species from a clade and accepts only one taxon name as input.", "\n")
				}
				stop("Input is length 1 and not in a good newick format.")
		}
	}
	cleaned.input <- stringr::str_trim(input, side = "both")
    if (use_tnrs) {
			# process names even if it's a "higher" taxon name:
			cleaned.input <- rotl::tnrs_match_names(cleaned.input)$unique_name
			# after some tests, decided to use above method instead of taxize::gnr_resolve, and just output the original input and the actual query for users to check out.
			# cleaned.input <- taxize::gnr_resolve(names = cleaned.input, data_source_ids=179, fields="all")$matched_name
    }
	cleaned_names <- gsub("_", " ", cleaned.input)
    if(any(get_spp_from_taxon)){
    	if(length(get_spp_from_taxon)==1) get_spp_from_taxon <- rep(get_spp_from_taxon,length(cleaned.input))
    	if(length(cleaned.input)!=length(get_spp_from_taxon)){
    		if(verbose) cat("Specify all taxa in input to get species names from.", "\n")
    		stop("input and get_spp_from_taxon arguments must have same length.")
    	}
    	species.names <- vector()
    	index <- 1
	    for (i in get_spp_from_taxon){
	    	if (i) {
	    		spp <- rphylotastic::taxon_get_species(taxon = cleaned_names[index], ...)
	    		if(length(spp)==0) {
	    			if(verbose) cat("\t", " No species names found for taxon ", cleaned_names[index], ".", "\n", sep="")
	    			if (!use_tnrs & verbose) cat("\t", "Setting use_tnrs = TRUE might change this, but it can be slow.", "\n")
	    			warning(paste("No species names available for input taxon '", cleaned_names[index], "'", sep=""))
	    		}
	    		species.names <- c(species.names, spp)
	    	} else {
	    		species.names <- c(species.names, cleaned_names[index])
	    	}
	    	index <- index + 1
	    }
		cleaned_names <- gsub("_", " ", species.names)
	}
	cleaned_names <- unique(cleaned_names)
  if(verbose) cat("OK.", "\n")
  cleaned_names.print <- paste(cleaned_names, collapse = " | ")
  if(verbose) cat("Working with the following taxa:", "\n", "\t", cleaned_names.print, "\n")
	datelife_query.return <- list(cleaned_names = cleaned_names, phy = phy.new)
	class(datelife_query.return) <- "datelifeQuery"
	return(datelife_query.return)
}

#' Takes a tree and fixes negative or zero length branches in several ways
#' @param tree A tree either as a newick character string or as a phylo object
#' @param fixing_criterion A character vector specifying the type of branch length to be fixed: "negative" or "zero"
#' @param fixing_method A character vector specifying the method to fix branch lengths: "bladj", "mrbayes" or a number to be assigned to all branches meeting fixing_criterion
#' @return A phylo object with fixed branch lengths
#' @export
tree_fix_brlen <- function(tree = NULL, fixing_criterion = "negative", fixing_method = 0){
	phy <- tree_check(tree = tree)
	fixing_criterion <- match.arg(arg = fixing_criterion, choices = c("negative", "zero"), several.ok = FALSE)
	if(fixing_criterion == "negative"){
		index <- which(phy$edge.length < 0)  # identifies edge numbers with negative edge lengths value
	} else {
		index <- which(phy$edge.length == 0)  # identifies edge numbers with null/zero edge lengths value
	}
	if(!is.numeric(fixing_method)){
		fixing_method <- match.arg(fixing_method, c("bladj", "mrbayes"))
	} else { # chunk for neg or zero br len to zero or any number determined by user
		for (i in index){
			# snode <- pos.phy$edge[i,1]
			# pool  <- pos.phy$edge[seq(nrow(pos.phy$edge))[-i], 1]
			# sisedge <- which(pool==snode) # determines position of sister edge
			# pos.phy$edge.length[sisedge] <- pos.phy$edge.length[sisedge] - pos.phy$edge.length[i]
			# adds neg branch length to sister branch, should add error to both sides???? or only to the daughter branches??
			cnode <- phy$edge[i,2]
			dauedge <- which(phy$edge[,1] == cnode)
			phy$edge.length[dauedge] <- phy$edge.length[dauedge] + phy$edge.length[i] + fixing_method[1]
			phy$edge.length[i] <- fixing_method[1]
			fixed.phy <- phy
		}
	}

	if(any(fixing_method == c("bladj", "mrbayes"))) { #chunk for bladj and mrbayes
		phy <- tree_add_nodelabels(tree = phy)  # all nodes need to be named
		cnode <- phy$edge[index,2]  # we assume that the negative edge length is the one that needs to be changed (but it could be the sister edge that should be shorter)
		tofix <- cnode-length(phy$tip.label)  # so, we take the crown node number of the negative branch lengths
		if(fixing_method == "bladj")
			fixed.phy <- make_bladj_tree(tree = phy, nodenames = phy$node.label[-tofix], nodeages = tree_get_node_data(tree = phy, node_data = "node_age")$node_age[-tofix])
		if(fixing_method == "mrbayes") {
			mrbayes.file <- paste0("phylo", "_brlen_fixed.nexus")  # make "phylo" an argument?
			ncalibration <- tree_get_node_data(tree = phy, node_data = c("descendant_tips_label", "node_age"))
			ncalibration <- lapply(ncalibration, "[", seq(phy$Nnode)[-tofix])
			phy <- tree_add_outgroup(tree = phy, outgroup = "fake_outgroup")
			fixed.phy <- make_mrbayes_tree(constraint = phy, ncalibration = ncalibration, mrbayes_output_file = mrbayes.file)
			fixed.phy <- ape::drop.tip(fixed.phy, "fake_outgroup")
		}
	}
	return(fixed.phy)
}

#' Gets ages, node numbers, node names and descendant tips number and label of all nodes from a dated tree
#' @inheritParams tree_fix_brlen
#' @param node_data A character vector containing one or all from: "node_number", "node_label", "node_age", "descendant_tips_number", "descendant_tips_label"
#' @return A list
#' @export
tree_get_node_data <- function(tree = NULL, node_data = c("node_number", "node_label", "node_age", "descendant_tips_number", "descendant_tips_label")){
	phy <- tree_check(tree = tree)
	nn <- phylo_get_node_numbers(phy)
	res <- vector(mode = "list")
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
		dt_num <- lapply(nn, function(x) tree_node_tips(tree = phy, node = x))# tip numbers stemming from node
		names(dt_num) <- nn
		dt_lab <- lapply(dt_num, function(x) phy$tip.label[x]) # tip labels corresponding to those tips numbers
		names(dt_lab) <- nn
		if("descendant_tips_number" %in% node_data) {
			res <- c(res, list(descendant_tips_number = dt_num))
		}
		if("descendant_tips_label" %in% node_data) {
			res <- c(res, list(descendant_tips_label = dt_lab))
		}
	}
	return(res)
}

#' Adds labels to nodes with no asigned label
#' @inheritParams tree_fix_brlen
#' @param node_prefix Character vector. If length 1, it will be used to name all nodes with no labels, followed by a number which can be the node_number or consecutive, as specified in node_number
#' @param node_index Character vector choosing one of "consecutive" or "node_number" as node label index. It will use consecutive numbers from 1 to total node number in the first case and node numbers in the second case.
#' @return A phylo object
#' @export
tree_add_nodelabels <- function(tree = NULL, node_prefix="n", node_index="node_number"){
	phy <- tree_check(tree = tree)
	node_index <- match.arg(arg = node_index, choices = c("consecutive","node_number"), several.ok = FALSE)
	if("node_number" %in% node_index){
		node_number <- phylo_get_node_numbers(phy = phy)
	}
	if("consecutive" %in% node_index){
		node_number <- seq(phy$Nnode)
	}
	if(is.null(phy$node.label)){
		phy$node.label <- paste0(node_prefix, node_index)
	} else {
		en <- which(phy$node.label == "")
		phy$node.label[en] <- paste0(node_prefix, en)
	}
	return(phy)
}

#' Gets node numbers from any phylogeny
#' @inheritParams phylo_check
#' @return A numeric vector with node numbers
#' @export
phylo_get_node_numbers <- function(phy){
	node_numbers <- (length(phy$tip.label) + 1):(length(phy$tip.label) + phy$Nnode)
	return(node_numbers)
}

#' Function to add an outgroup to any phylogeny, in phylo or newick format
#' @inheritParams tree_fix_brlen
#' @param outgroup A character vector with the name of the outgroup. If it has length>1, only first element will be used.
#' @return A phylo object.
#' @export
tree_add_outgroup <- function(tree = NULL, outgroup = "outgroup"){
		phy <- tree_check(tree = tree)
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
	phy <-  phytools::read.newick(text = phy)
	return(phy)
}

#' Checks if phy is a phylo class object and/or a chronogram
#' @param phy A phylo object
#' @param dated Boolean. If TRUE it checks if phylo object has branch lengths and is ultrametric.
#' @return A phylo object
#' @export
phylo_check <- function(phy = NULL, dated = FALSE){
	if (!inherits(phy, "phylo")){
		stop("tree is not a phylo object")
	}
	if(dated){
		if(!phylo_has_brlen(phy = phy)){
			stop("tree must have branch lengths")
		}
		if(!ape::is.ultrametric(phy, option = 2)){
			warning("branch lengths in tree should be proportional to time")
			stop("tree must be ultrametric")  # is this true?  # Think how to incorporate trees with extinct taxa
		}
	}
	# return(phy)
}

#' Checks if a tree is a phylo class object otherwise it uses input_process. Additionally it can check if tree is a chronogram with phylo_check
#' @inheritParams tree_fix_brlen
#' @inheritDotParams phylo_check
#' @return If tree is correctly formatted, it returns a phylo object
#' @export
tree_check <- function(tree = NULL, ...){
	if (!inherits(tree, "phylo")){
		tree <- input_process(input = tree, verbose = FALSE)
	}
	phylo_check(phy = tree, ...)
	return(tree)
}

#' To get tip numbers descending from any given node of a tree
#' @inheritParams phytools::getDescendants
#' @inheritParams tree_fix_brlen
#' @return A numeric vector with tip numbers descending from a node
#' @export
tree_node_tips <- function(tree = NULL, node = NULL, curr = NULL){
	phy <- tree_check(tree = tree, dated = FALSE)
	des <- phytools::getDescendants(tree = phy, node = node, curr = NULL)
	tips <- des[which(des <= length(phy$tip.label))]
	return(tips)
}

#' Takes a tree and uses bladj to estimate node ages and branch lengths given a set of fixed node ages and respective node names
#' @param nodenames A character vector with node names from tree with fixed ages
#' @param nodeages A numeric vector with known or fixed node ages from tree
#' @inheritParams tree_fix_brlen
#' @return A phylo tree
#' @export
make_bladj_tree <- function(tree = NULL, nodenames = NULL, nodeages = NULL){
	phy <- tree_check(tree = tree, dated = FALSE)
	if(is.null(phy$node.label)) {
		stop("phy must have node labels")
	}
	if(!is.null(phy$edge.length)) {
		phy$edge.length <- NULL
	}
	m <- match(nodenames, phy$node.label)
	if(any(is.na(m))) {
		stop("all nodenames must be in phy$node.label") # add a printed line saying which nodenames are not in phy$node.label
	}
	if(length(nodenames) != length(nodeages)) {
		stop("nodenames and nodeages must have the same length")
	}
	if(!is.character(nodenames)) {
		stop("nodenames must be a character vector")
	}
	if(!is.numeric(nodeages)) {
		stop("nodeages must be a numeric vector")
	}
	ages_df <- data.frame(
		a = nodenames,
		b = nodeages
	)
	new.phy <- phylocomr::ph_bladj(ages = ages_df, phylo = phy)
	attributes(new.phy) <- NULL
	new.phy <- ape::read.tree(text = new.phy)
	# plot(new.phy)
	return(new.phy)
}

#' Takes a constraint tree and uses mrBayes to get node ages and branch lengths given a set of node calibrations without any data.
# we can add the option to use data and no constraint tree.
#' @param constraint The constraint tree: a phylo object or a newick character string, with or without branch lengths.
#' @param ncalibration The node calibrations: a phylo object with branch lengths proportional to time; in this case all nodes from ncalibration will be used as calibration points. Alternatively, a list with two elements: the first is a character vector with node names from phy to calibrate; the second is a numeric vector with the corresponding ages to use as calibrations.
#' @inheritParams missing_taxa_check
#' @param mrbayes_output_file A character vector specifying the name of mrBayes run file and outputs (can specify directory too).
#' @return A phylo tree with branch lengths proportional to time. It will save all mrBayes outputs in the working directory.
#' @export
make_mrbayes_tree <- function(constraint = NULL, ncalibration = NULL, missing_taxa = NULL, mrbayes_output_file = "mrbayes_run.nexus"){
	make_mrbayes_runfile(constraint = constraint, ncalibration = ncalibration, missing_taxa = missing_taxa, mrbayes_output_file = mrbayes_output_file)
	message("Starting MrBayes run. This might take a while...")
	new.tree <- run_mrbayes(mrbayes_output_file = mrbayes_output_file)
	return(new.tree)
}

#' Makes a mrBayes run block file with a constraint topology and a set of node calibrations
#' @inheritParams make_mrbayes_tree
#' @return A MrBayes block run file in nexus format.
#' @export
make_mrbayes_runfile <- function(constraint = NULL, ncalibration = NULL, missing_taxa = NULL, mrbayes_output_file = "mrbayes_run.nexus"){
  if(!is.null(constraint)) {
		constraint <- tree_check(tree = constraint) # add outgroup = TRUE argument
  	constraint <- phylo_tiplabel_space_to_underscore(constraint)
		taxa <- constraint$tip.label
		constraints <- paleotree::createMrBayesConstraints(tree = constraint, partial = FALSE) # this works perfectly
		calibrations <- get_mrbayes_node_calibrations(constraint = constraint, ncalibration = ncalibration, ncalibration_type = "fixed")
		og <- tree_get_singleton_outgroup(tree = constraint) #if(outgroup)
	}
	ogroup <- c()
	if(!is.na(og)) {
		ogroup <- paste0("outgroup ", og, ";")
	}
	missing_taxa <- missing_taxa_check(missing_taxa)
	if(!is.null(missing_taxa)){
		if(is.vector(missing_taxa)){
			taxa <- c(taxa, missing_taxa)
		}
	}
	bayes_data <- c(paste("   Begin DATA; \nDimensions ntax=", length(taxa), "nchar = 1;"),
	"Format datatype = DNA gap=- missing=?;",
	"Matrix\n",
	paste(taxa, "?"),
	";")

    bayes_set <- c("   Begin MRBAYES;",
    	"unlink shape=(all) tratio=(all) statefreq=(all) revmat=(all) pinvar=(all);\n",
    	constraints[-length(constraints)],
    	ogroup, "",
    	constraints[length(constraints)], "",
    	"prset nodeagepr = calibrated;", "",
    	calibrations, "\n",
    	"   set usebeagle = no Beaglesse = no;", "",
    	paste("prset ", c("brlenspr = clock:birthdeath", "Extinctionpr = Fixed(0)",
    	"Speciationpr = exponential(1)", "clockvarpr = ibr", "ibrvarpr = exponential(10)"), ";", sep=""),
    	"mcmcp nruns = 1 nchains = 1 ngen = 50000000 samplefreq = 1000;",
    	"mcmc;", "",
    	paste0("sumt filename=", mrbayes_output_file, " burnin = 5000000 contype = halfcompat;\n"),
    	"end;"
    	)

	all <- c(bayes_data, "\n", bayes_set)
	write(all, mrbayes_output_file)
	return(all)
}

#' Fabricates dates of missing taxa (with no data) on an already dated tree.
#' @param dated_phy a tree (newick or phylo) with branch lengths proportional to absolute time
#' @inheritParams make_mrbayes_tree
#' @inheritParams datelife_search
#' @inheritParams missing_taxa_check
#' @return A phylo object
#' @export
tree_add_dates <- function(dated_phy = NULL, missing_taxa = NULL, dating_method = "mrbayes", mrbayes_output_file = "mrbayes_tree_add_dates.nexus"){
	dated_phy <- tree_check(tree = dated_phy, dated = TRUE)
	dating_method <- match.arg(dating_method, c("bladj", "mrbayes"))
	if(dating_method == "bladj"){
		dated_phy <- tree_add_nodelabels(tree = dated_phy)  # all nodes need to be named
		new.phy <- make_bladj_tree(tree = missing_taxa, nodenames = dated_phy$node.label, nodeages = tree_get_node_data(tree = dated_phy, node_data = "node_age")$node_age)
	}
	if(dating_method=="mrbayes"){
		dated_phy <- tree_add_outgroup(tree = dated_phy, outgroup = "an_outgroup")
		ncalibration <- tree_get_node_data(tree = dated_phy, node_data = c("node_age", "descendant_tips_label"))
		new.phy <- make_mrbayes_tree(constraint = dated_phy, ncalibration = ncalibration, missing_taxa = missing_taxa, mrbayes_output_file = mrbayes_output_file)
		new.phy <- ape::drop.tip(new.phy, "an_outgroup")
	}
	return(new.phy)
}

#' Checks that missing_taxa argument is ok to be used by make_mrbayes_runfile function
#' @param  missing_taxa either a tree (as phylo object or as a newick character string, with or without branch lengths) with all the taxa you want at the end, or a data frame assigning missing taxa to nodes in dated tree.
#' @inheritParams tree_add_dates
# # ' @param missing_taxa A phylo object, a newick character string or a dataframe with taxonomic assignations
#' @return A phylo object, a newick character string or a dataframe with taxonomic assignations
#' @export
missing_taxa_check <- function(missing_taxa = NULL, dated_phy = NULL){
	if(is.data.frame(missing_taxa)){ # or is.matrix??
		stop("not implemented yet")
		# checkPastisData
		return(missing_taxa)
	}
	if(is.vector(missing_taxa)){
		missing_taxa <- as.character(missing_taxa)
		return(missing_taxa)
	}
	if(is.null(missing_taxa)) {
		return(NULL)
	}
	missing_taxa <- input_process(missing_taxa)
	if(inherits(missing_taxa, "phylo")){
		phylo_check(phy = dated_phy, dated = TRUE)
		dtINmt <- dated_phy$tip.labels %in% missing_taxa$tip.labels
		mtINdt <- missing_taxa$tip.labels %in% dated_phy$tip.labels
		if (!all(dtINmt)) {
			stop("all taxa in dated_phy must be in missing_taxa tree too")
		}
		missing_taxa_pruned <- ape::drop.tip(missing_taxa, missing_taxa$tip.labels[mtINdt])
		# phylo_prune_missing_taxa(phy = , taxa = ) # use this one
		# dated_phy == missing_taxa_pruned # check that both trees are equal
	} else {
		stop("missing taxa must be NULL; a vector with species names;
		a dataframe with taxonomic assignations; a newick character string; or, a phylo object")
	}
	# IMPORTANT: Add a check that taxa in dated.trees is in reference_tree and viceversa
	return(missing_taxa)
	# if (is.null(reference_tree)){
		# if (is.null(add_taxon_distribution)) {
			# cat("specify a reference_tree or add_taxon_distribution to be added to the dated.trees")
			# stop("")
		# } else {
			# # construct a tree with phylotastic or take that from otol?
			# cat("Constructing a reference_tree with taxa from dated.tree and add_taxon_distribution")
		# }
	# } else {
		# if(!is.null(add_taxon_distribution)){
			# cat("A reference_tree was given, add_taxon_distribution argument is ignored")
		# }
}

#' Runs MrBayes from R.
#' @inheritParams make_mrbayes_tree
#' @return MrBayes outputs.
#' @export
run_mrbayes <- function(mrbayes_output_file = NULL){
	file <- mrbayes_output_file
	# code borrowed from phyloch::mrbayes()
	if(is.null(file)) {
		stop("You must provide a block file for MrBayes run")
	}
	if (.Platform$OS.type == "unix"){
		system(paste("mb > execute", file))
	} else {
		system(paste("mrbayes ", file, ".bayes", sep = ""))
	}
	tr <- ape::read.nexus(paste(file, ".con.tre", sep = ""))
	return(tr)
}

#' Makes a node calibrations block for a MrBayes run file from a list of taxa and ages or from a dated tree. It can follow a constraint tree.
#' @inheritParams make_mrbayes_tree
#' @inheritParams tree_fix_brlen
#' @param mrbayes_calibration_file NULL or a character vector indicating the name of mrbayes calibration block file.
#' @param ncalibration_type A character string specifying the type of calibration. Only "fixed" is implemented for now.
#' @return A set of MrBayes calibration commands printed in console as character strings or as a text file with name specified in file.
#' @export
# This function is set to match node names with constraints obtained from paleotree::GetMrBayesConstraints
get_mrbayes_node_calibrations <- function(constraint = NULL, ncalibration = NULL, ncalibration_type = "fixed", 	mrbayes_calibration_file = NULL){
	if(!is.null(constraint)){
		phy <- tree_check(tree = constraint)
	}
	if(length(ncalibration) == 2){  # if it is a list of descendant tips labels and node ages, from tree_get_node_data function
		if(!is.list(ncalibration)) {
			stop("ncalibration must be a newick character string, in phylo object or a list with taxon names and dates")
		}
		includes.ncalibration <- lapply(ncalibration$descendant_tips_label, function(x) gsub(" ", "_", x))
		nages <- ncalibration$node_age

	} else {  # if it is a tree
			ncalibration <- tree_check(tree = ncalibration, dated = TRUE)
	    phy <- phylo_tiplabel_space_to_underscore(phy)
	    ncalibration <- phylo_tiplabel_space_to_underscore(ncalibration)
	    splits.ncalibration <- ape::prop.part(ncalibration)
	    includes.ncalibration <- lapply(splits.ncalibration, function(x) ncalibration$tip.label[x])
			nages <- ape::branching.times(ncalibration)
	}
	nodes <- sapply(includes.ncalibration, function(tax)
				phytools::findMRCA(phy, tax, type="node")) - length(phy$tip.label)
	calibrations <- paste0("calibrate node", nodes, " = ", ncalibration_type, "(", nages, ");")
	root <- which(nodes == 1)  # tests for the presence of a root calibration, which should be implemented with treeagepr and not with calibrate
	if(length(root) != 0){
		nodes <- nodes[-root]
		nages <- nages[-root]
		calibrations <- c(calibrations, paste0("prset treeagepr = ", ncalibration_type, "(",
		nages[root], ");"))
	}
    if (!is.null(mrbayes_calibration_file)) {
        write(calibrations, mrbayes_calibration_file)
    }
    else {
        return(calibrations)
    }
}

#' Identifies the presence of a single lineage outgroup in a phylogeny.
#' @inheritParams make_mrbayes_tree
#' @inheritParams tree_fix_brlen
#' @return A character vector with the name of the single lineage outgroup. Returns NA if there is none.
#' @export
tree_get_singleton_outgroup <- function(tree = NULL){
  phy <- tree_check(tree = tree)
	phy <- phylo_tiplabel_space_to_underscore(phy)
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


#' Summarize a filtered results list from get_datelife_result function in various ways
#' @param datelife_result A list of patristic matrices; labels correspond to citations
#' @inheritParams datelife_search
#' @inherit datelife_search return details
#' @export
summarize_datelife_result <- function(datelife_result = NULL, summary_format = "phylo.all", input = NULL, partial = TRUE, update_cache = FALSE, cache = get("opentree_chronograms"), summary_print = c("citations", "taxa"), add_taxon_distribution = c("none", "summary", "matrix"), verbose = FALSE) {
		# if(!partial) {
		# 	datelife_result <- datelife_result[which(!sapply(datelife_result, anyNA))]
		# } # not necessary cause already filtered in get_datelife_result
	if(update_cache){
		cache <- update_datelife_cache(save = TRUE, verbose = verbose)
	}
	if(is.null(datelife_result) | !is.list(datelife_result)){
		stop("datelife_result argument must be a list from get_datelife_result function.")
	}
	summary_format.in <- match.arg(summary_format, choices = c("citations", "mrca", "newick.all", "newick.sdm", "newick.median", "phylo.sdm", "phylo.median", "phylo.median", "phylo.all", "html", "data.frame"))
	add_taxon_distribution.in <- match.arg(add_taxon_distribution, choices = c("none", "summary", "matrix"))
	summary_print.in <- match.arg(summary_print, c("citations", "taxa", "none"), several.ok = TRUE)
	if(is.null(input)){
		input.in <- unique(rapply(datelife_result, rownames))
		if(add_taxon_distribution.in !="none") warning("showing taxon summary from taxa found in at least one chronogram, this excludes input taxa not found in any chronogram")
	} else {
		if(!is.character(input)) stop("input must be a character vector")
		input.in <- input
	}
	results.index <- datelife_result_study_index(datelife_result, cache)
	return.object <- NA
	input.match <- unique(rapply(datelife_result, rownames))
	# if(any(!input.match %in% input)) warning("input does not contain all or any taxa from filteredresults object")
	absent.input <- input.in[!input.in %in% input.match]
	if(length(absent.input) <= 0) {
		if(is.null(input)){
			absent.input <- "NULL"
		} else {
			absent.input <- "None"
		}
	}
	if(add_taxon_distribution.in == "matrix"){
		taxon_distribution_list <- vector(mode = "list")
		# tax <- unique(rapply(datelife_result, rownames)) #rownames(datelife_result[[1]])
		for(result.index in sequence(length(datelife_result))){
			n <- rownames(datelife_result[[result.index]])
			m <- match(input.match,n)
			taxon_distribution_list[[result.index]] <- n[m]
		}
		taxon_distribution_matrix <- do.call(rbind, taxon_distribution_list) #transforms a list of names into a matrix of names
		taxon_distribution_matrix <- !is.na(taxon_distribution_matrix) # makes a boolean matrix
		colnames(taxon_distribution_matrix) <- input.match
		rownames(taxon_distribution_matrix) <- sequence(nrow(taxon_distribution_matrix))
	}
	if(add_taxon_distribution.in == "summary" | any(grepl("taxa", summary_print.in))){ # may add here another condition: | makeup_brlen ==TRUE
		# tax <- unique(rapply(datelife_result, rownames)) #rownames(datelife_result[[1]])
		x <- rapply(datelife_result, rownames)
		prop <- c()
		for (taxon in input.match){
			prop <- c(prop, paste(length(which(taxon == x)), "/", length(datelife_result), sep=""))
		}
		taxon_distribution_summary <- data.frame(Taxon = input.match, Chronograms = prop)
	}
	if(summary_format.in == "citations") {
		return.object <- names(datelife_result)
	}
	if(summary_format.in == "mrca") {
		return.object <- datelife_result_MRCA(datelife_result, partial = partial)
	}
	if(summary_format.in == "newick.all") {
		trees <- sapply(datelife_result, patristic_matrix_to_newick)
		return.object <- trees[which(!is.na(trees))]
	}
	if(summary_format.in == "newick.sdm") {
		local.results <- datelife_result_sdm(datelife_result)
		datelife_result <- local.results$datelife_result
		tree <- local.results$phy
		return.object <- ape::write.tree(tree)
	}
	if(summary_format.in == "phylo.sdm") {
		local.results <- datelife_result_sdm(datelife_result)
		datelife_result <- local.results$datelife_result
		tree <- local.results$phy
		return.object <- tree
	}
	if(summary_format.in == "newick.median") {
		patristic.array <- patristic_matrix_list_to_array(datelife_result)
		median.matrix <- summary_patristic_matrix_array(patristic.array)
		tree <- patristic_matrix_to_newick(median.matrix)
		return.object <- tree
	}
	if(summary_format.in == "phylo.median") {
		patristic.array <- patristic_matrix_list_to_array(datelife_result)
		median.matrix <- summary_patristic_matrix_array(patristic.array)
		tree <- patristic_matrix_to_phylo(median.matrix)
		return.object <- tree
	}
	if(summary_format.in == "phylo.all") {
		trees <- lapply(datelife_result, patristic_matrix_to_phylo)
		return.object <- trees[which(!is.na(trees))]
	}
	if(summary_format.in == "html") {
		out.vector1 <- "<table border='1'><tr><th>MRCA Age (MY)</th><th>Ntax</th><th>Citation</th><th>Newick"
		if(add_taxon_distribution.in == "matrix"){
			out.vector1 <- paste(out.vector1, paste("</th><th>", colnames(taxon_distribution_matrix), sep="", collapse=""), sep="")
		}
		out.vector1 <- paste(out.vector1, "</th></tr>", sep="")
		ages <- datelife_result_MRCA(datelife_result, partial = partial)
		trees <- sapply(datelife_result, patristic_matrix_to_newick)
		out.vector2 <- c()
		for(result.index in sequence(length(datelife_result))) {
			out.vector2 <- paste(out.vector2, "<tr><td>",ages[result.index],"</td><td>",sum(!is.na(diag(datelife_result[[result.index]]))), "</td><td>", names(datelife_result)[result.index], "</td><td>", trees[result.index],  sep="")
			if(add_taxon_distribution.in == "matrix"){
				out.vector2 <- paste(out.vector2, paste("</td><td>", taxon_distribution_matrix[result.index,], sep="", collapse=""), sep="")
			}
			out.vector2 <- paste(out.vector2, "</td></tr>", sep="")
		}
		out.vector <- paste(out.vector1, out.vector2, "</table>", sep="")
		if(add_taxon_distribution.in == "summary"){
			taxon_distribution_summary.html <- as.matrix(taxon_distribution_summary)
			out.vector3 <- "<p></p><table border='1'><tr><th>Taxon</th><th>Chronograms</th><tr>"
			for (summary.index in sequence(nrow(taxon_distribution_summary.html))){
				out.vector3 <- paste(out.vector3, paste("</td><td>", taxon_distribution_summary.html[summary.index,], sep="", collapse=""), "</td></tr>", sep="")
			}
			out.vector <- paste(out.vector, out.vector3, "</table>", sep="")
		}
		# out.vector4 <- c()
		if(add_taxon_distribution.in != "none") {
			out.vector4 <- "<p></p><table border='1'><tr><th> </th><th>Absent Taxa</th><tr>"
			for (i in 1:length(absent.input)){
				out.vector4 <- paste(out.vector4, "<tr><td>", i, "</td><td>", absent.input[i], "</td><tr>", sep="")
			}
			out.vector4 <- paste(out.vector4, "</table>", sep="")
			out.vector <- paste(out.vector, out.vector4, sep="")
		}
		return.object <- out.vector
	}
	if(summary_format.in == "data.frame") {
		out.df <- data.frame()
		ages <- datelife_result_MRCA(datelife_result, partial = partial)
		trees <- sapply(datelife_result, patristic_matrix_to_newick)
		for(result.index in sequence(length(datelife_result))) {
			out.line<- data.frame(Age = ages[result.index],Ntax = sum(!is.na(diag(datelife_result[[result.index]]))), Citation = names(datelife_result)[result.index], Newick= trees[result.index])
			if(result.index == 1) {
				out.df <- out.line
			} else {
				out.df <- rbind(out.df, out.line)
			}
		}
		if(add_taxon_distribution.in == "matrix"){
			out.df <- cbind(out.df, taxon_distribution_matrix)
		}
		rownames(out.df) <- NULL
		return.object <- out.df
	}
	if(add_taxon_distribution.in != "none" & summary_format.in !="html"){
		return.object <- list(return.object)
		if(add_taxon_distribution.in == "matrix") {
			if(summary_format.in !="data.frame") return.object <- c(return.object, list(taxon_distribution = taxon_distribution_matrix))
		}
		if(add_taxon_distribution.in == "summary") {
			return.object <- c(return.object, list(taxon_distribution = taxon_distribution_summary))
		}
		return.object <- c(return.object, list(absent.taxa = data.frame(Taxon = absent.input)))
		names(return.object)[1] <- summary_format.in
		if(summary_format.in =="data.frame"){
			names(return.object)[1] <- "results"
			if(add_taxon_distribution.in == "matrix") names(return.object)[1] <- "results_and_missing_taxa"
		}
	}

	if(any("citations" %in% summary_print.in) & !any(summary_format.in %in% c("citations", "html", "data.frame"))) {
		if(summary_format.in == "citations"){
			cat("Target taxa found in trees from:", "\n")
			print(names(datelife_result), quote = FALSE)
			cat("\n")
		} else {
			cat("Source chronograms from:", "\n")
			print(names(datelife_result), quote = FALSE)
			cat("\n")
		}
	}
	if(any(grepl("taxa", summary_print.in)) & add_taxon_distribution.in!="summary") {
		cat("Target taxa presence in source chronograms:", "\n")
		print(taxon_distribution_summary)
		cat("\n")
		cat("Target taxa completely absent from source chronograms:", "\n")
		print(data.frame(Taxon = absent.input))
		cat("\n")
	}
	return(return.object)
}
#' Function to compute the SDM supertree (Criscuolo et al. 2006) from a datelifeResult object.
#' @param datelife_result List of patristic matrices.
#' @param weighting A character vector. One of "flat", "taxa" or "inverse".
#' @return A list containing phy, a chronogram from SDM, and the original datelife_result that were actually used.
#' @export
#' @details
#' Weighting is how much weight to give each input tree.
#' flat = all trees have equal weighting
#' taxa = weight is proportional to number of taxa
#' inverse = weight is proportional to 1 / number of taxa
#' Criscuolo A, Berry V, Douzery EJ, Gascuel O. SDM: a fast distance-based approach for (super) tree building in phylogenomics. Syst Biol. 2006. 55(5):740–55. doi: 10.1080/10635150600969872.
datelife_result_sdm <- function(datelife_result, weighting="flat") {
	phy <- NA
	used.studies <- names(datelife_result)
	unpadded.matrices <- lapply(datelife_result, patristic_matrix_unpad)
	good.matrix.indices <- c()
	for(i in sequence(length(unpadded.matrices))) {
		test.result <- NA
		# Rationale here: some chronograms always cause errors with SDM, even when trying to get a consensus of them
		# with themselves. For now, throw out of synthesis.
		try(test.result <- mean(do.call(ape::SDM, c(unpadded.matrices[i], unpadded.matrices[i], rep(1, 2)))[[1]]), silent = TRUE)
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
	class(unpadded.matrices) <- "datelifeResult"
	return(list(phy = phy, datelife_result = unpadded.matrices))
}
