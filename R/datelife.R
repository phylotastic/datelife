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
#' @param input Target taxa names as a character vector, a newick character string, or a phylo object.
#' @param summary_format The desired output format for target chronograms (chronograms of target taxa). See details.
#' @param summary_print A character vector specifying type of summary information to be printed: "citations" for the references of chronograms from cache where target taxa are found, "taxa" for a summary of the number of chronograms where each target taxon is found, or "none" if nothing should be printed. Default to display both c("citations", "taxa").
#' @param add_taxon_distribution A character vector specifying if data on target taxa missing in source chronograms should be added to the output as a "summary" or as a presence/absence "matrix". Default to "none", no information on add_taxon_distribution added to the output.
#' @param partial If TRUE, use source chronograms even if they only match some of the desired taxa
#' @param use_tnrs If TRUE, use OpenTree's services to resolve names. This can dramatically improve the chance of matches, but also take much longer.
#' @param approximate_match If TRUE, use a slower TNRS to correct mispellings, increasing the chance of matches (including false matches).
#' @param update_cache default to FALSE
#' @param cache The cached set of chronograms and other info from data(opentree_chronograms).
#' @param dating_method The method used for tree dating.
#' @param get_spp_from_taxon boolean vector, default to FALSE. If TRUE, will get all species names from taxon names given in input. Must have same length as input. If input is a newick string , with some clades it will be converted to phylo object phy, and the order of get_spp_from_taxon will match phy$tip.label.
#' @param verbose Boolean. If TRUE, it gives printed updates to the user.
#' @param criterion Whether to get the grove with the most trees or the most taxa
#' @export
#' @details
#' Available output formats are:
#'
#' citations: A character vector of references where chronograms with some or all of the target taxa are published (source chronograms).
#'
#' mrca: A named numeric vector of most recent common ancestor (mrca) ages of target taxa defined in input, obtained from the source chronograms. Names of mrca vector are equal to citations.
#'
#' newick_all: A named character vector of newick strings corresponding to target chronograms derived from source chronograms. Names of newick_all vector are equal to citations.
#'
#' newick_sdm: Only if multiple source chronograms are available. A character vector with a single newick string corresponding to a target chronogram obtained with SDM supertree method (Criscuolo et al. 2006).
#'
#' newick_median: Only if multiple source chronograms are available. A character vector with a single newick string corresponding to a target chronogram from the median of all source chronograms.
#'
#' phylo_sdm: Only if multiple source chronograms are available. A phylo object with a single target chronogram obtained with SDM supertree method (Criscuolo et al. 2006).
#'
#' phylo_median: Only if multiple source chronograms are available. A phylo object with a single target chronogram obtained from source chronograms with median method.
#'
#' phylo_all: A named list of phylo objects corresponding to each target chronogram obtained from available source chronograms. Names of phylo_all list correspond to citations.
#'
#' phylo_biggest: The chronogram with the most taxa. In the case of a tie, the chronogram with clade age closest to the median age of the equally large trees is returned.
#'
#' html: A character vector with an html string that can be saved and then opened in any web browser. It contains a 4 column table with data on target taxa: mrca, number of taxa, citations of source chronogram and newick target chronogram.
#'
#' data_frame A data frame with data on target taxa: mrca, number of taxa, citations of source chronograms and newick string.
#'
#' For approaches that return a single synthetic tree, it is important that the trees leading to it form a grove (roughly, a sufficiently overlapping set of taxa between trees: see Ané et al. 2005, 10.1007/s00026-009-0017-x). In the rare case of multiple groves, should we take the one with the most trees or the most taxa?
#'
#' @examples
#' # obtain median ages from a set of source chronograms in newick format:
#' ages <- datelife_search(c("Rhea americana", "Pterocnemia pennata", "Struthio camelus",
#' 		"Mus musculus"), summary_format="newick_median")
#' # save the tree in newick format
#' write(ages, file="some.bird.ages.txt")
#'
#' # obtain median ages from a set of source chronograms in phylo format
#' # will produce same tree as above but in r phylo format:
#' ages.again <- datelife_search(c("Rhea americana", "Pterocnemia pennata", "Struthio camelus",
#' 		"Mus musculus"), summary_format="phylo_median")
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
		summary_format = "phylo_all", partial = TRUE, use_tnrs = FALSE, approximate_match = TRUE, update_cache = FALSE, cache = get("opentree_chronograms"), dating_method="PATHd8", summary_print= c("citations", "taxa"), add_taxon_distribution = c("none", "summary", "matrix"),  get_spp_from_taxon = FALSE, verbose = FALSE, criterion="taxa") {
			# find a way not to repeat partial and cache arguments, which are used in both get_datelife_result and summarize_datelife_result
			if(update_cache){
				cache <- update_datelife_cache(save = TRUE, verbose = verbose)
			}
			datelife_query <- make_datelife_query(input = input, use_tnrs = use_tnrs, approximate_match = approximate_match, get_spp_from_taxon = get_spp_from_taxon, verbose = verbose)
			datelife_result.here <- get_datelife_result(input = datelife_query, partial = partial, use_tnrs = use_tnrs, approximate_match = approximate_match, update_cache = FALSE, cache = cache, dating_method = dating_method, verbose = verbose)
			# print.datelife(datelife_result = datelife_result.here)
			# datelife <- list(datelife_query = datelife_query, datelife_result = datelife_result.here)
			# class(datelife) <- "datelife"
			# return(datelife)
			return(summarize_datelife_result(datelife_query = datelife_query, datelife_result = datelife_result.here, summary_format = summary_format, partial = partial, update_cache = FALSE, cache = cache, summary_print = summary_print, add_taxon_distribution = add_taxon_distribution, verbose = verbose, criterion=criterion))
}

# print.datelife <- function(datelife){
# 	datelife_result_check(datelife$datelife_result)
# 	cat("Number and which queried taxa found in chronogram database")
# 	cat("Number of queried taxa not found")
# 	cat("Number of chronograms with at least two queried taxa found in database")
# }

#' Go from a vector of species, newick string, or phylo object to a list of patristic matrices
#' @inheritParams datelife_search
#' @inheritParams make_datelife_query
# #' @inheritParams make_bold_otol_tree
# #' @inheritDotParams make_bold_otol_tree
#' @return A datelifeResult object (named list of patristic matrices).
#' @export
get_datelife_result <- function(input = c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"), partial = TRUE, use_tnrs = FALSE, approximate_match = TRUE, update_cache = FALSE, cache = get("opentree_chronograms"), dating_method="PATHd8", get_spp_from_taxon = FALSE, verbose = FALSE) {
	if(update_cache){
		cache <- update_datelife_cache(save = TRUE, verbose = verbose)
	}
	input <- datelife_query_check(datelife_query = input, use_tnrs = use_tnrs, approximate_match = approximate_match, get_spp_from_taxon = get_spp_from_taxon, verbose = verbose)
	# tree <- input$phy
	# cleaned_names <- input$cleaned_names
	datelife_query_length_check(cleaned_names = input$cleaned_names, get_spp_from_taxon = get_spp_from_taxon, verbose = verbose)
  results_list <- lapply(cache$trees, get_subset_array_dispatch, taxa = input$cleaned_names, phy = input$phy, dating_method = dating_method)
  datelife_result <- results_list_process(results_list, input$cleaned_names, partial)
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
#' @param datelife_query An object output of make_datelife_query function
#' @inheritParams datelife_search
#' @inheritDotParams make_datelife_query
#' @export
datelife_query_check <- function(datelife_query = NULL, ...){
	# if(is.null(input)){
	# 	input <- NA
	# }
	# if(length(input) == 1 & any(is.na(input))) {
	# 	stop("input argument is NULL or NA")
	# }
	if(missing(datelife_query) | is.null(datelife_query)){
		stop("datelife_query argument is missing, we have nothing to check")
	}
	badformat <- TRUE
	if(is.list(datelife_query) & "phy" %in% names(datelife_query) & "cleaned_names" %in% names(datelife_query)) {
		if(!inherits(datelife_query, "datelifeQuery")) {
			class(datelife_query) <- "datelifeQuery"
		}
		badformat <- FALSE
	}
	if(badformat){
		datelife_query <- make_datelife_query(input = datelife_query, ...)
		badformat <- FALSE  # useful for next block
	}
	# if(!badformat){
	# 	# merge datelife_query_length_check function here
	# }
	return(datelife_query)
}
#' checks that we have at least two taxon names to perform a search
#' @inheritParams datelife_search
#' @param cleaned_names A character vector; usually an output from make_datelife_query function
#' @export
datelife_query_length_check <- function(cleaned_names = NULL, get_spp_from_taxon = FALSE, verbose = FALSE){
	if(length(cleaned_names) == 1){
		if(verbose) {
			message("Cannot perform a search of divergence times with just one taxon.")
			if(get_spp_from_taxon) {
				message("\t", "Clade contains only one lineage.")
			} else {
				message("Performing a clade search?? set get_spp_from_taxon = TRUE")
			}
		}
		stop("Search is length 1.")
	}
}
#' checks if we obtained an empty search with the set of input taxon names
#' @inheritParams datelife_search
#' @param datelife_result An output of get_datelife_result function: A list of patristic matrices with names corresponding to original study citations.
#' @export
datelife_result_check <- function(datelife_result, use_tnrs, verbose = FALSE){
	if(length(datelife_result) < 1) {
		warning("Output is empty.", call. = FALSE)
		if(verbose) {
			message("Input species were not found in any chronograms available in cache.")
			if(!use_tnrs) {
				message("Setting use_tnrs = TRUE might change this, but it is time consuming.")
			}
		}
	}
}
#' Takes a phylo object or a character string and figure out if it's correct newick format or a list of species
#' @inheritParams datelife_search
#' @inheritParams make_datelife_query
#' @return A phylo object or NA if no tree
#' @export
input_process <- function(input, verbose = FALSE){
	if(inherits(input, "phylo")) {
		input <- ape::write.tree(input)
	}
  	if(length(input)>1) {
		message("There are multiple elements in input. Only the first one will be processed.")
		if(inherits(input, "multiPhylo")) {
			input <- ape::write.tree(input[[1]])
		} else {
			input <- input[1]
		}
	}
 	input <- gsub("\\+"," ",input)
  	input <- stringr::str_trim(input, side = "both")
  	phy.new.in <- NA
   	# if(length(input) == 1) {
    	# if(summary_print)cat("\t", "Input is length 1.", "\n")
	  	if(any(grepl("\\(.*\\).*;", input))) { #our test for newick
	    	phy.new.in <- ape::collapse.singles(phytools::read.newick(text = gsub(" ", "_", input)))
	    	if(verbose) {
					message("Input is a phylogeny and it is correcly formatted.")
				}
	  	} else {
	  		if(verbose) {
					message("Input is not a phylogeny.")
				} #not a warning nor stop, 'cause it is not a requirement for input to be a phylogeny at this step
	  	}
  	# }
	return(phy.new.in)
}

#' Cleans taxon names from input character vector, phylo object or newick character string. Process the two latter with input_process first.
#' @inheritParams datelife_search
#' @inheritDotParams rphylotastic::taxon_get_species -taxon
#' @return A list with the phy (or NA, if no tree) and cleaned vector of taxa
#' @export
make_datelife_query <- function(input = c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"), use_tnrs = FALSE, approximate_match = TRUE, get_spp_from_taxon = FALSE, verbose = FALSE, ...) {
	if(verbose) {
		message("Processing input...")
	}
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
		input <- strsplit(input, ',')[[1]] #  splits a character vector of comma separated names
		if(length(input) == 1){
			if(!get_spp_from_taxon[1]) {
					if(verbose) {
						message("Datelife needs at least two input taxon names to perform a search.")
						message("Setting get_spp_from_taxon = TRUE gets all species from a clade and accepts only one taxon name as input.")
					}
					stop("Input is length 1. If it is a taxon name, try with higher-taxon search. If it is a newick tree, check its format; it is not readable by DateLife.")
			}
		}
	}
	cleaned.input <- stringr::str_trim(input, side = "both")  # cleans the input of lingering unneeded white spaces
	ott_ids <- NA
  if (use_tnrs | any(get_spp_from_taxon)) {
		# process names even if it's a "higher" taxon name:
		cleaned.input_tnrs <- tnrs_match(names = cleaned.input)
		cleaned.input <- cleaned.input_tnrs$unique_name
		ott_ids <- cleaned.input_tnrs$ott_id
		# after some tests, decided to use rotl's method instead of taxize::gnr_resolve, and just output the original input and the actual query for users to check out.
		# cleaned.input <- taxize::gnr_resolve(names = cleaned.input, data_source_ids=179, fields="all")$matched_name
		if(!is.na(phy.new[1])) {
			phy.new$ott_ids <- cleaned.input_tnrs$ott_id
		}
  }
	cleaned_names <- gsub("_", " ", cleaned.input)
    if(any(get_spp_from_taxon)){
    	if(length(get_spp_from_taxon) == 1) {
				get_spp_from_taxon <- rep(get_spp_from_taxon, length(cleaned.input))
			}
    	if(length(cleaned.input)!= length(get_spp_from_taxon)){
    		if(verbose) {
					message("Specify all taxa in input to get species names from.")
				}
    		message("input and get_spp_from_taxon arguments must have same length.")
				return(NA)
    	}
			# using tol_subtree will give subspecies too \o/
			# so we might wanna stick to our own function only then...
			tip_labels <-  lapply(cleaned.input_tnrs$ott_id, function(x)
				res <- tryCatch(suppressWarnings(rotl::tol_subtree(x))$tip.label,
				error = function(e) NA))
			ott_ids <- lapply(tip_labels, function(x) as.numeric(gsub(".*_ott", "", x)))
			cleaned_names <- lapply(tip_labels, function(x) gsub("_ott.*", "", x))

			failures <- which(sapply(tip_labels, function(x)length(x) ==1))
			if(length(failures) >0){
				# cleaned.input_tnrs$unique_name[failures]
				ff <- get_ott_children(ott_id = cleaned.input_tnrs$ott_id[failures], ott_rank = "species")
				get_ott_children(ott_id = cleaned.input_tnrs$ott_id[1], ott_rank = "species")
				cleaned_names[failures] <- lapply(ff, rownames)
				ott_ids[failures] <- sapply(ff, "[", "ott_id")
				# data.frame(sapply(sapply(ff, "[", "rank"), length))
				# ff[11]
				# cleaned.input_tnrs$ott_id[failures][10]
				# length(ff) == length(failures)
			}
			# rphylotastic service includes trees absent from synthetic tree:
			# rphylotastic::taxon_get_species(taxon = cleaned.input[1])

    	# species.names <- vector()
    	# index <- 1
	    # for (i in get_spp_from_taxon){
	    # 	if (i) {
	    # 		# spp <- tryCatch(rphylotastic::taxon_get_species(taxon = cleaned_names[index], ...), error = function (e) NULL)
			# 		spp <- rotl::tol_subtree(gsub("ott", "", cleaned.input_tnrs$ott_ids)[1])$tip.label
	    # 		if(length(spp) == 0) {
	    # 			if(verbose) {
			# 				# message("\t", " No species names found for taxon ", cleaned_names[index], ".")
			# 				if (!use_tnrs) {
			# 					message("\t", "Setting use_tnrs = TRUE might change this, but it can be slow.")
			# 				}
			# 			}
			# 			spp <- cleaned_names[index]
	    # 			message(paste("No species names available for input taxon '", cleaned_names[index], "'", sep=""))
	    # 		}
	    # 		species.names <- c(species.names, spp)
	    # 	} else {
	    # 		species.names <- c(species.names, cleaned_names[index])
	    # 	}
	    # 	index <- index + 1
	    # }
			cleaned_names <- lapply(cleaned_names, function(x) gsub("_", " ", x))
			cleaned_names <- unlist(cleaned_names)
			ott_ids <- unlist(ott_ids)
		}
  if(verbose) {
		message("OK.")
	}
  cleaned_names.print <- paste(cleaned_names, collapse = " | ")
  if(verbose) {
		message("Working with the following taxa: ", "\n", "\t", cleaned_names.print)
	}
	datelife_query.return <- list(cleaned_names = cleaned_names, ott_ids = ott_ids, phy = phy.new)
	class(datelife_query.return) <- "datelifeQuery"
	return(datelife_query.return)
}
