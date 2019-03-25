
#' Cleans taxon names from input character vector, phylo object or newick character string. Process the two latter with input_process first.
#' @inheritParams datelife_search
#' @inheritDotParams rphylotastic::taxon_get_species -taxon
#' @return A list with the phy (or NA, if no tree) and cleaned vector of taxa
#' @details If input has length 1, get_spp_from_taxon is always set to TRUE from datelife_search, not in here bc other things depend on this function.
#' @export
make_datelife_query <- function(input = c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"),
use_tnrs = FALSE, approximate_match = TRUE, get_spp_from_taxon = FALSE, verbose = FALSE, ...) {
	if(verbose) {
		message("Processing input...")
	}
	# input <- rphylotastic::taxa_get_otol_tree(url_get_scientific_names(
	# 	URL="https://www.nps.gov/yell/learn/nature/upload/BirdChecklist2014.pdf"))
	# input <- threebirds_median
	# input$tip.label[3] <- "ttttttt"
	phy.new <- input_process(input = input, verbose = verbose)
	if(inherits(phy.new, "phylo")) {
	    # if(use_tnrs) {
	      	# phy.new$tip.label <- gsub("_", " ", cleaned_names)
	    # }
	  	# cleaned_names <- gsub("_", " ", phy.new$tip.label)
	  	cleaned.input <- phy.new$tip.label
	} else {
		if(!is.character(input)){
			message("Input must be a character vector, a newick character string, or a phylo object")
			return(NA)
		}
		if(length(input) == 1) { # if it is a character vector of length 1 with comma separated names
			cleaned.input <- strsplit(input, ',')[[1]] #  it splits it by the comma
		}
		cleaned.input <- stringr::str_trim(input, side = "both") # cleans the input of lingering unneeded white spaces
	}
	ott_ids <- NA
	use_tnrs_global <- FALSE
	if(use_tnrs | any(get_spp_from_taxon)){
		use_tnrs_global <- TRUE
	}
	if(inherits(phy.new, "phylo")){
		if(is.null(phy.new$ott_ids)){
			use_tnrs_global <- TRUE
		}
	}
  	if (use_tnrs_global) {
		# process names even if it's a "higher" taxon name:
		cleaned.input_tnrs <- clean_tnrs(tnrs_match(names = cleaned.input),
			remove_nonmatches = TRUE)
		# recover original names of invalid taxa and unmatched:
		ii <- !tolower(cleaned.input)%in%cleaned.input_tnrs$search_string
		cleaned.input <- c(cleaned.input_tnrs$unique_name, cleaned.input[ii])
		if(inherits(phy.new, "phylo")){
			if(is.null(phy.new$ott_ids)){
				cleaned.input <- gsub(" ", "_", cleaned.input)
				ii <- match(cleaned.input_tnrs$search_string, tolower(phy.new$tip.label))
				# after some tests, decided to use rotl's method instead of taxize::gnr_resolve, and just output the original input and the actual query for users to check out.
				# cleaned.input <- taxize::gnr_resolve(names = cleaned.input, data_source_ids=179, fields="all")$matched_name
				# rename the tip labels with tnrs matched names
				phy.new$tip.label[ii] <- cleaned.input[ii]
				ott_ids <- rep(NA, length(cleaned.input))
				ott_ids[ii] <- cleaned.input_tnrs$ott_id
				phy.new$ott_ids <- ott_ids
			}
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
		# rotl::tol_subtree is very fast but returns subspecies too \o/
		# it has no argument to restrict it to species only
		# so we are using our own function that wraps up theire services nicely
		# example: df <- get_ott_children(ott_ids = 698424, ott_rank = "species")
		df <- get_ott_children(ott_ids = cleaned.input_tnrs$ott_id, ott_rank = "species")
		# head(rownames(df[[1]])[grepl("species", df[[1]]$rank)])
		# the following does not work; it gives subspecies back
		# fixing it from get_ott_children function and here too
		cleaned_names <- lapply(df, function (x) rownames(x)[grepl("\\bspecies\\b", x$rank)])
		ott_ids <- lapply(df, function (x) x$ott_id[grepl("\\bspecies\\b", x$rank)])
		cleaned_names <- lapply(cleaned_names, function(x) gsub("_", " ", x))
		cleaned_names <- unlist(cleaned_names)
		ott_ids <- unlist(ott_ids)
		# enhance: create vector original_taxon with original names: rep(cleaned.input[i], length(cleaned_names[i]))
		original_taxa <- lapply(seq(nrow(cleaned.input_tnrs)), function(i)
			rep(cleaned.input_tnrs$unique_name[i], length(cleaned_names[[i]])))
	}
    if(verbose) {
		message("OK.")
	}
    cleaned_names.print <- paste(cleaned_names, collapse = " | ")
    if(verbose) {
		message("Working with the following taxa: ", "\n", "\t", cleaned_names.print)
	}
	# enhance: add original_taxa vector (from get_spp_from_taxon) to output
	datelife_query.return <- list(cleaned_names = cleaned_names, ott_ids = ott_ids,
		phy = phy.new)
	class(datelife_query.return) <- "datelifeQuery"
	return(datelife_query.return)
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
	if(inherits(input, "multiPhylo")) {
		input <- ape::write.tree(input[[1]])
		message("Input is a multiPhylo object. Only the first tree will be used.")
	}
	if(!is.character(input)){
		message("Input must be a character vector of names or a phylo object.")
		return(NA)
	}
 	input <- gsub("\\+"," ",input)
  	input <- stringr::str_trim(input, side = "both")
  	phy.new.in <- NA
	if(any(grepl("\\(.*\\).*;", input))) { #our test for newick
		input <- input[grepl("\\(.*\\).*;", input)] # leave only the elements that are newick strings
		if(length(input)>1){
			message("Input has several newick strings. Only the first one will be used.")
		}
  	phy.new.in <- ape::collapse.singles(phytools::read.newick(text = gsub(" ", "_", input[1])))
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

#' checks if input is a datelifeQuery object, otherwise it uses make_datelife_query to process it
#' @param datelife_query A datelifeQuery object, output of make_datelife_query function
#' @inheritParams datelife_search
#' @inheritDotParams make_datelife_query -input
#' @export
datelife_query_check <- function(datelife_query = NULL, ...){
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
		# badformat <- FALSE  # useful for next block
	}
	return(datelife_query)
}
