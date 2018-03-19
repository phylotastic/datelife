# Functions to perform a molecular dating analysis on a set of input taxa, using calibrations from chronograms in a database.

#' Use all calibrations from chronograms in a database to date a tree.
#' @param phy A phylo object
#' @param partial If TRUE, use source trees even if they only match some of the desired taxa
#' @param use_tnrs If TRUE, use OpenTree's services to resolve names. This can dramatically improve the chance of matches, but also take much longer
#' @param approximate_match If TRUE, use a slower TNRS to correct mispellings, increasing the chance of matches (including false matches)
#' @param cache The cached set of chronograms and other info from data(opentree_chronograms)
#' @param expand How much to expand by each step to get consistent calibrations
#' @param giveup How many expansion to try before giving up
#' @inheritParams datelife_search
#' @return list with chronogram, original calibrations, and expanded calibrations
#' @export
#' @details
#' This will try to use the calibrations as fixed ages.
#' If that fails (often due to conflict between calibrations), it will expand the range of the minage and maxage and try again. And repeat.
#' expand sets the expansion value: should be between 0 and 1
use_all_calibrations <- function(phy = make_bold_otol_tree(c("Rhea americana",  "Struthio camelus", "Gallus gallus"), chronogram = FALSE, verbose = FALSE), partial = TRUE, use_tnrs = FALSE, approximate_match = TRUE, update_cache = FALSE, cache = get("opentree_chronograms"), expand = 0.1, giveup = 100, verbose = FALSE) {
	# if(!geiger::is.phylo(phy)){
	# 	cat("phy argument must be a phylo object.", "\n")
	# 	stop()
	# }
	#	[code above] useful if we decide to allow newick as input, would need to add some newick to phylo code in here... OK! It is now implemented in input_process function:
	input <- input_process(input = phy, verbose = verbose)
	calibrations.df <- get_all_calibrations(input = gsub('_', ' ', phy$tip.label), partial = partial, use_tnrs = use_tnrs, approximate_match = approximate_match, update_cache = update_cache, cache = cache, verbose = verbose)
	phy$tip.label <- gsub(' ', '_', phy$tip.label) #underscores vs spaces: the battle will never end.
	calibrations.df$taxonA <- gsub(' ', '_', calibrations.df$taxonA)
	calibrations.df$taxonB <- gsub(' ', '_', calibrations.df$taxonB)
	calibrations.df <- calibrations.df[which(calibrations.df$taxonA %in% phy$tip.label),]
	calibrations.df <- calibrations.df[which(calibrations.df$taxonB %in% phy$tip.label),]
	original.calibrations.df <- calibrations.df
	chronogram <- NULL
	try(chronogram <- geiger::PATHd8.phylo(phy, calibrations.df), silent = TRUE)
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

#' Get all calibrations from chronograms in a database (specified in cache).
#' @param input vector of names, a newick string, or a phylo object
#' @param partial Boolean; default TRUE: use source trees even if they only match some of the desired taxa
#' @param use_tnrs Boolean; default False. If TRUE, use OpenTree's services to resolve names. This can dramatically improve the chance of matches, but also take much longer
#' @param approximate_match Boolean; default TRUE: use a slower TNRS to correct mispellings, increasing the chance of matches (including false matches)
#' @param cache The cached set of chronograms and other info from data(opentree_chronograms)
#' @inheritParams datelife_search
#' @return data.frame of calibrations
#' @export
get_all_calibrations <- function(input = c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"), partial = TRUE, use_tnrs = FALSE, approximate_match = TRUE, update_cache = FALSE, cache = get("opentree_chronograms"), verbose = FALSE) {
	datelife_phylo <- datelife_search(input = input, partial = partial, use_tnrs = use_tnrs, approximate_match = approximate_match, update_cache = update_cache, cache = cache, summary_format = "phylo.all", verbose = verbose)
	constraints.df <- data.frame()
	for (i in sequence(length(datelife_phylo))) {
		local.df <- geiger::congruify.phylo(reference = datelife_phylo[[i]], target = datelife_phylo[[i]], scale = NA)$calibrations
		local.df$reference <- names(datelife_phylo)[i]
		if(i == 1) {
			constraints.df <- local.df
		} else {
			constraints.df <- rbind(constraints.df, local.df)
		}
	}
	return(constraints.df)
}

#' Use Barcode of Life data to get branch lengths on the OToL tree of a set of taxa.
#' @inheritParams datelife_search
#' @param marker A character vector with the name of the gene from Barcode of Life Data Systems (BOLD) to be used for branch length estimation.
#' @param otol_version Version of OToL to use
#' @param chronogram Boolean; default to TRUE:  branch lengths represent time estimated with ape::chronoMPL. If FALSE, branch lengths represent relative substitution rates estimated with phangorn::acctran.
#' @param doML Boolean; if TRUE, does ML branch length optimization with phangorn::optim.pml
#' @inheritParams make_datelife_query
#' @return A phylogeny with branch lengths proportional to relative substitution rate.
#' @export
make_bold_otol_tree <- function(input = c("Rhea americana",  "Struthio camelus", "Gallus gallus"), use_tnrs = FALSE, approximate_match = TRUE, marker = "COI", otol_version = "v2", chronogram = TRUE, doML = FALSE, get_spp_from_taxon = FALSE, verbose = FALSE) {
	#otol returns error with missing taxa in v3 of rotl
	input <- datelife_query_check(datelife_query = input, use_tnrs = use_tnrs, approximate_match = approximate_match, get_spp_from_taxon = get_spp_from_taxon, verbose = verbose)
	input <- input$cleaned_names
	if (verbose) {
		cat("Searching", marker, "sequences for these taxa in BOLD...", "\n")
	}
	xx <- seq(1, length(input), 250)
	yy <- xx+249
	yy[length(xx)] <- length(input)
	# if(length(input)%%250 != 0) {
	# 	yy[length(xx)] <- length(input)
	# }
	sequences <- c()
	for (i in seq_len(length(xx))){
		sequences <- rbind(sequences, bold::bold_seqspec(taxon = input[xx[i]:yy[i]], marker = marker)) # bold::bold_seqspec function only allows searches of up to 335 names, ater that, it gives following Error: Request-URI Too Long (HTTP 414)
	}
	if(length(sequences) == 1) {  # it is length == 80 when there is at least 1 sequence available; if this is TRUE, it means there are no sequences in BOLD for the set of input taxa.
		if (verbose) cat("No sequences found in BOLD for input taxa...", "\n")
		# if (!use_tnrs) cat("Setting use_tnrs = TRUE might change this, but it can be slowish.", "\n")
		warning("Names in input do not match BOLD specimen records. No tree was constructed.")
		return(NA)
	}
	sequences$nucleotide_ATGC <- gsub("[^A,T,G,C]", "", sequences$nucleotides)  # preserve good nucleotide data, i.e., only A,T,G,C
	sequences$nucleotide_ATGC_length <- unlist(lapply(sequences$nucleotide_ATGC, nchar))  # add a column in data.frame, indicating the amount of good information contained in sequences#nucelotides (ATGC)
	if (verbose) cat("\t", "OK.", "\n")

	phy <- get_otol_synthetic_tree(input = input, otol_version = otol_version)

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
		taxon.index <- which(grepl(i, sequences$species_name))
		# if there are no sequences from any taxon, taxon.index is empty
		# but we make sure this is filtered steps before
		if (length(taxon.index) > 0){
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
		# if (use_tnrs == FALSE) cat("Setting use_tnrs = TRUE might change this, but it is time consuming.", "\n")
		return(NA)
	}
	if(length(taxa.to.drop) > 0) {
		if (verbose) {
			taxa.to.drop.print <- paste(taxa.to.drop, collapse = " | ")
			cat("No", marker, "sequences found for", taxa.to.drop.print, "...", "\n", "\t", "Dropping taxa from tree.", "\n")
		}
		#warning("No ", marker, " sequences found for ", taxa.to.drop.print, "...", "\n", "\t", "Taxa dropped from tree.")
		taxa.to.drop <- gsub(" ", "_", taxa.to.drop)
		phy <- ape::drop.tip(phy, taxa.to.drop)
	}
	if (verbose) {
		cat("Aligning with MAFFT...", "\n")
	}
	alignment <- ape::as.DNAbin(final.sequences)
	alignment <- phangorn::as.phyDat(ips::mafft(alignment))
	if (verbose) {
		cat( "\t", "OK.", "\n", "Estimating BOLD-OToL tree...", "\n")
	}
	pml.object <- phangorn::pml(phangorn::acctran(phy, alignment), data = alignment)
	phy <- pml.object$tree
	if(!ape::is.binary.tree(pml.object$tree)){
		if (verbose) {
			cat("\t", marker, " sequence data available generates a non-dichotomous tree...", "\n", "\t", "Resolving with multi2di...", "\n")
		}
		pml.object$tree <- ape::multi2di(pml.object$tree)
		phy <- pml.object$tree
	}
	if (verbose) {
		cat("\t", "OK.", "\n")
	}
	if (chronogram) {
		if (verbose) {
			cat("Dating BOLD-OToL tree with chronoMPL...", "\n")
		}
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
	if (verbose) {
		cat("Done.", "\n")
	}
	return(phy)
}

#' Taxon name resolution service applied to input by batches
#' @inheritParams datelife_search
#' @param reference_taxonomy A character vetor specifying the reference taxonomy to use for tnrs.
#' @return A data.frame from tnrs_match_names function
#' @export
input_tnrs <- function(input, reference_taxonomy = "otl"){  # we can add other reference taxonomies in the future
	if(reference_taxonomy == "otl"){
		xx <- seq(1, length(input), 250)
		yy <- xx+249
		yy[length(xx)] <- length(input)
		rr <- c()
		for (i in seq_len(length(xx))){
			rr <- rbind(rr, suppressWarnings(rotl::tnrs_match_names(names = input[xx[i]:yy[i]])))
		}
		# rr <- rotl::tnrs_match_names(names = input)
		# rr has the same order as input
		# when names are not matched it gives a warning: NAs introduced by coercion, so:
		rr <- rr[!is.na(rr$unique_name),]  # gets rid of names not matched with rotl::tnrs_match_names; otherwise rotl::tol_induced_subtree won't run
	}
	return(rr)
}
#' Gets Open Tree of Life synthetic tree of a set of lineages.
#' @inheritParams datelife_search
#' @inheritParams make_bold_otol_tree
#' @return A phylo object
#' @export
get_otol_synthetic_tree <- function(input, otol_version = "v2"){
	rr <- input_tnrs(input = input, reference_taxonomy = "otl")  # processes input with rotl::tnrs_match_names function by batches, so it won't choke
	phy <- ape::multi2di(rotl::tol_induced_subtree(ott_ids = rr$ott_id, label_format = "name",  otl_v = otol_version))
	phy$tip.label <- gsub(".*_ott","", phy$tip.label)  # leaves only the ott_id as tip.label, it's safer than matching by name
	# when there are synonyms among the input names, phy will conserve the accepted name (rr$uniqe_name) instead of the original query name from input (rr$search_string)
	# this produces an error downstream, while using phangorn::pml()
	# to avoid this error, we replace the unique name by the original query name in phy$tip.label:
	mm <- match(phy$tip.label, gsub(" ","_", rr$ott_id))  # this gets the order of tip labels in phy
	phy$tip.label <- gsub(" ","_", input[mm])  # this overlaps the original query over phy$tip.labels in the correct order
	return(phy)
}
