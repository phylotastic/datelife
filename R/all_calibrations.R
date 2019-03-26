# Functions to perform a molecular dating analysis on a set of input taxa, using calibrations from chronograms in a database.

#' Use all calibrations from chronograms in a database to date a tree.
#' @param phy A phylo object
#' @param all_calibrations A data frame of calibrations from get_all_calibrations function
#' @param partial If TRUE, use source trees even if they only match some of the desired taxa
#' @param use_tnrs If TRUE, use OpenTree's services to resolve names. This can dramatically improve the chance of matches, but also take much longer
#' @param approximate_match If TRUE, use a slower TNRS to correct mispellings, increasing the chance of matches (including false matches)
#' @param cache The cached set of chronograms and other info from data(opentree_chronograms)
#' @param expand How much to expand by each step to get consistent calibrations
#' @param giveup How many expansions to try before giving up
#' @inheritParams datelife_search
#' @return list with chronogram, original calibrations, and expanded calibrations
#' @export
#' @details
#' This will try to use the calibrations as fixed ages.
#' If that fails (often due to conflict between calibrations), it will expand the range of the minage and maxage and try again. And repeat.
#' expand sets the expansion value: should be between 0 and 1
use_all_calibrations <- function(phy = NULL,
	# enhance: use congruification to exapnd calibrations.
	all_calibrations = NULL, partial = TRUE, use_tnrs = FALSE, approximate_match = TRUE,
	update_cache = FALSE, cache = get("opentree_chronograms"), expand = 0.1,
	giveup = 100, verbose = FALSE) {
		if(is.null(phy)){ # just to run an example:
			phy <- make_bold_otol_tree(c("Rhea americana", "Struthio camelus", "Gallus gallus"), chronogram = TRUE, verbose = FALSE)
			if(!inherits(phy, "phylo")){
				phy <- get_dated_otol_induced_subtree(input = c("Rhea americana", "Struthio camelus", "Gallus gallus"))
			}
		} else {
			input <- input_process(input = phy, verbose = verbose)
			phy <- input
		}
		# phy must be a tree, check this
		if(!inherits(phy, "phylo")){
			message("phy is not a phylo object")
			return(NA)
		}
		# remove singleton nodes in phy:
		phy <- ape::collapse.singles(phy)
		if(is.null(all_calibrations)){
			calibrations.df <- get_all_calibrations(input = gsub('_', ' ', phy$tip.label),
			partial = partial, use_tnrs = use_tnrs, approximate_match = approximate_match,
			update_cache = update_cache, cache = cache, verbose = verbose)
		} else {
			# enhance: add a check for object structure
			# enhance: check that input names are in calibrations.df
			calibrations.df <- all_calibrations
		}
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
			while(is.null(chronogram) & attempts < giveup) {
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
				message("Dates are even more approximate than usual: had to expand constraints to have them agree.", "\n")
			}
		}
		# get_nodeage_distribution(ages = calibrations.df[, c("MinAge", "MaxAge")],
		# 						 taxonA = calibrations.df[, c("taxonA", "taxonA")],
		# 						 taxonB = calibrations.df[, c("taxonB", "taxonB")],
		# 						 phy = phy)
		return(list(phy = chronogram, calibrations.df = calibrations.df, original.calibrations.df = original.calibrations.df))
}
# figure out what part of patristic_matrix_to_phylo can be reused:
# get_nodeage_distribution() <- function(ages, taxonA, taxonB, phy){
#
# }
#' Get all calibrations from chronograms in a database (specified in cache).
#' @param input vector of names, a newick string, a phylo or multiPhylo object, a datelifeResult object
#' @param partial Boolean; default TRUE: use source trees even if they only match some of the desired taxa
#' @param use_tnrs Boolean; default False. If TRUE, use OpenTree's services to resolve names. This can dramatically improve the chance of matches, but also take much longer
#' @param approximate_match Boolean; default TRUE: use a slower TNRS to correct mispellings, increasing the chance of matches (including false matches)
#' @param cache The cached set of chronograms and other info from data(opentree_chronograms)
#' @inheritParams datelife_search
#' @return A data_frame of calibrations
#' @export
get_all_calibrations <- function(input = c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"), partial = TRUE, use_tnrs = FALSE, approximate_match = TRUE, update_cache = FALSE, cache = get("opentree_chronograms"), verbose = FALSE) {
	if(!inherits(input, "datelifeResult") & !inherits(input, "phylo") & !inherits(input, "multiPhylo")){
		datelife_phylo <- datelife_search(input = input, partial = partial, use_tnrs = use_tnrs, approximate_match = approximate_match, update_cache = update_cache, cache = cache, summary_format = "phylo_all", verbose = verbose)
	}
	# inherits(datelife_phylo, "datelifeResult")
	if(inherits(input, "datelifeResult")){
		datelife_phylo <- suppressMessages(summarize_datelife_result(datelife_result = input, summary_format = "phylo_all"))
	}
	if(inherits(input, "phylo")){
		datelife_phylo <- list(input)
	}
	if(inherits(input, "multiPhylo")){
		datelife_phylo <- input
	}
	constraints.df <- data.frame() # we cannot set an empty data frame because nrow depends on the number of nodes available on each tree
	for (i in seq(length(datelife_phylo))) {
		local.df <- suppressWarnings(geiger::congruify.phylo(reference = datelife_phylo[[i]], target = datelife_phylo[[i]], scale = NA))$calibrations
		# suppressedWarnings bc of meesage when running geiger::congruify.phylo(reference = datelife_phylo[[i]], target = datelife_phylo[[i]], scale = NA)
		# 		Warning message:
		# In if (class(stock) == "phylo") { :
		# the condition has length > 1 and only the first element will be used
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
#' @inheritParams get_otol_synthetic_tree
#' @param marker A character vector with the name of the gene from Barcode of Life Data Systems (BOLD) to be used for branch length estimation.
#' @param chronogram Boolean; default to TRUE:  branch lengths represent time estimated with ape::chronoMPL. If FALSE, branch lengths represent relative substitution rates estimated with phangorn::acctran.
#' @param doML Boolean; if TRUE, does ML branch length optimization with phangorn::optim.pml
#' @inheritParams make_datelife_query
#' @return A phylogeny with branch lengths proportional to relative substitution rate.
#' @details If input is a phylo object, it is used as backbone. If it is a character vector of taxon names, an induced OToL tree is used as backbone.
#' @export
make_bold_otol_tree <- function(input = c("Rhea americana",  "Struthio camelus", "Gallus gallus"), use_tnrs = FALSE, approximate_match = TRUE, marker = "COI", otol_version = "v3", chronogram = TRUE, doML = FALSE, get_spp_from_taxon = FALSE, verbose = FALSE) {
	# enhance: add an input check here to accept newick strings too
	if(inherits(input, "phylo")){
		phy <- input
		input <- phy$tip.label
	} else {
		phy <- get_otol_synthetic_tree(input = input, otol_version = otol_version, use_tnrs = use_tnrs, approximate_match = approximate_match, get_spp_from_taxon = get_spp_from_taxon)
		#otol returns error with missing taxa in v3 of rotl
		input <- phy$tip.label
	}
	if (verbose) {
		message("Searching ", marker, " sequences for these taxa in BOLD...")
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
		if (verbose) {
			message("No sequences found in BOLD for input taxa...")
		}
		# if (!use_tnrs) cat("Setting use_tnrs = TRUE might change this, but it can be slowish.", "\n")
		warning("Names in input do not match BOLD specimen records. No tree was constructed.")
		return(NA)
	}
	sequences$nucleotide_ATGC <- gsub("[^A,T,G,C]", "", sequences$nucleotides)  # preserve good nucleotide data, i.e., only A,T,G,C
	sequences$nucleotide_ATGC_length <- unlist(lapply(sequences$nucleotide_ATGC, nchar))  # add a column in data frame, indicating the amount of good information contained in sequences#nucelotides (ATGC)
	if (verbose) {
		message("\t", "OK.")
	}

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
		if (verbose) {
			message("BOLD sequences found only for one input name: ", input[which(!input %in% taxa.to.drop)], ".","\n","\t", "Cannot construct a tree." )
		}
		message("Not enough sequences available in BOLD. No tree was constructed.")
		# if (use_tnrs == FALSE) cat("Setting use_tnrs = TRUE might change this, but it is time consuming.", "\n")
		return(NA)
	}
	if(length(taxa.to.drop) > 0) {
		if (verbose) {
			taxa.to.drop.print <- paste(taxa.to.drop, collapse = " | ")
			message("No ", marker, " sequences found for ", taxa.to.drop.print, ".", "\n", "\t", "Dropping taxa from tree.")
		}
		#warning("No ", marker, " sequences found for ", taxa.to.drop.print, "...", "\n", "\t", "Taxa dropped from tree.")
		taxa.to.drop <- gsub(" ", "_", taxa.to.drop)
		phy <- ape::drop.tip(phy, taxa.to.drop)
	}
	if (verbose) {
		message("Aligning with MAFFT...")
	}
	alignment <- ape::as.DNAbin(final.sequences)
	alignment <- phangorn::as.phyDat(ips::mafft(alignment))
	if (verbose) {
		message( "\t", "OK.", "\n", "Estimating BOLD-OToL tree...")
	}
	xx <- phangorn::acctran(phy, alignment)
	pml.object <- phangorn::pml(xx, data = alignment)
	phy <- pml.object$tree
	if(!ape::is.binary.tree(pml.object$tree)){
		if (verbose) {
			message("\t", marker, " sequence data available generates a non-dichotomous tree.", "\n", "\t", "Resolving with multi2di...")
		}
		pml.object$tree <- ape::multi2di(pml.object$tree)
		phy <- pml.object$tree
	}
	if (verbose) {
		message("\t", "OK.")
	}
	if (chronogram) {
		if (verbose) {
			message("Dating BOLD-OToL tree with chronoMPL...")
		}
		pml.object$tree <- ape::chronoMPL(pml.object$tree, se = FALSE, test = FALSE)
		phy <- pml.object$tree
		if (verbose) {
			message("\t", "OK.")
		}
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
		message("Done.")
	}
	return(phy)
}

#' Gets Open Tree of Life synthetic tree of a set of lineages.
#' @inheritParams check_ott_input
#' @param otol_version Version of OToL to use
#' @inheritDotParams make_datelife_query -input
#' @return A phylo object
#' @export
get_otol_synthetic_tree <- function(input = NULL, ott_ids = NULL, otol_version = "v2", ...){
	# input <- birds_yellowstone
	input_ott_match <- suppressMessages(check_ott_input(input, ott_ids, ...))
	if(length(input_ott_match) < 2){
		message("At least two valid names or numeric ott_ids are needed to get a tree")
		return(NA)
	}
	# enhance: we might need a check of ott_id elements, are they all numeric, are there no NAs, etc.
	# also, another check here, are all ott_ids from valid taxa? this is checked with get_ott_children, but from other functions we should check This
	# enhance: add a class to get_ott_children outputs so its easier to check here if all ott ids are valid taxa, indicate if it has been cleaned, maybe within an atrribute
	# system.time({sapply(rotl::taxonomy_taxon_info(df$ott_id), "[", "flags")})
	# system.time({tnrs_match(rownames(df))}) # this one is faster
	phy <- tryCatch(suppressWarnings(rotl::tol_induced_subtree(ott_ids = input_ott_match, label_format = "name",  otl_v = otol_version)),
					error = function(e){
						message("Some or all input taxa are absent from OToL synthetic tree, look for invalid taxa and clean them from input")
						# this will happen if there are some extinct taxa, barren or any invalid taxa in taxonomy
						# enhance: to avoid it, clean invalid taxa at the beginning
						NA })
	if(length(phy) == 1){
		return(phy)
	}
	if(!ape::is.binary(phy)){
		message(paste0("OToL synthetic tree of input taxa is not fully resolved (",
		phy$Nnode, " nodes/", length(phy$tip.label), " tips)."))
		message("Resolving at random...")
		phy <- ape::multi2di(phy)
	}
	# example of weird behaviour on tip labeling from otol:
	# tnrs <- rotl::tnrs_match_names(c("Staphylococcus aureus", "Bacillus subtilis", "Neisseria meningitidis"))
	# tol_sub <- rotl::tol_induced_subtree(ott_ids = tnrs$ott_id)
	# curl -X POST https://api.opentreeoflife.org/v3/tree_of_life/induced_subtree -H "content-type:application/json" -d '{"ott_ids":[1090496, 1084928, 611778]}'
	# tol_subpost <- ape::read.tree(text = "((((((((((((((((((((((((Bacillus_subtilis_ott1084928)mrcaott6603ott41501)mrcaott6603ott219996)mrcaott6603ott51434)mrcaott6603ott7907)mrcaott1417ott6603)mrcaott1417ott130524)mrcaott609ott1417)mrcaott609ott695050)mrcaott609ott3174)mrcaott218ott609,(((((((mrcaott4291ott4292)mrcaott3204ott3205)mrcaott3204ott109201)mrcaott3204ott935729)mrcaott3204ott460129)Staphylococcus_ott720488)Staphylococcaceae_ott949800)mrcaott1420ott85485)mrcaott218ott1420)mrcaott218ott9340)mrcaott218ott1491)mrcaott217ott218)mrcaott217ott1456)mrcaott47ott53)mrcaott47ott63791)mrcaott47ott10726)mrcaott47ott184124)mrcaott47ott64165)mrcaott47ott2035)mrcaott47ott19087)mrcaott47ott63363,((((((((((((((((((((((((((((Neisseria_meningitidis_ott611778)mrcaott22944ott239089)mrcaott22944ott237408)mrcaott22944ott239091)mrcaott22944ott469870)mrcaott22944ott233491)mrcaott22944ott279619)mrcaott22944ott279613)mrcaott22944ott279628)mrcaott22944ott279610)mrcaott14134ott22944)mrcaott14134ott67965)mrcaott14134ott1011641)Neisseria_ott611812)mrcaott5074ott139204)Neisseriaceae_ott286853)Neisseriales_ott779197)mrcaott90ott5074)mrcaott90ott103)mrcaott90ott11872)mrcaott90ott191429)mrcaott89ott90)mrcaott89ott3892)mrcaott50ott89)mrcaott50ott21523)mrcaott50ott6117)mrcaott50ott107113)mrcaott50ott1100)mrcaott50ott73)mrcaott47ott50;")

	# now include ott_ids as phy element:
	if(is.null(phy$ott_ids)){
		phy$ott_ids <- input_ott_match[match(phy$tip.label, gsub(" ", "_",
		names(input_ott_match)))] # will give NAs where there are mrcaotts, but fixed after
	}
	mrca_index <- grep("mrcaott", phy$tip.label)
	if(length(mrca_index) > 0){
		mrcaott_names <- unlist(lapply(phy$tip.label[mrca_index], recover_mrcaott))
		phy$tip.label[mrca_index] <- names(mrcaott_names)
		phy$ott_ids[mrca_index] <- mrcaott_names
	}
	phy$ott_ids <- as.numeric(phy$ott_ids)
	return(phy)
}

#' Get an otol induced dated subtree from your set of queried taxa
#' @inheritParams check_ott_input
#' @inheritDotParams make_datelife_query -input
#' @return A phylo object with edge length proportional to time in Myrs. It will return NA if any ott_id is invalid.
#' @export
#' @details otol dated tree from Stephen Smith's otol scaling service.
# # ' @examples
#' # if you want to make an ltt plot of a dated OToL tree you'll need to get rid of singleton nodes with ape::collapse.singles
#' # and also probably do phytools::force.ultrametric
get_dated_otol_induced_subtree <- function(input = NULL, ott_ids = NULL, ...){
	# for debugging:
	# utils::data(threebirds.rda)
	# input <- threebirds_query$cleaned_names
	# ott_ids <- threebirds_query$ott_ids
	input_ott_match <- suppressMessages(check_ott_input(input, ott_ids, ...))
	if(length(input_ott_match) < 2){
		message("At least two valid names or numeric ott_ids are needed to get a tree")
		return(NA)
	}
  pp <- tryCatch(httr::POST("http://141.211.236.35:10999/induced_subtree",
						body = list(ott_ids = input_ott_match),
						encode = "json", httr::timeout(10)), error = function(e) NA)
	if(length(pp) > 1){ # this means it retrieved a tree succesfully
		pp <- httr::content(pp)
		rr <- httr::POST("http://141.211.236.35:10999/rename_tree",
		          body = list(newick = pp$newick),
		          encode = "json", httr::timeout(10))
		rr <- httr::content(rr)
		rr <- ape::read.tree(text = rr$newick)
		rr$ott_ids <- ape::read.tree(text = pp$newick)$tip.label
	  return(rr)
	}	else { # this means it errored and we return NA
		return(NA)
	}
}
