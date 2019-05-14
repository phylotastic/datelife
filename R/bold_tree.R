#' Use Barcode of Life Data (BOLD) to get branch lengths on a tree topology.
#' If input is not a tree, then it uses the OToL tree of the set of input taxa.
#' @inheritParams datelife_search
#' @param marker A character vector indicating the gene from BOLD system to be used for branch length estimation.
#' @param chronogram Boolean. If TRUE (default), branch lengths returned are estimated with ape::chronoMPL. If FALSE, branch lengths returned are estimated with phangorn::acctran and represent relative substitution rates .
#' @param doML Boolean; only relevant if chronogram = TRUE. If TRUE, it does ML branch length optimization with phangorn::optim.pml
#' @inheritParams get_otol_synthetic_tree
#' @inheritDotParams get_otol_synthetic_tree
#' @return A phylogeny with branch lengths proportional to relative substitution rate.
#' @details
#' If input is a phylo object, that is used as backbone topology.
#' If input is a character vector of taxon names, an induced OToL tree is used as backbone.
#' If there are not enough sequences to return a tree with branch lengths, it returns
#' a tree with no branch lengths (either the original input topology or the OToL tree
#' obtained for the set of input taxon names.
#' @export
# input <- phyloall[[1]]
make_bold_otol_tree <- function(input = c("Rhea americana",  "Struthio camelus", "Gallus gallus"),
marker = "COI", otol_version = "v3", chronogram = TRUE, doML = FALSE, verbose = FALSE, ...) {
	# enhance: add an input check here to accept newick strings too
	if(inherits(input, "phylo")){
		phy <- input
		input <- phy$tip.label
	} else {
		phy <- get_otol_synthetic_tree(input = input, otol_version = otol_version, ...)
		#otol returns error with missing taxa in v3 of rotl
		input <- phy$tip.label
	}
	if (verbose) {
		message("Searching ", marker, " sequences for these taxa in BOLD...")
	}
	# xx <- seq(1, length(input), 250)
	# yy <- xx+249
	# yy[length(xx)] <- length(input)
	# if(length(input)%%250 != 0) {
	# 	yy[length(xx)] <- length(input)
	# }
	phy$edge.length <- NULL # making sure there are no branch lengths in phy
	phy$tip.label <- gsub(" ", "_", phy$tip.label) # so phangorn::acctran works
	input <- gsub("_", " ", input) # so bold search works
	sequences <- c()
	progression <- utils::txtProgressBar(min = 0, max = length(input), style = 3)
	for (i in seq(length(input))){
		sequences <- rbind(sequences, bold::bold_seqspec(taxon = input[i]))
		# allows up to 335 names, then it gives Error: Request-URI Too Long (HTTP 414)
		# even if marker is specified, it will return other markers,
		# so in here we just retrieve all sequences and filter after
		utils::setTxtProgressBar(progression, i)
	}
	cat("\n") # just to make the progress bar look better
	sequences <- sequences[grepl(marker, sequences$markercode), ] # filter other markers
	if(length(sequences) == 1) {  # it is length == 80 when there is at least 1 sequence available; if this is TRUE, it means there are no sequences in BOLD for the set of input taxa.
		# if (!use_tnrs) cat("Setting use_tnrs = TRUE might change this, but it can be slowish.", "\n")
		message("Names in input do not match BOLD specimen records; no sequences
			were found in BOLD for the set of input taxa. Returning tree with no branch lengths.")
		return(phy)
	}
	sequences$nucleotide_ATGC <- gsub("[^A,T,G,C]", "", sequences$nucleotides)  # preserve good nucleotide data, i.e., only A,T,G,C
	sequences$nucleotide_ATGC_length <- unlist(lapply(sequences$nucleotide_ATGC, nchar))  # add a column in data frame, indicating the amount of good information contained in sequences#nucelotides (ATGC)
	if (verbose) {
		message("\t", "OK.")
	}

	final.sequences <- matrix("-", nrow = length(input), ncol = max(sapply(strsplit(sequences$nucleotides, ""), length)))
	final.sequences.names <- rep(NA, length(input))
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
		message("There are not enough sequences available in BOLD to reconstruct branch lengths. Returning tree with no branch lengths.")
		return(phy)
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
	# if there are only two sequences in the alignment phangorn::acctran will throw an error
	# this usually happens when the input/otol tree has only two tips
	if(length(alignment) <= 2){
		message("There are not enough sequences available in BOLD to reconstruct branch lengths. Returning tree with no branch lengths.")
		return(phy)
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
	phy$tip.label <- gsub(' ', '_', phy$tip.label)
	if (verbose) {
		message("Done.")
	}
	return(phy)
}
