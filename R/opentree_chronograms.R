#library(devtools)
#install_github("ropensci/rotl", dependencies = TRUE, build_vignette=FALSE)
#library(rotl)
#library(ape)
#library(knitcitations)

#' update Open Tree of Life cache
#' @param save Boolean; default TRUE: save all chronograms from Open Tree of Life to an RData file (default to opentree_chronograms.RData)
#' @inheritParams save_otol_chronograms
#' @inherit get_otol_chronograms return
#' @inherit save_otol_chronograms return
#' @export
update_datelife_cache <- function(save = TRUE, file = "opentree_chronograms.RData", verbose = TRUE){  #, new_studies_only = TRUE
	if (save) {
		cache_updated <- save_otol_chronograms(file = file, verbose = verbose)
	} else {
		cache_updated <- get_otol_chronograms(verbose = verbose)
	}
	return(cache_updated)
}

#' Update all cached files for the package
#'
#' For speed, datelife caches chronograms and other information. Running this (within the checked out version of datelife) will refresh these. Then git commit and git push them back
#' @return nothing
#' @export
update_all_cached <- function() {
	opentree_chronograms <- get_otol_chronograms()
	devtools::use_data(opentree_chronograms, overwrite=TRUE)
	contributor_cache <- make_contributor_cache(outputfile=paste0(tempdir(), '/contributor.rda'))
	devtools::use_data(contributor_cache, overwrite=TRUE)
	treebase_cache <- make_treebase_cache(outputfile=paste0(tempdir(), '/treebase.rda'))
	devtools::use_data(treebase_cache, overwrite=TRUE)
	depositor_cache <- make_all_associations(outputfile=NULL)
	devtools::use_data(depositor_cache, overwrite=TRUE)
}

#' Check for branch lengths in a tree
#' @inheritParams phylo_check
#' @return A TRUE or FALSE
#' @export
phylo_has_brlen <- function(phy) {
	brlen=TRUE
	if(is.null(phy$edge.length)) {
		brlen=FALSE
	}
	return(brlen)
}


# #' Allow trees with duplicate taxa to be read in, courtesy David Winter
# #' @param study_id Open Tree study id
# #' @param tree_id Open Tree tree id
# #' @param tip_label Which Open Tree tree label format you want
# #' @return A phylo object
# get_study_tree_with_dups <- function(study_id, tree_id, tip_label="ot:otttaxonname") {
# 	tr <- rotl:::.get_study_tree(study_id=study_id, tree_id=tree_id, tip_label=tip_label, format="newick")
# 	phy <- ape::read.tree(text=gsub(" ", "_", tr))
# 	phy$tip.label <- gsub("'", "", gsub("_", " ", phy$tip.label))
# 	return(	phy)
# }

#' Get all chronograms from Open Tree of Life
#' @param verbose If TRUE, give updates to the user
#' @return A list with elements for the trees, authors, curators, and study ids
#' @export
get_otol_chronograms <- function(verbose = FALSE) {
	if(verbose) {
		options(warn = 1)
	}
	chronogram.matches <- rotl::studies_find_trees(property = "ot:branchLengthMode", value = "ot:time", verbose = TRUE, detailed = TRUE)
	trees <- list()
	authors <- list()
	curators <- list()
	studies <- list()
	dois <- list()
	tree.count <- 0
	bad.ones <- c()
	# utils::opentree_chronograms

	for (study.index in sequence(dim(chronogram.matches)[1])) {
		if(verbose) {
			message("Downloading tree(s) from study ", study.index, " of ", dim(chronogram.matches)[1])
		}
		for(chrono.index in sequence((chronogram.matches$n_matched_trees[study.index]))) {
			study.id <- chronogram.matches$study_ids[study.index]
			# if(study.id %in% opentree_chronograms$studies)  # unlist(opentree_chronograms$studies) if new_studies_only
	#	new.tree <- get_study_tree(study_id=study.id, tree_id=tree.id, tip_label='ott_taxon_name')
			new.tree2 <- new.tree <- NULL
			tree.id <- strsplit(chronogram.matches$match_tree_ids[study.index], ", ")[[1]][chrono.index]
			potential.bad <- paste("tree_id='", tree.id, "', study_id='", study.id, "'", sep="")

			if(!grepl("\\.\\.\\.", tree.id) & !is.na(tree.id)) { #to deal with ellipsis bug
				#try(new.tree <- datelife:::get_study_tree_with_dups(study_id=study.id,tree_id=tree.id ))
				#try(new.tree <- rotl::get_study_subtree(study_id=study.id,tree_id=tree.id, tip_label="ott_taxon_name", subtree_id="ingroup")) #only want ingroup, as it's the one that's been lovingly curated.
				try(new.tree <- rotl::get_study_tree(study_id=study.id,tree_id=tree.id, tip_label="ott_taxon_name")) #would like to dedup; don't use get_study_subtree, as right now it doesn't take tip_label args
				#try(new.tree <- rotl::get_study_tree(study_id=study.id,tree_id=tree.id, tip_label="ott_taxon_name"))
				if(verbose) {
					message("tree_id = '", tree.id, "', study_id = '", study.id, "'")
				}
				if(!is.null(new.tree) & phylo_has_brlen(phy = new.tree)) {
					try.tree <- clean_ott_chronogram(new.tree)
					if(phylo_has_brlen(phy = try.tree)) {
						new.tree <- try.tree
						if(is_good_chronogram(new.tree)) {
							new.tree$tip.label <- gsub('_', ' ', new.tree$tip.label)
							if(verbose) {
								message("\t", "has tree with branch lengths")
							}
							doi <- NULL
							try(doi <- gsub('https?://(dx\\.)?doi.org/', '', attr(rotl::get_publication(rotl::get_study_meta(study.id)), "DOI")))
							authors <- append(authors, NA)
							if(length(doi) == 0){
								warning(paste(study.id, "has no DOI attribute, author names will not be retrieved."))
							} else {
								try(authors[length(authors)] <- list(paste(as.character(knitcitations::bib_metadata(doi)$author))))
							}
							curators <- append(curators, NA)
							try(curators[length(curators)] <- list(rotl::get_study_meta(study.id)[["nexml"]][["^ot:curatorName"]]))
							try(studies <- append(studies, study.id))
							tree.count <- tree.count + 1
							try(dois <- append(dois, chronogram.matches$study_doi[study.index]))
							try(new.tree2 <- rotl::get_study_tree(study_id=study.id,tree_id=tree.id, tip_label="ott_id")) #would like to dedup; don't use get_study_subtree, as right now it doesn't take tip_label args
							new.tree2 <- clean_ott_chronogram(new.tree2)
							new.tree$ott_ids <- gsub("_.*", "", new.tree2$tip.label)
							trees[[tree.count]] <- new.tree
							names(trees)[tree.count] <- rotl::get_publication(rotl::get_study_meta(study.id))[1]
							message("\t", "was good tree")
							potential.bad <- NULL
						}
					} else {
						# here goes the case with none or just one mapped taxa
						new.tree <- map_ott_chronogram(new.tree)
					}
				}
			} else {
				warning("Not all trees could be loaded from this study due to ellipsis bug, https://github.com/ropensci/rotl/issues/85")
			}
			if(!is.null(potential.bad)) {
				bad.ones <- append(bad.ones, potential.bad)
			}
			#save(list=ls(), file="opentree_chronograms.RData")
		}
	}
	if(verbose) {
		message("Problematic combos:")
		message(paste0(utils::capture.output(bad.ones), collapse = "\n"))
	}
	for(i in sequence(length(authors))){
		if(!any(is.na(authors[[i]])) & length(authors[[i]]) > 0){
			Encoding(authors[[i]]) <- "latin1"
			authors[[i]] <- iconv(authors[[i]], "latin1", "UTF-8")
		}
		for(j in sequence(length(curators[[i]]))){
			if(!any(is.na(curators[[i]][[j]])) & length(curators[[i]][[j]]) > 0){
				Encoding(curators[[i]][[j]]) <- "latin1"
				curators[[i]][[j]] <- iconv(curators[[i]][[j]], "latin1", "UTF-8")
			}
		}
	}
	# add taxonomic nodelabels to trees object here
	result <- list(trees = trees, authors = authors, curators = curators, studies = studies, dois = dois)
	return(result)
}

#' Save all chronograms from Open Tree of Life
#' @param file Path including file name
#' @param verbose If TRUE, give status updates to the user
#' @return None
#' @export
save_otol_chronograms <- function(file = "opentree_chronograms.RData", verbose = FALSE) {
	opentree_chronograms <- get_otol_chronograms(verbose = verbose)
	save(opentree_chronograms, file = file, compress = "xz")
}

#' Check to see that a tree is a valid chronogram
#' @inheritParams phylo_check
#' @return Boolean: TRUE if good tree
#' @export
is_good_chronogram <- function(phy) {
	passing <- TRUE
	if(class(phy) != "phylo") {
		passing <- FALSE
		warning("tree failed over not being class phylo")
	}
	if(ape::Ntip(phy) <= ape::Nnode(phy)) {
		passing <- FALSE
		warning("tree failed over not having more internal nodes than tips")
	}
	if(length(which(grepl("not mapped", phy$tip.label))) > 0) {
		warning("tree failed over having not mapped taxa that should have been purged")
		passing <- FALSE #not cleaned properly
	}
	if(any(is.na(phy$tip.label))) {
		passing <- FALSE
		warning("tree failed over having NA for tips")
	}	else {
		if(min(nchar(phy$tip.label)) <= 2) {
			passing <- FALSE
			warning("tree failed for having names of two or fewer characters")
		}
	}
	if(!ape::is.rooted(phy)) {
		passing <- FALSE
		warning("tree failed over not being rooted")
	}
	if(!ape::is.ultrametric(phy, tol = 0.01, option = 2)) {
		passing <- FALSE
		warning("tree failed over not being ultrametric (NOTE: this condition should be removed for paleo trees)")
	}
	return(passing)
}

#' Clean up some issues with OToL chronograms
#' @inheritParams phylo_check
#' @return A cleaned up phylo object
#' @export
clean_ott_chronogram <- function(phy) {
	original.phy <- phy
	if(class(phy) == "phylo") {
		if(ape::Ntip(phy) > ape::Nnode(phy)) {
			if(any(grepl("_", phy$tip.label))){
				NOT <- "not_mapped"
			} else {
				NOT <- "not mapped"
			}
			bad.taxa <- unique(c(which(nchar(phy$tip.label) <= 2), which(grepl(NOT, phy$tip.label))))
			if(length(bad.taxa) > 0) {
				if(length(bad.taxa) < (ape::Ntip(phy)-1)){#will return trees with as few as two mapped tips
					# but will do nothing to trees with all unmapped tips or just one mapped tip
					if(length(bad.taxa) > 0 & length(bad.taxa) < (ape::Ntip(phy)-1)) {
					phy <- try(ape::drop.tip(phy, bad.taxa))
					if(class(phy) == "try-error") {
						return(original.phy)
					}
				} else { # if none or just one tip is mapped
					# length(bad.taxa) - (ape::Ntip(phy)-1)
					return(NA) # we'll deal with these trees with batch_tnrs_match_names in map_ott_chronograms
				}
			}
			if(!ape::is.rooted(phy) & ape::is.ultrametric(phy, option = 2)) {
				phy$root.edge <- 0
			}
		}
	}
	return(phy)
}
#' Check for unmapped tip labels in opentree_chronograms cache
#' @param opentree_chronograms
#' @return A cleaned up opentree_chronograms object
#' @export


#' Get taxonomy for all species within opentree_chronograms object
#' @inheritParams phylo_check
#' @return A cleaned up phylo object
#' @export
@examples
data(opentree_chronograms)
names(opentree_chronograms)
get_them_all(opentree_chronograms$trees)
get_them_all <- function(trees){
	all <- unique(unlist(sapply(trees, "[", "ott_ids")))
	all[which(is.na(as.numeric(all)))]
	# problems <- opentree_chronograms$trees[sapply(sapply(trees, "[", "ott_ids"), function(x) any(grepl("tip", x)))]
	problems <- opentree_chronograms$trees[sapply(sapply(trees, "[", "tip.label"), function(x) any(grepl("not mapped", x)))]
	length(problems)
	names(problems)
	# [1] "Upham, N.S. & B.D. Patterson. 2015. Evolution of the caviomorph rodents: a complete phylogeny and timetree of living genera. Pp. 63–120 In: Biology of caviomorph rodents: diversity and evolution (A. I. Vassallo & D. Antenucci, eds.). SAREM Series A, Buenos Aires. "
	# [2] "Prashant P. Sharma, Caitlin M. Baker, Julia G. Cosgrove, Joanne E. Johnson, Jill T. Oberski, Robert J. Raven, Mark S. Harvey, Sarah L. Boyer, Gonzalo Giribet, 2018, 'A revised dated phylogeny of scorpions: Phylogenomic support for ancient divergence of the temperate Gondwanan family Bothriuridae', Molecular Phylogenetics and Evolution, vol. 122, pp. 37-45"
	# [3] "Yun-Yun Lv, Kai He, Sebastian Klaus, Rafe M. Brown, Jia-Tang Li, 2018, 'A comprehensive phylogeny of the genus  Kurixalus  (Rhacophoridae, Anura) sheds light on the geographical range evolution of frilled swamp treefrogs', Molecular Phylogenetics and Evolution, vol. 121, pp. 224-232"
	# [4] "Andressa Paladini, Daniela M. Takiya, Julie M. Urban, Jason R. Cryan, 2018, 'New World spittlebugs (Hemiptera: Cercopidae: Ischnorhininae): Dated molecular phylogeny, classification, and evolution of aposematic coloration', Molecular Phylogenetics and Evolution, vol. 120, pp. 321-334"
	# [5] "Renske E. Onstein, H. Peter Linder, 2016, 'Beyond climate: convergence in fast evolving sclerophylls in Cape and Australian Rhamnaceae predates the mediterranean climate', Journal of Ecology, pp. n/a-n/a"	sapply(problems,names)
	clean_ott_chronogram(phy = problems[[5]])
	sapply(problems, "[", "ott_ids")
	tips <- unique(unlist(sapply(problems, "[", "tip.label")))
	tips2 <- gsub(".* - ", "", tips)
	grep(" ", tips2)
	tips2[grepl("[.*[:blank:].*]{2,}", tips2)]
	# I'm useless with grep, so Brian suggested stringr
	tips2[stringr::str_count(tips2, " ")>=2]
	tips_tnrs <- batch_tnrs_match_names(tips2)
	tips_tnrs$unique_name
	sum(is.na(tips_tnrs$unique_name))
	names(tips_tnrs)
	tips_tnrs$search_string[is.na(tips_tnrs$unique)]
	length(tips_tnrs$unique_name)
	datelife_search(input = gsub(".* - ", "", sample(tips, 50)))
	all_got <- get_ott_lineage(ott_id = all)
}

#' Map opentree_chronograms tip_label's Open Tree of Life Taxonomy to nodes.
#' @inheritParams phylo_check
#' @return A cleaned up phylo object
#' @export
map_nodes_ott <- function(tree){
	got <- get_ott_lineage(tree)
}
