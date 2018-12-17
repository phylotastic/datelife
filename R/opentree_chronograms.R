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
#' @param max_tree_count Numeric indicating the mas number of trees to be cached. For testing purposes only.
#' @return A list with elements for the trees, authors, curators, and study ids
#' @export
get_otol_chronograms <- function(verbose = FALSE, max_tree_count = 500) {
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
	while(tree.count < max_tree_count){ # the only purpose for this conditional is to run the testhat for this function.
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
									data.frame(new.tree$tip.label)
									data.frame(try.tree$tip.label)
									if(verbose) {
										message("tree_id = '", tree.id, "', study_id = '", study.id, "'")
									}
									if(!is.null(new.tree) & phylo_has_brlen(phy = new.tree)) {
											try.tree <- clean_ott_chronogram(new.tree) # will give NA if just one or none tip labels are mapped to ott
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
											} # else {
											# 		# here goes the case with none or just one mapped taxa
											# 		new.tree <- map_ott_chronogram(new.tree)
											# }
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
	} # end of while(tree.count < max_tree_count) conditional
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
	# enhance: add taxonomic nodelabels to trees object here
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
	if(length(which(grepl("not.mapped", phy$tip.label))) > 0) {
		warning("tree failed over having not mapped taxa that should have been checked")
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
#' For now it 1) checks unmapped taxa and maps them with map_tiplabels_ott, 2) roots the chronogram if unrooted
#'
#' @inheritParams phylo_check
#' @return A cleaned up phylo object
#' @export
clean_ott_chronogram <- function(phy) {
	# phy <- problems[[5]]
	# homogenize everything to blanks, it speeds up tnrs_match_names bc no approximate matching needed
	phy.ori <- phy
	phy <- phylo_tiplabel_underscore_to_space(phy)
	# tip.label <- phy$tip.label
	bad.taxa <- unique(c(which(nchar(phy$tip.label) <= 2), which(grepl("not.mapped", phy$tip.label)))) # numeric of indices
	phy$tip.label[bad.taxa] <- sub(".*-.", "", phy$tip.label[bad.taxa])  # this gets the original label and gets rid of the not.mapped tag
	phy$tip.label[bad.taxa] <- gsub("aff ", "", phy$tip.label[bad.taxa])  # removes aff tag
    phy$tip.label[bad.taxa][stringr::str_count(phy$tip.label[bad.taxa], " ")>=2] <- gsub("^([^ ]* [^ ]*) .*$", "\\1", phy$tip.label[bad.taxa][stringr::str_count(phy$tip.label[bad.taxa], " ")>=2])
    # grepl("^([^_]*_[^_]*)_.*$", tip.labels[stringr::str_count(tip.labels, "_")>=2])  # this works with underscores
	# drop duplicated tips now:
	tipstodrop <- c()
	# identify not mapped taxa that have mapped duplicates first
	# in this case, drop the not mapped tips and keep the mapped ones.
	# first identify the not mapped tips that have mapped duplicates:
	cond <- match(unique(phy$tip.label[bad.taxa]), phy$tip.label[-bad.taxa])
	if(any(!is.na(cond))){
		# the following only gets the first match:
		# mm <- match(unique(phy$tip.label[bad.taxa])[!is.na(cond)], phy$tip.label[bad.taxa])
		# this gets all bad.taxa that match:
		mm <- match(phy$tip.label[bad.taxa], unique(phy$tip.label[bad.taxa])[!is.na(cond)])
		tipstodrop <- c(tipstodrop, bad.taxa[!is.na(mm)])
		bad.taxa <- bad.taxa[is.na(mm)]
	}
	# now identify all remaining duplicated and not mapped
	dd <- duplicated(phy$tip.label[bad.taxa])
	tipstodrop <- c(tipstodrop, bad.taxa[dd])
	phy <- ape::drop.tip(phy, tip = unique(tipstodrop))
	bad.taxa <- bad.taxa[!dd]
	phy$mapped <- rep("ott", length(phy$tip.label))
	if(length(bad.taxa) > 0) {
		phy <- tnrs_match.phylo(phy, tip = bad.taxa)
	}
	if(!ape::is.rooted(phy) & ape::is.ultrametric(phy, option = 2)) {
		phy$root.edge <- 0
	}
	return(phy)
}



#' Map opentree_chronograms tip_label's Open Tree of Life Taxonomy to nodes.
#' @inheritParams tree_fix_brlen
#' @return A cleaned up phylo object
#' @export
map_nodes_ott <- function(tree){
	phy <- tree_check(tree = tree)
	got <- get_ott_lineage(tree)
	return(got)
}
