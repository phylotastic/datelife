#library(devtools)
#install_github("ropensci/rotl", dependencies = TRUE, build_vignette=FALSE)
#library(rotl)
#library(ape)
#library(knitcitations)

#' Check for branch lengths in a tree
#' @param x A phylo object
#' @return A TRUE or FALSE
#' @export
HasBrlen <- function(x) {
	brlen=TRUE
	if(is.null(x$edge.length)) {
		brlen=FALSE
	}
	return(brlen)
}


#' Allow trees with duplicate taxa to be read in, courtesy David Winter
#' @param study_id Open Tree study id
#' @param tree_id Open Tree tree id
#' @param tip_label Which Open Tree tree label format you want
#' @return A phylo object
get_study_tree_with_dups <- function(study_id, tree_id, tip_label="ot:otttaxonname") {
	tr <- rotl:::.get_study_tree(study_id=study_id, tree_id=tree_id, tip_label=tip_label, format="newick")
	phy <- ape::read.tree(text=gsub(" ", "_", tr))
	phy$tip.label <- gsub("'", "", gsub("_", " ", phy$tip.label))
	return(	phy)
}

#' Get all chronograms from Open Tree of Life
#' @param verbose If TRUE, give updates to the user
#' @return A list with elements for the trees, authors, curators, and study ids
#' @export
GetOToLChronograms <- function(verbose=FALSE) {
	chronogram.matches <- rotl::studies_find_trees(property="ot:branchLengthMode", value="ot:time")
	trees <- list()
	authors <- list()
	curators <- list()
	studies <- list()
	tree.count <- 0
	for (study.index in sequence(dim(chronogram.matches)[1])) {
		if(verbose) {
			print(paste("Downloading tree(s) from study", study.index, "of",dim(chronogram.matches)[1]))
		}
		for(chrono.index in sequence((chronogram.matches$n_matched_trees[study.index]))) {
			study.id <- chronogram.matches$study_ids[study.index]
	#	new.tree <- get_study_tree(study_id=study.id, tree_id=tree.id, tip_label='ott_taxon_name')
			new.tree <- NULL
			try(new.tree <- get_study_tree_with_dups(study_id=study.id, tree_id=strsplit(chronogram.matches$match_tree_ids[study.index], ", ")[[1]][chrono.index]))
			if(!is.null(new.tree) & HasBrlen(new.tree)) {
				new.tree <- CleanChronogram(new.tree)
				if(HasBrlen(new.tree)) {
					if(IsGoodChronogram(new.tree)) {
						if(verbose) {
							print("has tree with branch lengths")
						}
						doi <- NULL
						try(doi <- gsub('http://dx.doi.org/', '', attr(rotl::get_publication(rotl::get_study_meta(study.id)), "DOI")))
						authors <- append(authors, NA)
						try(authors[length(authors)] <- list(paste(as.character(knitcitations::bib_metadata(doi)$author))))
						curators <- append(curators, NA)
						try(curators[length(curators)] <- list(rotl::get_study_meta(study.id)[["nexml"]][["^ot:curatorName"]]))
						try(studies <- append(studies, study.id))
						tree.count <- tree.count+1
						trees[[tree.count]] <-new.tree
						names(trees)[tree.count] <- rotl::get_publication(rotl::get_study_meta(study.id))[1]
					}
				}
			}
			#save(list=ls(), file="opentree_chronograms.RData")
		}
	}
	result <- list(trees=trees, authors=authors, curators=curators, studies=studies)
	return(result)
}

#' Save all chronograms from Open Tree of Life
#' @param file Path including file name
#' @param verbose If TRUE, give status updates to the user
#' @return None
#' @export
SaveOToLChronograms <- function(file="opentree_chronograms.rda", verbose=FALSE) {
	datelife.cache <- GetOToLChronograms(verbose=verbose)
	save(datelife.cache, file=file)
}

#' Check to see that a chronogram is valid
#' @param phy Input phylo object
#' @return Boolean: TRUE if good tree
#' @export
IsGoodChronogram <- function(phy) {
	passing <- TRUE
	if(class(phy) != "phylo") {
		passing <- FALSE
	}
	if(ape::Ntip(phy)<=ape::Nnode(phy)) {
		passing <- FALSE
	}
	if(length(which(grepl("not mapped", phy$tip.label)))>0) {
		passing <- FALSE #not cleaned properly
	}
	if(min(nchar(phy$tip.label))<=2) {
		passing <- FALSE
	}
	if(!ape::is.rooted(phy)) {
		passing <- FALSE
	}
	return(passing)
}

#' Clean up some issues with OToL chronograms
#' @param phy Input phylo object
#' @return A cleaned up phylo object
#' @export
CleanChronogram <- function(phy) {
	if(class(phy)=="phylo") {
		if(ape::Ntip(phy)>ape::Nnode(phy)) {
			bad.taxa <- c(which(nchar(phy$tip.label)<=2), which(grepl("not mapped", phy$tip.label)))
			if(length(bad.taxa)>0 & length(bad.taxa) < (ape::Ntip(phy)-3)) {
				phy <- try(ape::drop.tip(phy, bad.taxa))
			}
			if(!ape::is.rooted(phy) & ape::is.ultrametric(phy)) {
				phy$root.edge <- 0
			}
		}
	}
	return(phy)
}
