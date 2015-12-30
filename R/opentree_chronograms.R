library(devtools)
#install_github("ropensci/rotl", dependencies = TRUE, build_vignette=FALSE)
library(rotl)
library(ape)
library(knitcitations)
chronogram.matches <- studies_find_trees(property="ot:branchLengthMode", value="ot:time")[[1]]
trees <- list()
authors <- list()
curators <- list()

HasBrlen <- function(x) {
	brlen=TRUE
	if(is.null(x$edge.length)) {
		brlen=FALSE	
	}	
	return(brlen)
}

#temporary hack to allow trees with duplicate taxa to be read in, courtesy David Winter
#note lack of underscores in ottaxonname in this fn
get_study_tree_with_dups <- function(study_id, tree_id, tip_label="ot:otttaxonname") {
	tr <- rotl:::.get_study_tree(study_id=study_id, tree_id=tree_id, tip_label=tip_label, format="newick")
	phy <- ape::read.tree(text=gsub(" ", "_", tr))
	phy$tip.label <- gsub("_", " ", phy$tip.label)
	return(	phy)
}



for (chrono.index in sequence(length(chronogram.matches))) {
	print(paste("tree number", chrono.index, "of", length(chronogram.matches)))
	study.id <- unlist(unname(chronogram.matches[[chrono.index]][1]))
	tree.id <- unname(unlist(chronogram.matches[[chrono.index]][[2]][[1]][2]))
#	new.tree <- get_study_tree(study_id=study.id, tree_id=tree.id, tip_label='ott_taxon_name')
	new.tree <- get_study_tree_with_dups(study_id=study.id, tree_id=tree.id)
	print(new.tree)
	if(HasBrlen(new.tree)) {
		doi <- NULL
		try(doi <- gsub('http://dx.doi.org/', '', attr(get_publication(get_study_meta(study.id)), "DOI")))
		try(authors <- append(authors, bib_metadata(doi)$author))
		try(curators <- unlist(append(curators, get_study_meta(study.id)[["nexml"]][["^ot:curatorName"]])))
		trees[[chrono.index]] <-new.tree
		names(trees)[chrono.index] <- get_publication(get_study_meta(study.id))[1]
	}
	save(list=ls(), file="opentree_chronograms.RData")
}

#since some of the trees come in with missing brlen, and so all brlen are cut, delete these trees

print(t(t(sort(table(as.character(authors)), decreasing=TRUE))))
print(t(t(sort(table(as.character(curators)), decreasing=TRUE))))


trees2 <- trees[sapply(trees, HasBrlen)]
save(list=ls(), file="opentree_chronograms.RData")

#print(get_study_tree(study_id="pg_2853", tree_id="tree6624"))
