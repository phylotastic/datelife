library(devtools)
#install_github("ropensci/rotl", dependencies = TRUE, build_vignette=FALSE)
library(rotl)
library(ape)
library(knitcitations)
chronogram.matches <- studies_find_trees(property="ot:branchLengthMode", value="ot:time")
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


tree.count <- 0
for (study.index in sequence(dim(chronogram.matches)[1])) {
	for(chrono.index in sequence(length(chronogram.matches$n_matched_trees[study.index]))) {
		study.id <- chronogram.matches$study_ids[study.index]
#	new.tree <- get_study_tree(study_id=study.id, tree_id=tree.id, tip_label='ott_taxon_name')
		new.tree <- get_study_tree_with_dups(study_id=study.id, tree_id=strsplit(chronogram.matches$match_tree_ids[study.index], ", ")[[1]][chrono.index])
		if(HasBrlen(new.tree)) {
			doi <- NULL
			try(doi <- gsub('http://dx.doi.org/', '', attr(get_publication(get_study_meta(study.id)), "DOI")))
			try(authors <- append(authors, bib_metadata(doi)$author))
			try(curators <- unlist(append(curators, get_study_meta(study.id)[["nexml"]][["^ot:curatorName"]])))
			tree.count <- tree.count+1
			trees[[tree.count]] <-new.tree
			names(trees)[tree.count] <- get_publication(get_study_meta(study.id))[1]
		}
		save(list=ls(), file="opentree_chronograms.RData")
	}
}

#since some of the trees come in with missing brlen, and so all brlen are cut, delete these trees

print(t(t(sort(table(as.character(authors)), decreasing=TRUE))))
print(t(t(sort(table(as.character(curators)), decreasing=TRUE))))


trees2 <- trees[sapply(trees, HasBrlen)]
save(list=ls(), file="opentree_chronograms.RData")

#print(get_study_tree(study_id="pg_2853", tree_id="tree6624"))
