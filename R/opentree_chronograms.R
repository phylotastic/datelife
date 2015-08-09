#install.packages("rotl")
library(rotl)
library(ape)
chronogram.matches <- studies_find_trees(property="ot:branchLengthMode", value="ot:time")[[1]]
trees <- list()

HasBrlen <- function(x) {
	brlen=TRUE
	if(is.null(x$edge.length)) {
		brlen=FALSE	
	}	
	return(brlen)
}

for (chrono.index in sequence(length(chronogram.matches))) {
	print(paste("tree number", chrono.index, "of", length(chronogram.matches)))
	study.id <- unlist(unname(chronogram.matches[[chrono.index]][1]))
	tree.id <- unname(unlist(chronogram.matches[[chrono.index]][[2]][[1]][2]))
	new.tree <- get_study_tree(study_id=study.id, tree_id=tree.id, tip_label=ott_taxon_name)
	print(new.tree)
	if(HasBrlen(new.tree)) {
		trees[[chrono.index]] <-new.tree
		names(trees)[chrono.index] <- get_publication(get_study_meta(study.id))[1]
	}
	save(list=ls(), file="opentree_chronograms.RData")
}

#since some of the trees come in with missing brlen, and so all brlen are cut, delete these trees



trees2 <- trees[sapply(trees, HasBrlen)]
save(list=ls(), file="opentree_chronograms.RData")

#print(get_study_tree(study_id="pg_2853", tree_id="tree6624"))
