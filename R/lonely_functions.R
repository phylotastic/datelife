GetQuantiles <- function(ages,probs = c(0.5,0,0.025,0.975,1) ) {
  # just utility wrapper function with different defaults
  return(stats::quantile(ages,probs))
}

VectorToTableRow <- function(x,digits = 2) {
  return(paste(paste("<td>",round(x,digits),sep = ""),"</td>",sep = "",collapse = ""))
}

SamplePatristicMatrix <- function(patristic.matrix.array, uncertainty) {
 # if (dim(patristic.matrix.array)[3] == 1) {
 # 	patristic.matrix<-patristic.matrix.array[,,1]
 # 	#need order of node depths, from just the upper triangular and diagonal part of the matrix
 # 	element.order<-order(patristic.matrix[upper.tri(patristic.matrix,diag = FALSE)],decreasing = TRUE)
 # 	new.patristic.matrix<-patristic.matrix*0
 # 	cur.val<-patristic.matrix[upper.tri(patristic.matrix,diag = FALSE)][element.order[1]]
 #   new.patristic.matrix[upper.tri(new.patristic.matrix,diag = FALSE)][element.order[1]] <- cur.val + runif(1, -cur.val*uncertainty/100, cur.val*uncertainty/100)
#	element.order<-element.order[-1]
#  	for (i in sequence(length(element.order))) {
#  		cur.val<-patristic.matrix[upper.tri(patristic.matrix,diag = FALSE)][element.order[i]]
#  		new.patristic.matrix[upper.tri(new.patristic.matrix,diag = FALSE)][element.order[i]] <- cur.val + runif(1, -cur.val*uncertainty/100, min(cur.val*uncertainty/100, min( ))
#  	}
#  }
#  else {
  	return(patristic.matrix <- patristic.matrix.array[,,sample.int(1, size = dim(patristic.matrix.array)[3] )] )
 # }
}

GetSubsetMatrix <- function(patristic.matrix, taxa, phy4 = NULL) {
  #gets a subset of the patristic.matrix. If you give it a phylo4 object, it can check to see if taxa are a clade
  patristic.matrix.new <- patristic.matrix[ rownames(patristic.matrix) %in% taxa,colnames(patristic.matrix) %in% taxa ]
  problem.new <- "none"
  final.size <- sum(rownames(patristic.matrix.new) %in% taxa) # returns number of matches
  if (final.size < length(taxa)) {
    problem.new <- "some of the queried taxa are not on this chronogram, so this is probably an underestimate" # fewer taxa on final matrix than we asked for
    if (final.size < 2 ) {
      problem.new <- "insufficient coverage" # we either have one species or zero. Not enough for an MRCA
      patristic.matrix.new <- NA # to make sure no one uses the zero by mistake
    }
  }
  if(!is.null(phy4)) {
    if (length(phylobase::descendants(phy4, phylobase::MRCA(phy4, taxa), type = "tips")) > taxa) {
       problem <- "set of taxa not a clade, so this is probably an overestimate"
    }
  }
  return(list(patristic.matrix= patristic.matrix.new,problem= problem.new))
}

#' Return the relevant authors for a set of studies
#' @param results.index A vector from FindMatchingStudyIndex() with the indices of the relevant studies
#' @param cache The cache
#' @return A vector with counts of each author, with names equal to author names
#' @export
TabulateRelevantAuthors <- function(results.index, cache = get("opentree_chronograms")) {
	authors <- cache$authors[results.index]
	return(table(unlist(authors)))
}

#' Return the relevant curators for a set of studies
#' @param results.index A vector from FindMatchingStudyIndex() with the indices of the relevant studies
#' @param cache The cache
#' @return A vector with counts of each curator, with names equal to curator names
#' @export
TabulateRelevantCurators <- function(results.index, cache = get("opentree_chronograms")) {
	curators <- cache$curators[results.index]
	return(table(unlist(curators)))
}
