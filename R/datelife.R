#skeleton of code to make open access info on dating life
# library(ape)
# library(abind)
# library(phangorn)
# library(compare)
# library(geiger) #note this uses the next version of geiger, which has congruifier code. The relevant code is in auteur-congruify.R
# library(datelife2) #Klaus Schliep's code for pruning efficiently. Eventually, this will be moved into phangorn
# library(parallel)
# library(doMC)
# library(stringr)
# library(ggplot2)
# library(taxize)
# library(plyr)
# library(rotl)
# source("/Library/WebServer/Sites/datelife.org/datelife/R/cleaning.r")

#' Take input string, figure out if it's newick or list of species
#' @param input A newick string or vector of taxa
#' @param usetnrs Whether to use OpenTree's TNRS for the input
#' @param approximatematch IDK
#' @return A list with the phy (or NA, if no tree) and cleaned vector of taxa
#' @export
ProcessInput <- function(input=c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"), usetnrs="no", approximatematch="yes") {
 input<-gsub("\\+"," ",input)
  input<-str_trim(input, side = "both")
  phy <- NA
  do_approximate_matching = TRUE
  if (approximatematch=="no") {
    do_approximate_matching = FALSE
  }
  prune_na = TRUE
  if (prunenonmatch =="no") {
    prune_na = FALSE
  }
  
  if(grepl('\\(', input) & grepl('\\)', input) & (substr(input,nchar(input),nchar(input))==";")) { #our test for newick
    phy<-read.tree(text=input)
  }
  cleaned.names<-""
  if(!is.null(phy)) {
    if(usetnrs=="yes") {
      phy <- tnrs_OToL_phylo(phylo=phy, do_approximate_matching, prune_na= prune_na)
    }
    	cleaned.names<-phy$tip.label
    } else {
      #cleaned.names<-strsplit( gsub("\\s","",input), ",")[[1]]
      cleaned.names <- input
      if (usetnrs=="yes") {
        phy <- tnrs_OToL_names(names=cleaned.names, do_approximate_matching, prune_na= prune_na)

      }
    }
    return(list(phy=phy, cleaned.names=cleaned.names))
}

AllMatching <- function(patristic.matrix, taxa) {
	return(sum(!(taxa %in% rownames(patristic.matrix) ))==0)
}


#' Find the index of relevant studies in a datelife.cache object
#' @param results.list A list returned from using GetSubsetArrayDispatch on datelife.cache$trees
#' @param taxa A vector of taxa to match
#' @param partial If TRUE, return matrices that have only partial matches
#' @return A vector with the indices of studies that have relevant info
#' @export
FindMatchingStudyIndex <- function(results.list, taxa=NULL, partial=FALSE) {
  results.list <- lapply(datelife.cache$trees,GetSubsetArrayDispatch, taxa=taxa, phy=NULL)
  patristic.matrices <- lapply(results.list, "[[", "patristic.matrix.array")
  results <- which(!is.na(patristic.matrices))
  if(!partial) {
  	results <- which(sapply(patristic.matrices, AllMatching, taxa=taxa))
  }
  return(results)
}

#' Return the relevant authors for a set of studies 
#' @param results.index A vector from FindMatchingStudyIndex() with the indices of the relevant studies
#' @param datelife.cache The cache
#' @return A vector with counts of each author, with names equal to author names
#' @export
TabulateRelevantAuthors <- function(results.index, datelife.cache) {
	authors <- datelife.cache$authors[results.index]
	return(table(unlist(authors)))
}

#' Return the relevant curators for a set of studies 
#' @param results.index A vector from FindMatchingStudyIndex() with the indices of the relevant studies
#' @param datelife.cache The cache
#' @return A vector with counts of each curator, with names equal to curator names
#' @export
TabulateRelevantCurators <- function(results.index, datelife.cache) {
	curators <- datelife.cache$curators[results.index]
	return(table(unlist(curators)))
}

#' Take results.list and process it
#' @param results.list A list returned from using GetSubsetArrayDispatch on datelife.cache$trees
#' @param taxa A vector of taxa to match
#' @param partial If TRUE, return matrices that have only partial matches
#' @return A list with the patristic.matrices that are not NA
#' @export
ProcessResultsList <- function(results.list, taxa=NULL, partial=FALSE) {
	if(is.null(taxa)) {
		taxa <- unique(unname(unlist(lapply(final.matrices, rownames))))
	}
	patristic.matrices <- lapply(results.list, "[[", "patristic.matrix.array")	
	
	final.matrices <- patristic.matrices[!is.na(patristic.matrices)]
	if(!partial) {
		final.matrices <- final.matrices[sapply(final.matrices, AllMatching, taxa=taxa)]
	}
	if(length(final.matrices)>0) {
		final.matrices <- lapply(final.matrices, PadMatrix, all.taxa=taxa)	
	}
	return(final.matrices)
}

#Note that originally trees were stored as patristic matrices. This was intended
#to make subsetting fast. The downside is large memory usage. Klaus Schliep wrote
#fast tree subsetting for phylo and multiphylo objects, so now trees are stored
#internally as objects of this type, but with the final output after pruning
#going through patristic matrices.

#Some trees are so large that they can't be stored as patristic distance matrices. For all others, 
#patristic matrices are better. For example, for the 20,000 HeathEtAl2012 trees of 35 taxa,
#getting a subset down to two taxa takes 0.0475 seconds just for the pruning, 0.0504 seconds
#for pruning and getting a subset, for a single tree (run times go up linearly with number of trees:
#pruning and converting 1000 trees takes 3 seconds). Subsetting 1000 trees from the patristic
#distance matrix takes just 0.0013 seconds.

ReadRDFTree <- function(identifier, format = "phylo") {
	# reads a tree in RDF format from the hypothetical data store
	tree<-read.file( identifier ) #do magic
	if (format == "phylo") {
		return(as.phylo(tree)) #assumes we have a converter
	}
}

#in case we want to cache. Not clear we do
ComputePatristicDistance <- function(phy, test=TRUE,tol=0.0001) {
	# stores the distance between taxa
	if (test) {
		if (!ape::is.ultrametric(phy,tol)) {
			stop("currently we require that chronograms be ultrametric") # can pad them so that terminals all reach to present
		}
	}
	patristic.matrix<-stats::cophenetic(phy)
	return(patristic.matrix)
}

GetSubsetMatrix <- function(patristic.matrix, taxa, phy4=NULL) {
  #gets a subset of the patristic.matrix. If you give it a phylo4 object, it can check to see if taxa are a clade
  patristic.matrix <- patristic.matrix[ rownames(patristic.matrix) %in% taxa,colnames(patristic.matrix) %in% taxa ]
  problem <- "none"
  final.size <- sum(rownames(patristic.matrix) %in% taxa) # returns number of matches
  if (final.size < length(taxa)) {
    problem <- "some of the queried taxa are not on this chronogram, so this is probably an underestimate" # fewer taxa on final matrix than we asked for
    if (final.size < 2 ) {
      problem <- "insufficient coverage" # we either have one species or zero. Not enough for an MRCA 
      patristic.matrix <- NA # to make sure no one uses the zero by mistake
    }
  }
  if(!is.null(phy4)) {
    if (length(phylobase::descendants(phy4, phylobase::MRCA(phy4, taxa), type="tips")) > taxa) {
       problem <- "set of taxa not a clade, so this is probably an overestimate"
    }
  }
  return(list(patristic.matrix=patristic.matrix,problem=problem))
}

#' Summarize a filtered results list in various ways
#' @param filtered.results A list of patristic matrices; labels correspond to citations
#' @param output.format The desired output format
#' @return Depends on output format
#' @export
SummarizeResults <- function(filtered.results, output.format, partial=TRUE) {
	if(!partial) {
		filtered.results <- filtered.results[which(!sapply(filtered.results, anyNA))]
	}
	output.format <- match.arg(output.format, choices=c("citations", "mrca", "newick"))
	if(output.format=="citations") {
		return(names(filtered.results))	
	}
	if(output.format=="mrca") {
		return(GetAges(filtered.results, partial=partial))	
	}
	if(output.format=="newick") {
		trees <- sapply(filtered.results, PatristicMatrixToNewick)
		return(trees[which(!is.na(trees))])
	}
	
}

#' Figure out which subset function to use
#' @param study.element The thing being passed in: an array or a phylo to serve as reference
#' @param taxa Vector of taxon names to get a subset for
#' @param phy A user tree to congruify in phylo format (ape)
#' @param phy4 A user tree to congruify in phylo4 format (phylobase)
#' @return A patristic matrix with for the taxa.
#' @export
GetSubsetArrayDispatch <- function(study.element, taxa, phy=NULL, phy4=NULL) {
  if(class(study.element)=="array") {
    return(GetSubsetArrayBoth(study.element, taxa, phy, phy4))
  } else {
    return(GetSubsetArrayBothFromPhylo(study.element, taxa, phy, phy4))
  }
}

GetSubsetArrayBothFromPhylo <- function(reference.tree, taxa, phy=NULL, phy4=NULL) {
#COMMENTING OUT: OpenTree gives single trees, let's just standardize on those
#  if (class(reference.tree)=="phylo") {
#    reference.tree<-c(reference.tree) #from here in, assumes multiphylo object, even if a single tree 
#  }
  if (is.null(phy)) {
    return(GetSubsetArrayFromPhylo(reference.tree=reference.tree, taxa=taxa, phy4=phy4)) 
  }
  else { #congruify
    return(GetSubsetArrayCongruifyFromPhylo(reference.tree=reference.tree, taxa=taxa, phy=phy)) 
  }

}

PruneTree <- function(phy, taxa) {
	return(ape::drop.tip(phy, tip=phy$tip.label[-(which(phy$tip.label %in% taxa))]))
}

GetSubsetArrayFromPhylo <- function(reference.tree, taxa, phy4=NULL) {
  final.size<-sum(reference.tree$tip.label %in% taxa) # returns number of matches
  if(final.size>=2) { #it's worth doing the pruning
    reference.tree<-PruneTree(reference.tree, taxa) 
    #reference.tree<-pruneTrees(reference.tree, taxa) #pruneTrees is the new, fast fn from Klaus Schliep. Eventually will be in phangorn, currently in datelife2
  }
  problem <- "none"
  patristic.matrix.array <- NA
  if (final.size < length(taxa)) {
    problem <- "missing some taxa on chronogram, so this is probably an underestimate" # fewer taxa on final matrix than we asked for
    if (final.size < 2 ) {
      problem <- "insufficient coverage" # we either have one species or zero. Not enough for an MRCA 
      patristic.matrix.array <- NA # to make sure no one uses the zero by mistake
    }
  }
  if (final.size >= 2) {
  	patristic.matrix.array <- ComputePatristicDistance(reference.tree)
  }
  if(!is.null(phy4)) {
    if (length(phylobase::descendants(phy4, phylobase::MRCA(phy4, taxa), type="tips")) > taxa) {
      problem <- "set of taxa not a clade, so this is probably an overestimate"
    }
  }
  return(list(patristic.matrix.array=patristic.matrix.array,problem=problem))
}

GetSubsetArrayCongruifyFromPhylo <- function(reference.tree, taxa, phy=NULL) {
  final.size<-sum(reference.tree$tip.label %in% taxa) # returns number of matches
  if(final.size>=2) { #it's worth doing the pruning
   reference.tree<-PruneTree(reference.tree, taxa) 
   #reference.tree<-pruneTrees(reference.tree, taxa) #pruneTrees is the new, fast fn from Klaus Schliep. Eventually will be in phangorn, currently in datelife2
  }
  problem <- "none"
  patristic.matrix.array <- NA
  if (final.size < length(taxa)) {
    problem <- "missing some taxa on chronogram, so this is probably an underestimate" # fewer taxa on final matrix than we asked for
    if (final.size < 3 ) {
      problem <- "insufficient coverage" # we either have one species or zero. Not enough for an MRCA 
      patristic.matrix.array <- NA # to make sure no one uses the zero by mistake
      return(list(patristic.matrix.array=patristic.matrix.array,problem=problem))
      
    }
  }
  if (final.size >= 3) {
  	patristic.matrix.array <-   CongruifyTreeFromPhylo(reference.tree, query.tree=phy)
  }
  return(list(patristic.matrix.array=patristic.matrix.array,problem=problem))
}


GetSubsetArrayBoth <- function(patristic.matrix.array, taxa, phy=NULL, phy4=NULL) {
  if (is.null(phy)) {
    return(GetSubsetArray(patristic.matrix.array=patristic.matrix.array, taxa=taxa, phy4=phy4)) 
  }
  else { #congruify
    return(GetSubsetArrayCongruify(patristic.matrix.array=patristic.matrix.array, taxa=taxa, phy=phy)) 
  }
}

CongruifyTree <- function(patristic.matrix, query.tree) {
  result.matrix<-matrix(nrow=dim(patristic.matrix)[1], ncol=dim(patristic.matrix)[2])
  if(is.null(query.tree$edge.length)) {
    query.tree$edge.length<-numeric(nrow(query.tree$edge))
  }
  try(result.matrix<-ComputePatristicDistance(geiger::congruify.phylo(PatristicMatrixToTree(patristic.matrix), query.tree, NULL, 0, scale="PATHd8")$phy))
  return(result.matrix)
}

CongruifyTreeFromPhylo <- function(reference.tree, query.tree) {
  result.matrix<-matrix(nrow=ape::Ntip(reference.tree), ncol=ape::Ntip(reference.tree))
  if(is.null(query.tree$edge.length)) {
    query.tree$edge.length<-numeric(nrow(query.tree$edge)) #makes it so that branches that don't match reference tree get zero length
  }
  try(result.matrix<-ComputePatristicDistance(geiger::congruify.phylo(reference.tree, query.tree, NULL, 0, scale="PATHd8")$phy))
  return(result.matrix)
}

GetSubsetArrayCongruify <- function(patristic.matrix.array, taxa, phy=NULL) {
  #gets a subset of the patristic.matrix.array. 
  patristic.matrix.array <- patristic.matrix.array[ rownames(patristic.matrix.array) %in% taxa,colnames(patristic.matrix.array) %in% taxa,  ]
  problem <- "none"
  final.size <- sum(rownames(patristic.matrix.array) %in% taxa) # returns number of matches
  if (final.size < length(taxa)) {
    problem <- "missing some taxa on chronogram, so this is probably an underestimate" # fewer taxa on final matrix than we asked for
    if (final.size < 3 ) {
      problem <- "insufficient coverage" # we either have one species or zero. Not enough for an MRCA 
      patristic.matrix.array <- NA # to make sure no one uses the zero by mistake
      return(list(patristic.matrix.array=patristic.matrix.array,problem=problem))
      
    }
  }
  patristic.matrix.list <- SplitArray(patristic.matrix.array)
  patristic.matrix.array<-BindMatrices(lapply(patristic.matrix.list, CongruifyTree, query.tree=phy)) #yes, this should be parallel
  return(list(patristic.matrix.array=patristic.matrix.array,problem=problem))
}

GetSubsetArray <- function(patristic.matrix.array, taxa, phy4=NULL) {
  #gets a subset of the patristic.matrix.array. If you give it a phylo4 object, it can check to see if taxa are a clade
  patristic.matrix.array <- patristic.matrix.array[ rownames(patristic.matrix.array) %in% taxa,colnames(patristic.matrix.array) %in% taxa,  ]
  problem <- "none"
  final.size <- sum(rownames(patristic.matrix.array) %in% taxa) # returns number of matches
  if (final.size < length(taxa)) {
    problem <- "missing some taxa on chronogram, so this is probably an underestimate" # fewer taxa on final matrix than we asked for
    if (final.size < 2 ) {
      problem <- "insufficient coverage" # we either have one species or zero. Not enough for an MRCA 
      patristic.matrix.array <- NA # to make sure no one uses the zero by mistake
    }
  }
  if(!is.null(phy4)) {
    if (length(phylobase::descendants(phy4, phylobase::MRCA(phy4, taxa), type="tips")) > taxa) {
      problem <- "set of taxa not a clade, so this is probably an overestimate"
    }
  }
  return(list(patristic.matrix.array=patristic.matrix.array,problem=problem))
}

#' Get time of MRCA from patristic matrix
#' @param patristic.matrix A patristic matrix
#' @param partial If TRUE, drop NA from the patristic matrix; if FALSE, will return NA if there are missing entries
#' @return The depth of the MRCA
#' @export
GetAge <- function(patristic.matrix, partial=TRUE) {
  # 0.5 since patristic distance is down to the root and back up
  return(0.5 * max(patristic.matrix, na.rm=partial)) 
}

#' Get vector of MRCA from list of patristic matrices
#' @param filtered.results List of patristic matrices
#' @param partial If TRUE, drop NA from the patristic matrix; if FALSE, will return NA if there are missing entries
#' @return Vector of MRCA ages with names same as in filtered.results
#' @export
GetAges <- function(filtered.results, partial=TRUE) {
	ages <- sapply(filtered.results, GetAge, partial=partial)
	return(ages)
}

#' Function to reorder a matrix so that row and column labels are in alphabetical order
#' @param patristic.matrix A patristic matrix with row and column names for taxa
#' @return patristic.matrix A patristic matrix with row and column names for taxa in alphabetial order
#' @export
ReorderMatrix <- function(patristic.matrix) {
  return(patristic.matrix[order(rownames(patristic.matrix)),order(colnames(patristic.matrix))]) 
}

#' Function to fill in empty cells in a patristic matrix for missing taxa
#' @param patristic.matrix A patristic matrix with row and column names for taxa
#' @param all.taxa A vector of the names of all taxa you want, including ones not in the patristic matrix
#' @return Patristic.matrix for all.taxa, with NA for entries between taxa where at least one was not in the original patristic.matrix
#' @export
PadMatrix <- function(patristic.matrix, all.taxa) {
	final.matrix <- rbind(patristic.matrix, rep(NA, length(all.taxa) - dim(patristic.matrix)[1]))
	final.matrix <- cbind(final.matrix, rep(NA, length(all.taxa) - dim(patristic.matrix)[1]))
 	rownames(final.matrix) <- c(rownames(patristic.matrix), all.taxa[-which(all.taxa %in% rownames(patristic.matrix))])
 	colnames(final.matrix) <- c(colnames(patristic.matrix), all.taxa[-which(all.taxa %in% colnames(patristic.matrix))])
	return(ReorderMatrix(final.matrix))	
}

TestNameOrder <- function(patristic.matrix, standard.rownames, standard.colnames) {
  if (compare::compare(rownames(patristic.matrix),standard.rownames)$result!=TRUE) {
    return(FALSE) 
  }
  if (compare::compare(colnames(patristic.matrix),standard.colnames)$result!=TRUE) {
    return(FALSE)
  }
  return(TRUE)
}

BindMatrices <- function(patristic.matrix.list) {
  original.size<-length(patristic.matrix.list)
  patristic.matrix.list<-lapply(patristic.matrix.list,ReorderMatrix)
  if(length(patristic.matrix.list)<1) {
    stop(paste("The patristic matrices you are trying to bind are too few; input was ", original.size, " and current length is ", length(patristic.matrix.list), sep=""))
  }
  standard.rownames<-rownames(patristic.matrix.list[[1]])
  standard.colnames<-colnames(patristic.matrix.list[[1]])
  matching.names<-sapply(patristic.matrix.list,TestNameOrder,standard.rownames,standard.colnames)
  if (sum(matching.names)!=length(matching.names)) {
    stop("The patristic matrices you are trying to bind do not have the same taxa")
  }
  return(abind(patristic.matrix.list, along=3 ))
}

SplitArray <- function(patristic.matrix.array) {
  return(lapply(sequence(dim(patristic.matrix.array)[3]),AsubForLapply,patristic.matrix.array))
}

AsubForLapply <- function(idx, x, dims=3) {
  return(asub(x, idx, dims)) 
}

GetQuantiles <- function(ages,probs=c(0.5,0,0.025,0.975,1) ) {
  # just utility wrapper function with different defaults
  return(quantile(ages,probs))
}

VectorToTableRow <- function(x,digits=2) {
  return(paste(paste("<td>",round(x,digits),sep=""),"</td>",sep="",collapse=""))
}

#' Convert patristic matrix to a phylo object
#' @param patristic.matrix A patristic matrix
#' @return A rooted phylo object
#' @export
PatristicMatrixToTree <- function(patristic.matrix) {
  if(anyNA(patristic.matrix)) {
  	patristic.matrix <- patristic.matrix[rowSums(is.na(patristic.matrix)) != ncol(patristic.matrix),colSums(is.na(patristic.matrix)) != nrow(patristic.matrix)]	
  }
  if(dim(patristic.matrix)[1] < 3) {
  	return(NA)	
  }
  tree <- ape::nj(patristic.matrix)
  if(ape::Ntip(tree)>2) {
    tree <- 	phangorn::midpoint(tree)
  }
  return(tree)
}

#' Convert patristic matrix to a newick string
#' @param patristic.matrix A patristic matrix
#' @return A newick string
#' @export
PatristicMatrixToNewick <- function(patristic.matrix) {
  tree <- PatristicMatrixToTree(patristic.matrix)
  if(class(tree)=="phylo") {
  	return(ape::write.tree(tree))
  }
  return(NA)
}



SummaryPatristicMatrix <- function(patristic.matrix.array,fn=median) {
  return(apply(patristic.matrix.array,MARGIN=c(1,2),fn))
}

SamplePatristicMatrix <- function(patristic.matrix.array, uncertainty) {
 # if (dim(patristic.matrix.array)[3] == 1) {
 # 	patristic.matrix<-patristic.matrix.array[,,1]
 # 	#need order of node depths, from just the upper triangular and diagonal part of the matrix
 # 	element.order<-order(patristic.matrix[upper.tri(patristic.matrix,diag=FALSE)],decreasing=TRUE)
 # 	new.patristic.matrix<-patristic.matrix*0
 # 	cur.val<-patristic.matrix[upper.tri(patristic.matrix,diag=FALSE)][element.order[1]]
 #   new.patristic.matrix[upper.tri(new.patristic.matrix,diag=FALSE)][element.order[1]] <- cur.val + runif(1, -cur.val*uncertainty/100, cur.val*uncertainty/100)
#	element.order<-element.order[-1]
#  	for (i in sequence(length(element.order))) {
#  		cur.val<-patristic.matrix[upper.tri(patristic.matrix,diag=FALSE)][element.order[i]]
#  		new.patristic.matrix[upper.tri(new.patristic.matrix,diag=FALSE)][element.order[i]] <- cur.val + runif(1, -cur.val*uncertainty/100, min(cur.val*uncertainty/100, min( ))
#  	}
#  }
#  else {
  	return(patristic.matrix<-patristic.matrix.array[,,sample.int(1, size=dim(patristic.matrix.array)[3] )] )
 # }
}

