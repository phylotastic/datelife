#skeleton of code to make open access info on dating life
library(ape)
library(phylobase)
library(abind)
library(phangorn)
library(compare)
library(geiger) #note this uses the next version of geiger, which has congruifier code
library(datelife2) #Klaus Schliep's code for pruning efficiently. Eventually, this will be moved into phangorn

#Note that originally trees were stored as patristic matrices. This was intended
#to make subsetting fast. The downside is large memory usage. Klaus Schliep wrote
#fast tree subsetting for phylo and multiphylo objects, so now trees are stored
#internally as objects of this type, but with the final output after pruning
#going through patristic matrices. Both sets of functions are maintained for
#now, but eventually the ones that take a patristic.matrix.array and then 
#subset it will be deprecated for those that first take a multiphylo or phylo
#object.

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
		if (!is.ultrametric(phy,tol)) {
			stop("currently we require that chronograms be ultrametric") # can pad them so that terminals all reach to present
		}
	}
	patristic.matrix<-cophenetic(phy)
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
    if (length(descendants(phy4, MRCA(phy4, taxa), type="tips")) > taxa) {
       problem <- "set of taxa not a clade, so this is probably an overestimate"
    }
  }
  return(list(patristic.matrix=patristic.matrix,problem=problem))
}

GetSubsetArrayBothFromPhylo <- function(reference.tree, taxa, phy=NULL, phy4=NULL) {
  if (is.null(phy)) {
    return(GetSubsetArrayFromPhylo(reference.tree=reference.tree, taxa=taxa, phy4=phy4)) 
  }
  else { #congruify
    return(GetSubsetArrayCongruifyFromPhylo(reference.tree=reference.tree, taxa=taxa, phy=phy)) 
  }

}


GetSubsetArrayFromPhylo <- function(reference.tree, taxa, phy4=NULL) {
  reference.tree<-pruneTrees(reference.tree, taxa) #pruneTrees is the new, fast fn from Klaus Schliep. Eventually will be in phangorn, currently in datelife2
  problem <- "none"
  final.size <- 0
  patristic.matrix.array <- NA
  if (!is.null(reference.tree)) {
  	if(class(reference.tree)=="phylo") {
  		final.size<-Ntip(reference.tree)
  	} else {
  		final.size<-Ntip(reference.tree[[1]]) #since Ntip only works on phylo objects
  	}
  }
  if (final.size < length(taxa)) {
    problem <- "missing some taxa on chronogram, so this is probably an underestimate" # fewer taxa on final matrix than we asked for
    if (final.size < 2 ) {
      problem <- "insufficient coverage" # we either have one species or zero. Not enough for an MRCA 
      patristic.matrix.array <- NA # to make sure no one uses the zero by mistake
    }
  }
  if (final.size >= 2) {
  	if(class(reference.tree)=="phylo") {
  		patristic.matrix.array <- BindMatrices(list(ComputePatristicDistance(reference.tree)))
  	} else {
  		patristic.matrix.array <- BindMatrices(lapply(reference.tree,ComputePatristicDistance))
  	}
  }
  if(!is.null(phy4)) {
    if (length(descendants(phy4, MRCA(phy4, taxa), type="tips")) > taxa) {
      problem <- "set of taxa not a clade, so this is probably an overestimate"
    }
  }
  return(list(patristic.matrix.array=patristic.matrix.array,problem=problem))
}

GetSubsetArrayCongruifyFromPhylo <- function(reference.tree, taxa, phy=NULL) {
  reference.tree<-pruneTrees(reference.tree, taxa) #pruneTrees is the new, fast fn from Klaus Schliep. Eventually will be in phangorn, currently in datelife2
  problem <- "none"
  final.size <- 0
  patristic.matrix.array <- NA
  if (!is.null(reference.tree)) {
  	if(class(reference.tree)=="phylo") {
  		final.size<-Ntip(reference.tree)
  	} else {
  		final.size<-Ntip(reference.tree[[1]]) #since Ntip only works on phylo objects
  	}
  }
  if (final.size < length(taxa)) {
    problem <- "missing some taxa on chronogram, so this is probably an underestimate" # fewer taxa on final matrix than we asked for
    if (final.size < 3 ) {
      problem <- "insufficient coverage" # we either have one species or zero. Not enough for an MRCA 
      patristic.matrix.array <- NA # to make sure no one uses the zero by mistake
      return(list(patristic.matrix.array=patristic.matrix.array,problem=problem))
      
    }
  }
  if (final.size >= 3) {
  	if(class(reference.tree)=="phylo") {
  		patristic.matrix.array <- BindMatrices(list(CongruifyTreeFromPhylo(reference.tree, query.tree=phy)))
  	} else {
  		patristic.matrix.array <- BindMatrices(lapply(reference.tree,CongruifyTreeFromPhylo, query.tree=phy)) #think about mclapply
  	}
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
  try(result.matrix<-ComputePatristicDistance(congruify.phylo(PatristicMatrixToTree(patristic.matrix), query.tree, NULL, 0, scale="PATHd8")$phy))
  return(result.matrix)
}

CongruifyTreeFromPhylo <- function(reference.tree, query.tree) {
  result.matrix<-matrix(nrow=Ntip(reference.tree), ncol=Ntip(reference.tree))
  if(is.null(query.tree$edge.length)) {
    query.tree$edge.length<-numeric(nrow(query.tree$edge)) #makes it so that branches that don't match reference tree get zero length
  }
  try(result.matrix<-ComputePatristicDistance(congruify.phylo(reference.tree, query.tree, NULL, 0, scale="PATHd8")$phy))
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
    if (length(descendants(phy4, MRCA(phy4, taxa), type="tips")) > taxa) {
      problem <- "set of taxa not a clade, so this is probably an overestimate"
    }
  }
  return(list(patristic.matrix.array=patristic.matrix.array,problem=problem))
}

GetAge <- function(patristic.matrix) {
  # 0.5 since patristic distance is down to the root and back up
  return(0.5 * max(patristic.matrix)) 
}

GetAges <- function(patristic.matrix.array) {
  if (length(dim(patristic.matrix.array))==2) {
    return(GetAge(patristic.matrix.array))
  }
  else if (dim(patristic.matrix.array)[3]==1) {
    return(GetAge(patristic.matrix.array))
  }
  else {
    return( sapply(SplitArray(patristic.matrix.array), GetAge ))
  }
}

ReorderMatrix <- function(patristic.matrix) {
  return(patristic.matrix[order(rownames(patristic.matrix)),order(colnames(patristic.matrix))]) 
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
  patristic.matrix.list<-lapply(patristic.matrix.list,ReorderMatrix)
  standard.rownames<-rownames(patristic.matrix.list[[1]])
  standard.colnames<-colnames(patristic.matrix.list[[1]])
  matching.names<-sapply(patristic.matrix.list,TestNameOrder,standard.rownames,standard.colnames)
  if (sum(matching.names)!=length(matching.names)) {
    stop("The patristic matrices you are trying to bind to not have the same taxa")
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

PatristicMatrixToTree <- function(patristic.matrix) {
  return(midpoint(nj(patristic.matrix)))
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

