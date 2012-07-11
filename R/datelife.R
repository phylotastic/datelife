#skeleton of code to make open access info on dating life
library(ape)
library(phylobase)
library(abind)
library(compare)
library(phangorn)

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
  if (compare(rownames(patristic.matrix),standard.rownames)$result!=TRUE) {
    return(FALSE) 
  }
  if (compare(colnames(patristic.matrix),standard.colnames)$result!=TRUE) {
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
