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

#' Core function to go from a vector of species, newick string, or phylo object get a chronogram or dates back
#' @aliases datelife
#' @param input A vector of names, a newick string, or a phylo object
#' @param output.format The desired output format. See details.
#' @param partial If TRUE, use source trees even if they only match some of the desired taxa
#' @param usetnrs If TRUE, use OpenTree's services to resolve names. This can dramatically improve the chance of matches, but also take much longer
#' @param approximatematch If TRUE, use a slower TNRS to correct mispellings, increasing the chance of matches (including false matches)
#' @param cache The cached set of chronograms and other info from data(opentree_chronograms)
#' @param method The method used for congruification. PATHd8 only right now, r8s and treePL later.
#' @return Varies depending on the chosen format
#' @export
#' @details
#' The output formats are citations, mrca, newick.all, newick.sdm, newick.median, phylo.sdm, phylo.median, phylo.all, html
#' @examples
#' ages <- EstimateDates(c("Rhea americana", "Pterocnemia pennata", "Struthio camelus", "Mus musculus"), output.format="mrca")
EstimateDates <- function(input=c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"), output.format="phylo.sdm", partial=TRUE, usetnrs=FALSE, approximatematch=TRUE, cache=get("opentree_chronograms"), method="PATHd8") {
	filtered.results.in <- GetFilteredResults(input, partial, usetnrs, approximatematch, cache)
	output.format.in <- output.format
	cache.in <- cache
	return(SummarizeResults(filtered.results=filtered.results.in, output.format=output.format.in, cache=cache.in))
}

#' Go from a vector of species, newick string, or phylo object to a list of patristic matrices
#' @param input A vector of names, a newick string, or a phylo object
#' @param partial If TRUE, use source trees even if they only match some of the desired taxa
#' @param usetnrs If TRUE, use OpenTree's services to resolve names. This can dramatically improve the chance of matches, but also take much longer
#' @param approximatematch If TRUE, use a slower TNRS to correct mispellings, increasing the chance of matches (including false matches)
#' @param cache The cached set of chronograms and other info from data(opentree_chronograms)
#' @param method The method used for congruification. PATHd8 only right now, r8s and treePL later.
#' @return List of patristic matrices
#' @export
GetFilteredResults <- function(input=c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"), partial=TRUE, usetnrs=FALSE, approximatematch=TRUE, cache=get("opentree_chronograms"), method="PATHd8") {
    input.processed <- ProcessInput(input, usetnrs, approximatematch)
    tree <- input.processed$phy
    cleaned.names <- input.processed$cleaned.names
    results.list<-lapply(cache$trees,GetSubsetArrayDispatch, taxa=cleaned.names, phy=tree, method=method)
    filtered.results <- ProcessResultsList(results.list, cleaned.names, partial)
    return(filtered.results)
}

#' Take input string, figure out if it's newick or list of species
#' @param input A newick string or vector of taxa, or a phylo object
#' @param usetnrs Whether to use OpenTree's TNRS for the input
#' @param approximatematch If using TNRS, use approximate matching
#' @return A list with the phy (or NA, if no tree) and cleaned vector of taxa
#' @export
ProcessInput <- function(input=c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"), usetnrs=FALSE, approximatematch=TRUE) {
	if(class(input) == "phylo") {
		input <- ape::write.tree(input)
	}
 input<-gsub("\\+"," ",input)
  input<-stringr::str_trim(input, side = "both")
  phy.new <- NA
   if(length(input)==1) {
	  if(grepl('\\(', input) & grepl('\\)', input) & (substr(input,nchar(input),nchar(input))==";")) { #our test for newick
	    phy.new <-ape::read.tree(text=gsub(" ", "_", input))
	  }
  }
  cleaned.names<-""
  if(!is.na(phy.new[1])) {
    if(usetnrs) {
			cleaned.names <- rotl::tnrs_match_names(phy.new$tip.label)$unique_name
      phy.new$tip.label <- gsub("_", " ", cleaned.names)
    }
  	cleaned.names<-gsub("_", " ", phy.new$tip.label)
  } else {
    #cleaned.names<-strsplit( gsub("\\s","",input), ",")[[1]]
		if(length(input)==1) {
			input <- strsplit(input, ',')[[1]]
		}
		cleaned.names <- stringr::str_trim(input, side = "both")
    #cleaned.names <- input
    if (usetnrs) {
      cleaned.names <- gsub("_", " ", rotl::tnrs_match_names(cleaned.names)$unique_name)
    }
  }
  cleaned.names <- gsub("_", " ", cleaned.names)
  return(list(phy=phy.new, cleaned.names=cleaned.names))
}

#' Are all desired taxa in the patristic.matrix?
#' @param patristic.matrix A patristic matrix, rownames and colnames must be taxa
#' @param taxa A vector of taxon names
#' @return A Boolean
#' @export
AllMatching <- function(patristic.matrix, taxa) {
	return(sum(!(taxa %in% rownames(patristic.matrix) ))==0)
}


#' Find the index of relevant studies in a opentree_chronograms object
#' @param filtered.results The patristic.matrices that will be used
#' @param cache The cache of studies
#' @return A vector with the indices of studies that have relevant info
#' @export
FindMatchingStudyIndex <- function(filtered.results, cache=get("opentree_chronograms")) {
    return(which(names(cache$trees) %in% names(filtered.results)))
}

#' Return the relevant authors for a set of studies
#' @param results.index A vector from FindMatchingStudyIndex() with the indices of the relevant studies
#' @param cache The cache
#' @return A vector with counts of each author, with names equal to author names
#' @export
TabulateRelevantAuthors <- function(results.index, cache=get("opentree_chronograms")) {
	authors <- cache$authors[results.index]
	return(table(unlist(authors)))
}

#' Return the relevant curators for a set of studies
#' @param results.index A vector from FindMatchingStudyIndex() with the indices of the relevant studies
#' @param cache The cache
#' @return A vector with counts of each curator, with names equal to curator names
#' @export
TabulateRelevantCurators <- function(results.index, cache=get("opentree_chronograms")) {
	curators <- cache$curators[results.index]
	return(table(unlist(curators)))
}

#' Take results.list and process it
#' @param results.list A list returned from using GetSubsetArrayDispatch on opentree_chronograms$trees
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
		to.delete <- c()
		for (i in sequence(length(final.matrices))) {
			if(all(is.na(final.matrices[[i]]))) {
				to.delete <- c(to.delete, i)
			}
		}
		if(length(to.delete)>0) {
			final.matrices <- final.matrices[-to.delete]
		}
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


#in case we want to cache. Not clear we do
ComputePatristicDistance <- function(phy, test=TRUE,tol=0.01) {
	# stores the distance between taxa
	patristic.matrix <- NA
	if(class(phy)=="phylo") {
		if (test) {
			if (!ape::is.ultrametric(phy,tol)) {
				stop("currently we require that chronograms be ultrametric") # can pad them so that terminals all reach to present
			}
		}
		patristic.matrix<-stats::cophenetic(phy)
	}
	return(patristic.matrix)
}

GetSubsetMatrix <- function(patristic.matrix, taxa, phy4=NULL) {
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
    if (length(phylobase::descendants(phy4, phylobase::MRCA(phy4, taxa), type="tips")) > taxa) {
       problem <- "set of taxa not a clade, so this is probably an overestimate"
    }
  }
  return(list(patristic.matrix= patristic.matrix.new,problem= problem.new))
}

#' Summarize a filtered results list in various ways
#' @param filtered.results A list of patristic matrices; labels correspond to citations
#' @param output.format The desired output format
#' @param partial If TRUE, use source trees even if they only match some of the desired taxa
#' @param cache The cached trees and other information
#' @param suppress.citations If using a format that would normally print() citations, turn this off
#' @return Depends on output format
#' @export
SummarizeResults <- function(filtered.results, output.format, partial=TRUE, cache=get("opentree_chronograms"), suppress.citations=FALSE) {
	if(!partial) {
		filtered.results <- filtered.results[which(!sapply(filtered.results, anyNA))]
	}
	results.index <- FindMatchingStudyIndex(filtered.results, cache)
	output.format <- match.arg(output.format, choices=c("citations", "mrca", "newick.all", "newick.sdm", "newick.median", "phylo.sdm", "phylo.median", "phylo.median", "phylo.all", "html", "data.frame"))
	return.object <- NA
	if(output.format=="citations") {
		return.object <- names(filtered.results)
	}
	if(output.format=="mrca") {
		return.object <- GetAges(filtered.results, partial=partial)
	}
	if(output.format=="newick.all") {
		trees <- sapply(filtered.results, PatristicMatrixToNewick)
		return.object <- trees[which(!is.na(trees))]
	}
	if(output.format=="newick.sdm") {
		local.results <- RunSDM(filtered.results)
		filtered.results <- local.results$filtered.results
		tree <- local.results$phy
		return.object <- ape::write.tree(tree)
	}
	if(output.format=="phylo.sdm") {
		local.results <- RunSDM(filtered.results)
		filtered.results <- local.results$filtered.results
		tree <- local.results$phy
		return.object <- tree
	}
	if(output.format=="newick.median") {
		patristic.array <- BindMatrices(filtered.results)
		median.matrix <- SummaryPatristicMatrixArray(patristic.array)
		tree <- PatristicMatrixToNewick(median.matrix)
		return.object <- tree
	}
	if(output.format=="phylo.median") {
		patristic.array <- BindMatrices(filtered.results)
		median.matrix <- SummaryPatristicMatrixArray(patristic.array)
		tree <- PatristicMatrixToTree(median.matrix)
		return.object <- tree
	}
	if(output.format=="phylo.all") {
		trees <- lapply(filtered.results, PatristicMatrixToTree)
		return.object <- trees[which(!is.na(trees))]
	}
	if(output.format=="html") {
		out.vector <- "<table border='1'><tr><th>MRCA Age (MY)</th><th>Ntax</th><th>Citation</th><th>Newick</th></tr>"
		ages <- GetAges(filtered.results, partial=partial)
		trees <- sapply(filtered.results, PatristicMatrixToNewick)
		for(result.index in sequence(length(filtered.results))) {
			out.vector <- paste(out.vector, paste("<tr><td>",ages[result.index],"</td><td>",sum(!is.na(diag(filtered.results[[result.index]]))), "</td><td>", names(filtered.results)[result.index], "</td><td>", trees[result.index], "</td></tr>", sep=""), sep="")
		}
		out.vector <- paste(out.vector, "</table>")
		return.object <- out.vector
	}
	if(output.format=="data.frame") {
		out.df <- data.frame()
		ages <- GetAges(filtered.results, partial=partial)
		trees <- sapply(filtered.results, PatristicMatrixToNewick)
		for(result.index in sequence(length(filtered.results))) {
			out.line<- data.frame(Age=ages[result.index],Ntax=sum(!is.na(diag(filtered.results[[result.index]]))), Citation=names(filtered.results)[result.index], Newick= trees[result.index])
			if(result.index==1) {
				out.df <- out.line
			} else {
				out.df <- rbind(out.df, out.line)
			}
		}
		rownames(out.df) <- NULL
		return.object <- out.df
	}
	if((output.format != "citations") & !suppress.citations) {
		print("Using trees from:")
		print(names(filtered.results))
	}
	return(return.object)
}

#' Figure out which subset function to use
#' @param study.element The thing being passed in: an array or a phylo to serve as reference
#' @param taxa Vector of taxon names to get a subset for
#' @param phy A user tree to congruify in phylo format (ape)
#' @param phy4 A user tree to congruify in phylo4 format (phylobase)
#' @param method Which method to use for congruification
#' @return A patristic matrix with for the taxa.
#' @export
GetSubsetArrayDispatch <- function(study.element, taxa, phy=NULL, phy4=NULL, method="PATHd8") {
  if(class(study.element)=="array") {
    return(GetSubsetArrayBoth(study.element, taxa, phy, phy4, method))
  } else {
    return(GetSubsetArrayBothFromPhylo(study.element, taxa, phy, phy4, method))
  }
}

GetSubsetArrayBothFromPhylo <- function(reference.tree.in, taxa.in, phy.in=NULL, phy4.in=NULL, method.in="PATHd8") {
#COMMENTING OUT: OpenTree gives single trees, let's just standardize on those
#  if (class(reference.tree)=="phylo") {
#    reference.tree<-c(reference.tree) #from here in, assumes multiphylo object, even if a single tree
#  }
	congruify=FALSE
	if(!is.null(phy.in[1])) {
		congruify=TRUE
		if(is.na(phy.in[1])) {
			congruify=FALSE
		}
	}
  if (!congruify) {
    return(GetSubsetArrayFromPhylo(reference.tree=reference.tree.in, taxa=taxa.in, phy4=phy4.in, method.in))
  }
  else { #congruify
    return(GetSubsetArrayCongruifyFromPhylo(reference.tree=reference.tree.in, taxa=taxa.in, phy=phy.in, method.in))
  }

}

PruneTree <- function(phy, taxa) {
	return(ape::drop.tip(phy, tip=phy$tip.label[-(which(phy$tip.label %in% taxa))]))
}

GetSubsetArrayFromPhylo <- function(reference.tree, taxa, phy4=NULL, method="PATHd8") {
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

GetSubsetArrayCongruifyFromPhylo <- function(reference.tree, taxa, phy=NULL, method="PATHd8") {
  final.size<-sum(reference.tree$tip.label %in% taxa) # returns number of matches
  if(final.size>=2) { #it's worth doing the pruning
   reference.tree<-PruneTree(reference.tree, taxa)
   #reference.tree<-pruneTrees(reference.tree, taxa) #pruneTrees is the new, fast fn from Klaus Schliep. Eventually will be in phangorn, currently in datelife2
  }
  problem.new <- "none"
  patristic.matrix.array.new <- NA
  if (final.size < length(taxa)) {
    problem.new <- "missing some taxa on chronogram, so this is probably an underestimate" # fewer taxa on final matrix than we asked for
    if (final.size < 3 ) {
      problem.new <- "insufficient coverage" # we either have one species or zero. Not enough for an MRCA
      patristic.matrix.array.new <- NA # to make sure no one uses the zero by mistake
      return(list(patristic.matrix.array= patristic.matrix.array.new,problem= problem.new))

    }
  }
  if (final.size >= 3) {
  	patristic.matrix.array.new <-   CongruifyTreeFromPhylo(reference.tree, query.tree=phy)
  }
  return(list(patristic.matrix.array= patristic.matrix.array.new,problem= problem.new))
}


GetSubsetArrayBoth <- function(patristic.matrix.array, taxa, phy=NULL, phy4=NULL, method="PATHd8") {
  if (is.null(phy)) {
    return(GetSubsetArray(patristic.matrix.array=patristic.matrix.array, taxa=taxa, phy4=phy4))
  }
  else { #congruify
    return(GetSubsetArrayCongruify(patristic.matrix.array=patristic.matrix.array, taxa=taxa, phy=phy, method))
  }
}

CongruifyTree <- function(patristic.matrix, query.tree, method="PATHd8", attempt.fix=TRUE) {
  result.matrix<-matrix(nrow=dim(patristic.matrix)[1], ncol=dim(patristic.matrix)[2])
  if(is.null(query.tree$edge.length)) {
    query.tree$edge.length<-numeric(nrow(query.tree$edge))
  }
#	try(result.matrix<-ComputePatristicDistance(ConvertUnderscoresToSpaces(geiger::congruify.phylo(ConvertSpacesToUnderscores(PatristicMatrixToTree(patristic.matrix)), ConvertSpacesToUnderscores(query.tree), NULL, 0, scale=method)$phy)))
try(result.matrix<-ComputePatristicDistance(CongruifyAndCheck(reference=PatristicMatrixToTree(patristic.matrix), target=query.tree, scale=method, attempt.fix=attempt.fix)))
  return(result.matrix)
}

CongruifyTreeFromPhylo <- function(reference.tree, query.tree, method="PATHd8", attempt.fix=TRUE) {
  result.matrix<-matrix(nrow=ape::Ntip(reference.tree), ncol=ape::Ntip(reference.tree))
  if(is.null(query.tree$edge.length)) {
    query.tree$edge.length<-numeric(nrow(query.tree$edge)) #makes it so that branches that don't match reference tree get zero length
  }
	try(result.matrix<-ComputePatristicDistance(CongruifyAndCheck(reference=reference.tree, target=query.tree, scale=method, attempt.fix=attempt.fix)))
  return(result.matrix)
}

CongruifyAndCheck <- function(reference, target, taxonomy=NULL, tol=0.01, scale="pathd8", attempt.fix=TRUE) {
  if(!ape::is.ultrametric(reference, tol=tol)) {
    return(NA)
  }
	new.tree <- ConvertUnderscoresToSpaces(suppressWarnings(geiger::congruify.phylo(ConvertSpacesToUnderscores(reference), ConvertSpacesToUnderscores(target), taxonomy=taxonomy, tol=tol, scale=scale)$phy)) #suppressing warnings b/c geiger ignores tolerance
	if(anyNA(new.tree$edge.length) & attempt.fix) {
		warning("Congruification resulted in NA edge lengths. Resolving polytomies and making up starting branch lengths")
		new.tree <- ConvertUnderscoresToSpaces(geiger::congruify.phylo(ConvertSpacesToUnderscores(reference), ConvertSpacesToUnderscores(ape::compute.brlen(ape::multi2di(target))), taxonomy, tol, scale)$phy)
		if(anyNA(new.tree$edge.length)) {
			new.tree <- NA
		}
	} 
	return(new.tree)
}

#' Convert spaces to underscores in trees
#' @param phy A phylo object
#' @return A phylo object
#' @export
ConvertSpacesToUnderscores <- function(phy) {
	phy$tip.label <- gsub(" ", "_", phy$tip.label)
	return(phy)
}

#' Convert underscores to spaces in trees
#' @param phy A phylo object
#' @return A phylo object
#' @export
ConvertUnderscoresToSpaces <- function(phy) {
	phy$tip.label <- gsub("_", " ", phy$tip.label)
	return(phy)
}

GetSubsetArrayCongruify <- function(patristic.matrix.array, taxa, phy=NULL, method="PATHd8") {
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
  patristic.matrix.array<-BindMatrices(lapply(patristic.matrix.list, CongruifyTree, query.tree=phy, scale=method)) #yes, this should be parallel
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
	number.missing <- length(all.taxa) - dim(patristic.matrix)[1]
	final.matrix <- patristic.matrix
	if(number.missing>0) {
		final.matrix <- rbind(patristic.matrix, matrix(nrow=number.missing, ncol=dim(patristic.matrix)[2]))
		final.matrix <- cbind(final.matrix, matrix(ncol=number.missing, nrow=dim(final.matrix)[1]))
		rownames(final.matrix) <- c(rownames(patristic.matrix), all.taxa[-which(all.taxa %in% rownames(patristic.matrix))])
	 	colnames(final.matrix) <- c(colnames(patristic.matrix), all.taxa[-which(all.taxa %in% colnames(patristic.matrix))])
	}
 	return(ReorderMatrix(final.matrix))
}

#' Function to remove missing taxa
#' @param patristic.matrix A patristic matrix with row and column names for taxa
#' @return Patristic.matrix for all.taxa
#' @export
UnpadMatrix <- function(patristic.matrix) {
	bad.ones <- which(apply(is.na(patristic.matrix),2,all))
	if(length(bad.ones)>0) {
		patristic.matrix <- patristic.matrix[-bad.ones, -bad.ones]
	}
	return(patristic.matrix)
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

#' Convert list of patristic matrices to a 3D array
#' @param patristic.matrix.list List of patristic matrices
#' @return A 3d array of patristic matrices
#' @export
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
  return(abind::abind(patristic.matrix.list, along=3 ))
}

SplitArray <- function(patristic.matrix.array) {
  return(lapply(sequence(dim(patristic.matrix.array)[3]),AsubForLapply,patristic.matrix.array))
}

AsubForLapply <- function(idx, x, dims=3) {
  return(abind::asub(x, idx, dims))
}

GetQuantiles <- function(ages,probs=c(0.5,0,0.025,0.975,1) ) {
  # just utility wrapper function with different defaults
  return(stats::quantile(ages,probs))
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


#' Summarize patristic matrix array (by default, median)
#' @param patristic.matrix.array 3D array of patristic matrices
#' @param fn The function to use ot summarize
#' @return A 2d array with the median (or max, or mean, etc) of the input array
#' @export
SummaryPatristicMatrixArray <- function(patristic.matrix.array,fn=stats::median) {
  return(apply(patristic.matrix.array,MARGIN=c(1,2),fn, na.rm=TRUE))
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


#' Get all calibrations given trees in database
#' @param input A vector of names, a newick string, or a phylo object
#' @param partial If TRUE, use source trees even if they only match some of the desired taxa
#' @param usetnrs If TRUE, use OpenTree's services to resolve names. This can dramatically improve the chance of matches, but also take much longer
#' @param approximatematch If TRUE, use a slower TNRS to correct mispellings, increasing the chance of matches (including false matches)
#' @param cache The cached set of chronograms and other info from data(opentree_chronograms)
#' @return data.frame of calibrations
#' @export
GetAllCalibrations <- function(input=c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"), partial=TRUE, usetnrs=FALSE, approximatematch=TRUE, cache=get("opentree_chronograms")) {
	phylo.results <- EstimateDates(input=input, partial=partial, usetnrs=usetnrs, approximatematch=approximatematch, cache=cache, output.format="phylo.all")
	constraints.df <- data.frame()
	for (i in sequence(length(phylo.results))) {
		local.df <- geiger::congruify.phylo(reference=phylo.results[[i]], target=phylo.results[[i]], scale=NA)$calibrations
		local.df$reference <- names(phylo.results)[i]
		if(i==1) {
			constraints.df <- local.df
		} else {
			constraints.df <- rbind(constraints.df, local.df)
		}
	}
	return(constraints.df)
}

#' Use all calibrations given trees in database to date a tree
#' @param phy A phylo object
#' @param partial If TRUE, use source trees even if they only match some of the desired taxa
#' @param usetnrs If TRUE, use OpenTree's services to resolve names. This can dramatically improve the chance of matches, but also take much longer
#' @param approximatematch If TRUE, use a slower TNRS to correct mispellings, increasing the chance of matches (including false matches)
#' @param cache The cached set of chronograms and other info from data(opentree_chronograms)
#' @param expand How much to expand by each step to get consistent calibrations
#' @param giveup How many expansion to try before giving up
#' @return list with chronogram, original calibrations, and expanded calibrations
#' @export
#' @details
#' This will try to use the calibrations as fixed ages.
#' If that fails (often due to conflict between calibrations), it will expand the range of the minage and maxage and try again. And repeat.
#' expand sets the expansion value: should be between 0 and 1
UseAllCalibrations <- function(phy=GetBoldOToLTree(c("Rhea americana",  "Struthio camelus", "Gallus gallus")), partial=TRUE, usetnrs=FALSE, approximatematch=TRUE, cache=get("opentree_chronograms"), expand=0.1, giveup=100) {
	calibrations.df <- GetAllCalibrations(input=gsub('_', ' ', phy$tip.label), partial=partial, usetnrs=usetnrs, approximatematch=approximatematch, cache=cache)
	phy$tip.label <- gsub(' ', '_', phy$tip.label) #underscores vs spaces: the battle will never end.
	calibrations.df$taxonA <- gsub(' ', '_', calibrations.df$taxonA)
	calibrations.df$taxonB <- gsub(' ', '_', calibrations.df$taxonB)
	calibrations.df <- calibrations.df[which(calibrations.df$taxonA %in% phy$tip.label),]
	calibrations.df <- calibrations.df[which(calibrations.df$taxonB %in% phy$tip.label),]
	original.calibrations.df <- calibrations.df
	chronogram <- NULL
	try(chronogram <- geiger::PATHd8.phylo(phy, calibrations.df), silent=TRUE)
	attempts=0
	if(expand!=0) {
		while(is.null(chronogram) & attempts<giveup) {
			calibrations.df <- original.calibrations.df
			calibrations.df$MaxAge <- calibrations.df$MaxAge * ((1+expand)^attempts)
			calibrations.df$MinAge <- calibrations.df$MinAge * ((1-expand)^attempts)

			# We will have no fixed ages. Pathd8 just quietly gives up. So instead, we add a tiny branch with a zero calibration
			# between it and its sister.
			made.up.edgelength <- min(1e-9, .001*min(phy$edge.length))
			phy2 <- phytools::bind.tip(ape::reorder.phylo(phy), "tinytip", edge.length=made.up.edgelength, where=1, position=made.up.edgelength) #bind tip has weird behavior for non-reordered trees
			calibrations.df[dim(calibrations.df)[1]+1,]<- c("fixed", 0, 0, phy$tip.label[1], "tinytip", "none")
			try(chronogram <- geiger::PATHd8.phylo(phy2, calibrations.df))
			if(!is.null(chronogram)) {
				chronogram <- ape::drop.tip(chronogram, "tinytip")
			}
			attempts <- attempts+1
		}
		if(attempts>0) {
			print("Dates are even more approximate than usual: had to expand constraints to have them agree")
		}
	}
	return(list(phy=chronogram, calibrations.df=calibrations.df, original.calibrations.df=original.calibrations.df))
}

#' Use Barcode of Life data to get branch lengths on the OToL tree
#' @param input A vector of names
#' @param marker Gene name to select
#' @param otol_version Version of OToL to use
#' @param doML Boolean; if TRUE, does ML brlen optimization
#' @return A phylogeny with ML branch lengths
#' @export
GetBoldOToLTree <- function(input=c("Rhea americana",  "Struthio camelus","Gallus gallus"), marker="COI", otol_version="v2", doML=FALSE) {
	#otol returns error with missing taxa in v3 of rotl
	phy <- ape::multi2di(rotl::tol_induced_subtree(ott_ids=rotl::tnrs_match_names(names=input)$ott_id, label_format="name",  otl_v=otol_version))
	phy$tip.label <- gsub("_ott.*","", phy$tip.label)
	sequences <- bold::bold_seqspec(taxon=input, marker=marker)
	final.sequences <- matrix("-", nrow=length(input), ncol=max(sapply(strsplit(sequences$nucleotides, ""), length)))
	final.sequences.names <- rep(NA, length(input))
	for (i in sequence(dim(sequences)[1])) {
		taxon <- sequences$species_name[i]
		if(!(taxon %in% final.sequences.names)) {
			seq <- strsplit(sequences$nucleotide[i],"")[[1]]
			matching.index <- 1+sum(!is.na(final.sequences.names))
			final.sequences[matching.index, sequence(length(seq))] <- seq
			final.sequences.names[matching.index] <- taxon
		}
	}
	rownames(final.sequences) <- gsub(" ", "_", final.sequences.names)
	final.sequences <- final.sequences[!is.na(final.sequences.names),]

	alignment <- ape::as.DNAbin(final.sequences)
	alignment <- phangorn::as.phyDat(ips::mafft(alignment))
	taxa.to.drop <- phy$tip.label[which(!phy$tip.label %in% rownames(final.sequences))]
	if(length(taxa.to.drop)>0) {
		phy <- ape::drop.tip(phy, taxa.to.drop)
	}
	pml.object <- phangorn::pml(phangorn::acctran(phy, alignment), data=alignment)
	pml.object$tree <- ape::chronoMPL(pml.object$tree, se=FALSE, test=FALSE)
	phy <- pml.object$tree
	if(doML) {
		phy <- phangorn::optim.pml(pml.object, data=alignment, rearrangement="none", optRooted=TRUE, optQ=TRUE)$tree
	}
	phy$tip.label <- gsub('_', ' ', phy$tip.label)
	return(phy)
}




#
# ReadDistance <- function(file) {
# 	data <- readLines(file, n=-1)[-1] #read in phylip distance, perhaps from SDM
# 	data <- strsplit(data, " +")
# 	data[sapply(data, length)>0] #trim trailing
# 	final.matrix <- matrix(nrow=length(data), ncol=length(data))
# 	all.names <- rep(NA, length(data))
# 	for (data.index in sequence(length(data))) {
# 		local.line <- data[[data.index]]
# 		all.names[data.index] <- local.line[1]
# 		local.line <- as.numeric(local.line[-1])
#
# 	}
# }

#
# SDM <- function(filtered.results, weights=NULL)) {
#
# 	patristic.array <- BindMatrices(filtered.results)
# 	k <- length(filtered.results)
# 	if(is.null(weights)) {
# 		weights <- rep(1/k, k)
# 	} else {
# 		weights <- weights/sum(weights) #just to make sure total weight is 1
# 	}
#
# 	#Use their appendix and stick it in solve
# 	#a*x=b
# 	#The a matrix has a vector of alpha values, then a_i for matrix 1, a_i for matrix 2..., then
#
# }
#
# #rewrite of code in ape by Andrei Popescu niteloserpopescu@gmail.com to allow for debugging
# apeSDM <- function(...) {
# 	st <- list(...)
# 	k <- length(st)/2
# 	ONEtoK <- seq_len(k)
# 	for (i in ONEtoK) st[[i]] <- as.matrix(st[[i]])
# 	ROWNAMES <- lapply(st[ONEtoK], rownames)
# 	NROWS <- sapply(ROWNAMES, length)
# 	tot <- sum(NROWS)
# 	labels <- unique(unlist(ROWNAMES))
# 	sp <- unlist(st[k + ONEtoK])
# 	astart <- numeric(tot)
# 	astart[1] <- k
# 	for (i in 2:k) astart[i] <- astart[i - 1] + NROWS[i - 1]
# 	n <- length(labels)
# 	miustart <- k + tot
# 	niustart <- miustart + n
# 	lambstart <- niustart + k - 1
# 	X <- matrix(0, n, n, dimnames = list(labels, labels))
# 	V <- w <- X
# 	tmp <- 2 * k + tot + n
# 	col <- numeric(tmp)
# 	for (i in 1:(n - 1)) {
# 		for (j in (i + 1):n) {
# 			for (p in ONEtoK) {
# 				if (is.element(labels[i], ROWNAMES[[p]]) && is.element(labels[j],
# 					ROWNAMES[[p]])) {
# 						w[i, j] <- w[j, i] <- w[i, j] + sp[p]
# 					}
# 				}
# 			}
# 		}
# 		ONEtoN <- seq_len(n)
# 		Q <- matrix(0, tmp, tmp)
# 		for (p in ONEtoK) {
# 			d_p <- st[[p]]
# 			for (l in ONEtoK) {
# 				d <- st[[l]]
# 				sum <- 0
# 				dijp <- -1
# 				if (l == p) {
# 					for (i in ONEtoN) {
# 						for (j in ONEtoN) {
# 							if (i == j)
# 							next
# 							pos <- match(labels[c(i, j)], ROWNAMES[[l]])
# 							if (all(!is.na(pos))) {
# 								ipos <- pos[1L]
# 								jpos <- pos[2L]
# 								dij <- d[ipos, jpos]
# 								sum <- sum + dij * dij - sp[l] * dij *
# 								dij/w[i, j]
# 								tmp2 <- dij - sp[l] * dij/w[i, j]
# 								Q[p, astart[l] + ipos] <- Q[p, astart[l] +
# 								ipos] + tmp2
# 								Q[p, astart[l] + jpos] <- Q[p, astart[l] +
# 								jpos] + tmp2
# 							}
# 						}
# 					}
# 				}
# 				else {
# 					for (i in ONEtoN) {
# 						for (j in ONEtoN) {
# 							if (i == j)
# 							next
# 							pos <- match(labels[c(i, j)], ROWNAMES[[l]])
# 							posp <- match(labels[c(i, j)], ROWNAMES[[p]])
# 							if (all(!is.na(pos)) && all(!is.na(posp))) {
# 								ipos <- pos[1L]
# 								jpos <- pos[2L]
# 								dij <- d[ipos, jpos]
# 								dijp <- d_p[posp[1L], posp[2L]]
# 								sum <- sum - sp[l] * dij * dijp/w[i, j]
# 								tmp2 <- sp[l] * dijp/w[i, j]
# 								Q[p, astart[l] + ipos] <- Q[p, astart[l] +
# 								ipos] - tmp2
# 								Q[p, astart[l] + jpos] <- Q[p, astart[l] +
# 								jpos] - tmp2
# 							}
# 						}
# 					}
# 				}
# 				Q[p, l] <- sum
# 			}
# 			Q[p, lambstart + 1] <- 1
# 		}
# 		r <- k
# 		for (p in ONEtoK) {
# 			dp <- st[[p]]
# 			for (i in ONEtoN) {
# 				if (is.element(labels[i], ROWNAMES[[p]])) {
# 					r <- r + 1
# 					for (l in ONEtoK) {
# 						d <- st[[l]]
# 						if (l == p) {
# 							ipos <- match(labels[i], ROWNAMES[[p]])
# 							for (j in ONEtoN) {
# 								if (i == j)
# 								next
# 								jpos <- match(labels[j], ROWNAMES[[p]])
# 								if (!is.na(jpos)) {
# 									dij <- d[ipos, jpos]
# 									Q[r, l] <- Q[r, l] + dij - sp[l] * dij/w[i,
# 									j]
# 									tmp2 <- 1 - sp[l]/w[i, j]
# 									Q[r, astart[l] + ipos] <- Q[r, astart[l] +
# 									ipos] + tmp2
# 									Q[r, astart[l] + jpos] <- Q[r, astart[l] +
# 									jpos] + tmp2
# 								}
# 							}
# 						}
# 						else {
# 							for (j in ONEtoN) {
# 								if (i == j)
# 								next
# 								if (!is.element(labels[j], rownames(dp)))
# 								next
# 								pos <- match(labels[c(i, j)], ROWNAMES[[l]])
# 								if (all(!is.na(pos))) {
# 									ipos <- pos[1L]
# 									jpos <- pos[2L]
# 									dij <- d[ipos, jpos]
# 									Q[r, l] <- Q[r, l] - sp[l] * dij/w[i,
# 									j]
# 									tmp2 <- sp[l]/w[i, j]
# 									Q[r, astart[l] + ipos] <- Q[r, astart[l] +
# 									ipos] - tmp2
# 									Q[r, astart[l] + jpos] <- Q[r, astart[l] +
# 									jpos] - tmp2
# 								}
# 							}
# 						}
# 					}
# 					if (p < k)
# 					Q[r, ] <- Q[r, ] * sp[p]
# 					Q[r, miustart + i] <- 1
# 					if (p < k)
# 					Q[r, niustart + p] <- 1
# 				}
# 			}
# 		}
# 		r <- r + 1
# 		col[r] <- k
# 		Q[r, ONEtoK] <- 1
# 		for (i in ONEtoN) {
# 			r <- r + 1
# 			for (p in ONEtoK) {
# 				ipos <- match(labels[i], ROWNAMES[[p]])
# 				if (!is.na(ipos))
# 				Q[r, astart[p] + ipos] <- 1
# 			}
# 		}
# 		for (p in 1:(k - 1)) {
# 			r <- r + 1
# 			for (i in ONEtoN) {
# 				ipos <- match(labels[i], ROWNAMES[[p]])
# 				if (!is.na(ipos))
# 				Q[r, astart[p] + ipos] <- 1
# 			}
# 		}
# 		a <- solve(Q, col, 1e-19)
# 		for (i in ONEtoN) {
# 			for (j in ONEtoN) {
# 				if (i == j) {
# 					X[i, j] <- V[i, j] <- 0
# 					next
# 				}
# 				sum <- 0
# 				sumv <- 0
# 				for (p in ONEtoK) {
# 					d <- st[[p]]
# 					pos <- match(labels[c(i, j)], ROWNAMES[[p]])
# 					if (all(!is.na(pos))) {
# 						ipos <- pos[1L]
# 						jpos <- pos[2L]
# 						dij <- d[ipos, jpos]
# 						sum <- sum + sp[p] * (a[p] * dij + a[astart[p] +
# 							ipos] + a[astart[p] + jpos])
# 							sumv <- sumv + sp[p] * (a[p] * dij)^2
# 						}
# 					}
# 					X[i, j] <- sum/w[i, j]
# 					V[i, j] <- sumv/(w[i, j])^2
# 				}
# 	}
# 	list(X, V)
# }
#
#
#
#
#
#' Function to compute the SDM supertree Criscuolo et al. 2006
#' @param filtered.results List of patristic matrices
#' @param weighting flat, taxa, inverse
#' @return A list containing phy (a chronogram), and filtered.results that were actually used
#' @export
#' @details
#' Weighting is how much weight to give each input tree.
#'    flat = all trees have equal weighting
#'    taxa = weight is proportional to number of taxa
#'    inverse = weight is proportional to 1 / number of taxa
#' Criscuolo A, Berry V, Douzery EJ, Gascuel O. SDM: a fast distance-based approach for (super) tree building in phylogenomics. Syst Biol. 2006;55(5):740–55. doi: 10.1080/10635150600969872.
RunSDM <- function(filtered.results, weighting="flat") {
	phy <- NA
	used.studies <- names(filtered.results)
	unpadded.matrices <- lapply(filtered.results, UnpadMatrix)
	good.matrix.indices <- c()
	for(i in sequence(length(unpadded.matrices))) {
		test.result <- NA
		# Rationale here: some chronograms always cause errors with SDM, even when trying to get a consensus of them
		# with themselves. For now, throw out of synthesis.
		try(test.result <- mean(do.call(ape::SDM, c(unpadded.matrices[i], unpadded.matrices[i], rep(1, 2)))[[1]]), silent=TRUE)
		if(is.finite(test.result)) {
			good.matrix.indices <- append(good.matrix.indices,i)
		}
	}
	if(length(good.matrix.indices)>0) {
		unpadded.matrices <- unpadded.matrices[good.matrix.indices]
		used.studies <- used.studies[good.matrix.indices]
		weights = rep(1, length(unpadded.matrices))
		if (weighting=="taxa") {
			weights = unname(sapply(unpadded.matrices, dim)[1,])
		}
		if (weighting=="inverse") {
			weights = 1/unname(sapply(unpadded.matrices, dim)[1,])
		}
		SDM.result <- do.call(ape::SDM, c(unpadded.matrices, weights))[[1]]
		#agnes in package cluster has UPGMA with missing data; might make sense here
		try(phy <- phangorn::upgma(SDM.result))
	}
	return(list(phy=phy, filtered.results=unpadded.matrices))
}
