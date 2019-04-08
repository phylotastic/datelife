
#' Convert patristic matrix to a phylo object. Used inside: summarize_datelife_result, CongruiyTree.
#' @param patristic_matrix A patristic matrix
#' @param clustering_method A character vector indicating the method to construct the tree. Options are "nj" for Neighbor-Joining and "upgma" for Unweighted Pair Group Method with Arithmetic Mean.
# We might add the option to insert a function as clustering_method.
# Before, we hard coded it to try Neighbor-Joining first; if it errors, it will try UPGMA.\
# Now, it uses nj for phylo_all summary,
#' @param fix_negative_brlen Boolean indicating whether to fix negative branch lengths in resulting tree or not. Default to TRUE.
#' @inheritParams tree_fix_brlen
#' @return A rooted phylo object
#' @export
patristic_matrix_to_phylo <- function(patristic_matrix, clustering_method = "nj", fix_negative_brlen = TRUE, fixing_method = 0, ultrametric = TRUE) {
    # patristic_matrix <- threebirds_result[[5]]
    if(!inherits(patristic_matrix, "matrix") & !inherits(patristic_matrix, "data.frame")){
        message("patristic_matrix argument is not a matrix")
        return(NA)
    }
    # has to be matrix not data frame:
    if(inherits(patristic_matrix, "data.frame")){
        patristic_matrix <- as.matrix(patristic_matrix)
    }
    clustering_method <- match.arg(arg = clustering_method, choices = c("nj", "upgma"), several.ok = FALSE)
    if(anyNA(patristic_matrix)) {
  	   patristic_matrix <- patristic_matrix[rowSums(is.na(patristic_matrix)) != ncol(patristic_matrix),colSums(is.na(patristic_matrix)) != nrow(patristic_matrix)]
   }   # Get rid of rows and columns with all NA or NaNs, leaves the ones with some NA or NaNs
    if(dim(patristic_matrix)[1] == 2) {
        phy <- ape::rtree(n = 2, rooted = TRUE, tip.label = rownames(patristic_matrix), br = 0.5 * patristic_matrix[1,2])
        phy$clustering_method <- "ape::rtree"
        phy$citation <- names(patristic_matrix)
        return(phy)
    }
    phycluster <- cluster_patristicmatrix(patristic_matrix)
    phy <- choose_cluster(phycluster, clustering_method)
    if(!inherits(phy, "phylo")){
      message("Clustering patristic matrix to phylo failed.")
      phy$citation <- names(patristic_matrix)
      return(phy)
    }
    phy$negative_brlen <- NA
    mess1 <- "Converting from patristic distance matrix to a tree resulted in some negative branch lengths;"
    mess2 <- "the largest by magnitude is"
    # when original tree IS ultrametric and has negative edges:
    if(is.null(phy$edge.length.original) & any(phy$edge.length < 0)){
      warning(paste(mess1, mess2, min(phy$edge.length)))
    }
    # when original tree is NOT ultrametric and has negative edges:
    if(!is.null(phy$edge.length.original) & any(phy$edge.length.original < 0)){
      warning(paste(mess1, mess2, min(phy$edge.length.original)))
      # and when tree was forced ultrametric and there still are neg edges:
      if(any(phy$edge.length < 0)){
        warning(paste("After phytools::forcing.ultrametric there still are negative branch lengths;", mess2, min(phy$edge.length)))
      }
    }
    if(any(phy$edge.length < 0)){
        phy$negative_brlen <- list(edge_number = which(phy$edge.length < 0))
        # phy$edge.length[which(phy$edge.length < 0)] <- 0 #sometimes NJ returns tiny negative branch lengths. https://github.com/phylotastic/datelife/issues/11
        if(fix_negative_brlen){
            phy$negative_brlen <- list(edge_number = which(phy$edge.length < 0))
            phy <- tree_fix_brlen(tree = phy, fixing_criterion = "negative", fixing_method = fixing_method)
            fixing_method_called <- as.list(environment())$fixing_method
            phy$negative_brlen <- c(phy$negative_brlen, list(fixing_method = fixing_method_called))
            warning(paste0("Negative branch lengths were fixed with tree_fix_brlen, fixing_method = ", fixing_method_called))
        }
    }
    # for cases when there are no neg branch lengths to fix (or we don't want them fixed)
    # and we still want the final tree to be ultrametric:
    if(ultrametric){
      if(is.null(phy$edge.length.original)){
          phy <- force_ultrametric(phy)
      }
    }
    phy$citation <- names(patristic_matrix)
    class(phy) <- c(class(phy), "datelifeTree")
    return(phy)
}
#' Force a non ultrametric phylo object to be ultrametric
#' @inheritParams phylo_check
#' @return A phylo object
#' @export
force_ultrametric <- function(phy){
      # enhance: check if there is an edge.length.original already There
      # something like how many grepl("edge.length.original") in names(phy) and add an integer after it.
      if(!inherits(phy, "phylo")){
          message("phy argument is not a phylo object.")
          return(NA)
      }
      if(!ape::is.ultrametric(phy)){
        phy$edge.length.original <- phy$edge.length
        phy <- phytools::force.ultrametric(tree = phy, method = "extend")
        phy$force.ultrametric <- "extend"
      }
      return(phy)
}
#' Cluster a patristic matrix with several methods
#'
#' @inheritParams patristic_matrix_to_phylo
#' @return A list of results from clustering with NJ and UPGMA
#' @export
cluster_patristicmatrix <- function(patristic_matrix){
    if(!inherits(patristic_matrix, "matrix") & !inherits(patristic_matrix, "data.frame")){
        message("patristic_matrix argument is not a matrix")
        return(NA)
    }
    # has to be matrix not data frame:
    if(inherits(patristic_matrix, "data.frame")){
        patristic_matrix <- as.matrix(patristic_matrix)
    }
  if(dim(patristic_matrix)[1] < 2) {
     return(NA)
  }
  if(dim(patristic_matrix)[1] == 2) {
      message("patristic_matrix is dimension 2, you don't need to cluster.")
      return(NA)
  } else {
        phyclust <- vector(mode = "list", length = 4)
        names(phyclust) <- c("nj", "njs", "upgma", "upgma_daisy")
        phyclust$nj <- tryCatch(ape::nj(patristic_matrix), error = function (e) NULL)
        if (is.null(phyclust$nj)){ # case when we have missing data (NA) on patristic_matrix and regular nj does not work; e.g. clade thraupidae SDM.results have missing data, and upgma chokes
            # njs appears to be the only option for missing data with NJ
            # but see Criscuolo and Gascuel. 2008. Fast NJ-like algorithms to deal with incomplete distance matrices. BMC Bioinformatics 9:166
            phyclust$njs <- tryCatch(ape::njs(patristic_matrix), error = function(e) NULL)
            if(!is.null(phyclust$njs)){
              phyclust$njs <- tryCatch(phytools::midpoint.root(phyclust$njs),
                          error = function(e) NULL)
            }
        } else {
            phyclust$nj <- tryCatch(phytools::midpoint.root(phyclust$njs),
                      error = function(e) NULL)
        }
        # root the tree on the midpoint (only for trees with ape::Ntip(phy) > 2):
        # phy <- tryCatch(phangorn::midpoint(phy), error = function(e) NULL)
        # using phytools::midpoint.root instead of phangorn::midpoint bc it's less error prone.
        # sometimes, nj and njs do not work if patristic matrices come from sdm. why? it was probably the midpoint function from phangorn. Using phytools one now.
        # use regular upgma (or implemented with daisy and hclust) when nj or midpoint.root fail
        phyclust$upgma <- tryCatch(phangorn::upgma(patristic_matrix), error = function (e) NULL)
        if (is.null(phyclust$upgma)){ # case when we have missing data (NA) on patristic_matrix and regular upgma does not work; e.g. clade thraupidae SDM.results have missing data, and upgma chokes
            phyclust$upgma_daisy <- tryCatch({
              # using daisy to calculate dissimilarity matrix instead of as.dist (which is used in phangorn::upgma) when there are NAs in the matrix; agnes does not work with NAs either.
              patristic_matrix <- patristic_matrix*0.25 # doing this because it's giving ages that are too old, but it's a random number
              DD <- cluster::daisy(x = patristic_matrix, metric = "euclidean")
              hc <- stats::hclust(DD, method = "average") # original clustering method from phangorn::upgma. Using agnes() instead hclust() to cluster gives the same result.
              phy <- ape::as.phylo(hc)
              phy <- phylobase::reorder(phy, "postorder")
              phy
            }, error = function(e) NULL)
        }
        return(phyclust)
  }
}
#' Choose a phylo object from cluster_patristicmatrix. If not ultrametric, it does not force it.
#'
#' @inheritParams patristic_matrix_to_phylo
#' @param phycluster An output from cluster_patristicmatrix
#' @return A phylo object or NA
#' @export
choose_cluster <- function(phycluster, clustering_method){
    if(!inherits(phycluster, "list")){
        message("phycluster arguments is not a list; check it out.")
        return(NA)
    }
  phy_return <- NA
  if(length(phycluster) == 0){
    message("phycluster argument is length 0")
    return(NA)
  }
  if(inherits(phycluster, "phylo")){ # it is a tree of two tips
    return(phycluster)
  } else { # it is a list of results from cluster_patristicmatrix
    fail <- sapply(phycluster, is.null)
    if(all(fail)){
        message("The patristic matrix could not be transformed into a tree with any of the default methods (NJ, UPGMA)")
        return(NA)
    }
    phycluster <- phycluster[!fail] # take out the fails or unattempted from cluster_patristicmatrix
    if(length(phycluster) == 1){
      phy <- phycluster[[1]]
      phy$clustering_method <- names(phycluster)
      # if(!ape::is.ultrametric(phy)){
      #   phy$edge.length.original <- phy$edge.length
      #   phy <- phytools::force.ultrametric(tree = phy, method = "extend")
      #   phy$force.ultrametric <- "extend"
      # }
      return(phy)
    } else {
      ultram <- sapply(phycluster, ape::is.ultrametric)
      ultram2 <- sapply(phycluster, ape::is.ultrametric, 2)
      if(length(ultram) == 0 & length(ultram2) == 0){
        message("The patristic matrix could not be transformed into an ultrametric tree with any of the default methods (NJ, UPGMA)")
        # return(NA)
      }
      choice <- grepl(clustering_method, names(phycluster)) # choice can only be one
      ff <- which(choice & ultram) # if the chosen method gives an ultrametric tree
      if(length(ff) == 1){
        phy <- phycluster[[ff]]
        phy$clustering_method <- names(phycluster)[ff]
        return(phy)
      }
      ff <- which(!choice & ultram) # if not, take the not chosen but ultrametric
      if(length(ff) == 1){
        phy <- phycluster[[ff]]
        phy$clustering_method <- names(phycluster)[ff]
        return(phy)
      }
      ff <- which(choice & ultram2) # if not, take the chosen one but less ultrametric
      if(length(ff) == 1){
        phy <- phycluster[[ff]]
        # phy$edge.length.original <- phy$edge.length
        # phy <- phytools::force.ultrametric(tree = phy, method = "extend")
        # phy$force.ultrametric <- "extend"
        phy$clustering_method <- names(phycluster)[ff]
        return(phy)
      }
      ff <- which(!choice & ultram2) # if not, take the not chosen one but less ultrametric
      if(length(ff) == 1){
        phy <- phycluster[[ff]]
        # phy$edge.length.original <- phy$edge.length
        # phy <- phytools::force.ultrametric(tree = phy, method = "extend")
        # phy$force.ultrametric <- "extend"
        phy$clustering_method <- names(phycluster)[ff]
        return(phy)
      }
    }
  }
}
#' Go from a smmary matrix to an ultrametric phylo object.
#' @param summ_matrix A summary patristic distance matrix from sdm or median. See details.
#' @param total_distance Boolean. If TRUE it will divide the matrix byhalf, if FALSE it will take iy as is.
#' @param use A character vector indicating what type of age to use for summary. One of the following
#' \describe{
#'	\item{mean}{ It will use the mean of the node age distributions.
#'	}
#'	\item{min}{ It will use the minimum age from the node age distrbutions.
#'	}
#'	\item{max}{ Choose this if you wanna be conservative; it will use the maximum age from the node age distrbutions.
#'	}
#' }
#' @param target_tree A phylo object. Use this in case you want a particular backbone for the output tree.
#' @inheritDotParams get_otol_synthetic_tree
#' @return An ultrametric phylo object.
#' @details It can take a regular patristic distance matrix, but there are simpler methods for that implemented in patristic_matrix_to_phylo.
#' @export
summary_matrix_to_phylo <- function(summ_matrix, total_distance = TRUE, use = "mean", target_tree = NULL, ...){ # enhance: allow other methods, not only bladj.
    # for debugging here
    # summ_matrix <- subset2_sdm_matrix
    # summ_matrix <- median.matrix
    use_age <- match.arg(use, c("mean", "median", "min", "max"))
    if(!inherits(summ_matrix, "matrix") & !inherits(summ_matrix, "data.frame")){
        message("summ_matrix argument is not a matrix")
        return(NA)
    }
    # summ_matrix <- data.frame(summ_matrix)
    # everything up to patristic_matryx_to_phylo ok if it is a data frame too
    if(inherits(summ_matrix, "data.frame")){
        summ_matrix <- as.matrix(summ_matrix)
        colnames(summ_matrix) <- gsub("\\.", " ", colnames(summ_matrix))
    }
    if(total_distance){
      summ_matrix <- summ_matrix * 0.5  # bc it's total distance tip to tip
    }
	ages <- tA <- tB <- c()
    # to compute the final length of the data frame do: ncol(xx)^2 - sum(1:(ncol(xx)-1))
	# calibrations <- matrix(nrow = ncol(xx)^2 - sum(1:(ncol(xx)-1)), ncol = 3)
	# identify if SDM matrix has some negative values; extract taxon names:
	negs <- which(summ_matrix < 0)
	neg_names <- rownames(summ_matrix)[ceiling(negs/nrow(summ_matrix))]
	# extract unique ages from summ_matrix:
	for(i in seq(ncol(summ_matrix))){
		# calibrations[start:start+i,1] <- xx[1:i,i]
		# calibrations[start:start+i,2] <- rownames(xx)[1:i]
		# calibrations[start:start+i,3] <- rep(colnames(xx)[i], i)
		# start <- sum(!is.na(calibrations[,1])) +1
		ages <- c(ages, summ_matrix[1:i,i])
		tA <- c(tA, rownames(summ_matrix)[1:i])
		tB <- c(tB, rep(colnames(summ_matrix)[i], i))
	}
    tA <- gsub(" ", "_", tA)
    tB <- gsub(" ", "_", tB)
	calibrations <- data.frame(Age = ages, taxonA = tA, taxonB = tB, stringsAsFactors = FALSE)
	# calibrations$taxonA <- as.character(calibrations$taxonA)
	# calibrations$taxonB <- as.character(calibrations$taxonB)
	calibrations <- calibrations[!is.na(calibrations[,"Age"]), ] # get rid of NaN and NAs
	calibrations <- calibrations[calibrations[,"Age"] != 0, ] # get rid of 0's
	calibrations <- calibrations[calibrations[,"Age"] > 0, ] # get rid of negative values too
	if(any(is.na(calibrations[,"Age"]))){
		warning("for some reason there are still NAs in the matrix")}
	# enhance: where does this negative values come from in SDM?
	# get a backbone tree:
	# chronogram <- geiger::PATHd8.phylo(phy_target, calibrations)
	# try(chronogram <- geiger::PATHd8.phylo(phy_target, calibrations), silent = TRUE)
  if(!inherits(target_tree, "phylo")){
    target_tree <- suppressMessages(get_otol_synthetic_tree(...))
    if(!inherits(target_tree, "phylo")){
      # we should find a better way to do this, but it should be ok for now:
      target_tree <- suppressWarnings(suppressMessages(patristic_matrix_to_phylo(summ_matrix, ultrametric = TRUE)))
      # target_tree <- consensus(phyloall, p = 0.5) # can't use consensus here: not all trees have the same number of tips, duh
    }
    target_tree <- ape::collapse.singles(target_tree)
      # ape::is.ultrametric(target_tree)
      # ape::is.binary(target_tree)
      # plot(target_tree, cex = 0.5)
  }
  if(!inherits(target_tree, "phylo")){
		message("target_tree is not phylo object; returning NA")
        message("Hint: Was summ_matrix constructed from an object with no good groves? try running get_best_grove first")
		# enhance: add a more formal test of best grove
		return(NA)
  }
  target_tree$edge.length <- NULL
  target_tree$edge.length.original <- NULL
  target_tree$tip.label <- gsub(" ", "_", target_tree$tip.label)
  # test that taxonA and taxonB are all in target tree tip labels
  if(any(is.na(match(unique(c(calibrations$taxonA, calibrations$taxonB)), target_tree$tip.label)))){
    message("taxa in summ_matrix are not in target_tree; returning NA")
    return(NA)
  }
	# get the coincident node numbers:
  # ape::is.binary(target_tree)
	target_tree_nodes <- sapply(seq(nrow(calibrations)), function(i)
		phytools::findMRCA(tree = target_tree,
			tips = as.character(calibrations[i,c("taxonA", "taxonB")]),
			type = "node"))
	target_tree_nodes <- target_tree_nodes - ape::Ntip(target_tree)
	all_nodes <- sort(unique(target_tree_nodes))
	# get the node age distribution:
	all_ages <- lapply(all_nodes, function(i) calibrations[target_tree_nodes == i, "Age"])
	# any(sapply(all_ages, is.null)) # if FALSE, all nodes have at least one calibration.
	calibrations2 <- data.frame(MRCA = paste0("n", all_nodes), MinAge = sapply(all_ages, min), MaxAge= sapply(all_ages, max))
	# calibrations2$MRCA is a factor so have to be made as.character to work with bladj
	if(all(all_nodes < ape::Ntip(target_tree))){
        all_nodes_numbers <- all_nodes + ape::Ntip(target_tree)
		node_index <- "consecutive"
	} else {
    all_nodes_numbers <- all_nodes
		node_index <- "node_number"
	}
	target_tree$node.label <- NULL # make sure its nullso we can rename all nodes of interest to match our labels
	target_tree <- tree_add_nodelabels(tree = target_tree, node_index = node_index)  # all nodes need to be named so make_bladj_tree runs properly
  if("mean" %in% use){
    node_ages <- sapply(seq(nrow(calibrations2)), function(i) sum(calibrations2[i,c("MinAge", "MaxAge")])/2)
  }
  if("min" %in% use){
    node_ages <- calibrations2[,c("MinAge")]
  }
  if("max" %in% use){
    node_ages <- calibrations2[,c("MaxAge")]
  }
  new_phy <- make_bladj_tree(tree = target_tree, nodenames = as.character(calibrations2$MRCA), nodeages = node_ages)
  new_phy$clustering_method <- "datelife"
  new_phy$calibrations_distribution <- stats::setNames(all_ages, all_nodes)
  new_phy$calibrations_MIN <- calibrations2$MinAge
  new_phy$calibrations_MAX <- calibrations2$MaxAge
  new_phy$calibrations_MRCA <- all_nodes_numbers
  new_phy$ott_ids <- NULL
  if(!is.null(target_tree$ott_ids)){
      tt <- match(new_phy$tip.label, target_tree$tip.label)
      # match(c("a", "b", "c", "d"), c("c", "d", "a", "a", "a", "b"))
      new_phy$ott_ids <- target_tree$ott_ids[tt]
  }
  return(new_phy)
}
