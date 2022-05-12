#' Convert a patristic matrix to a `phylo` object.
#'
#' Function `patristic_matrix_to_phylo` is used inside [summarize_datelife_result()].
#'
#' @param patristic_matrix A patristic matrix
#' @param clustering_method A character vector indicating the method to construct
#' the tree. Options are:
#' \describe{
#' 	\item{nj}{Neighbor-Joining method applied with [ape::nj()].}
#' 	\item{upgma}{Unweighted Pair Group Method with Arithmetic Mean method applied
#'        with [phangorn::upgma()].}
#' 	\item{bionj}{An improved version of the Neighbor-Joining method applied with
#'       [ape::bionj()].}
#' 	\item{triangle}{Triangles method applied with [ape::triangMtd()]}
#' 	\item{mvr}{Minimum Variance Reduction method applied with [ape::mvr()].}
#' }
#' @details
#' We might add the option to insert a function as `clustering_method` in the future.
#' Before, we had hard-coded the function to try Neighbor-Joining (NJ) first; if it
#'  errors, it will try UPGMA.
#' Now, it uses NJ for a "phylo_all" summary, and we are using our own algorithm to
#'  get a tree from a summary matrix.
#' @param fix_negative_brlen Boolean indicating whether to fix negative branch
#'  lengths in resulting tree or not. Default to `TRUE`.
#' @param variance_matrix A variance matrix from a `datelifeResult` object,
#'  usually an output from [datelife_result_variance_matrix()].
#'  Only used if `clustering_method = "mvr"`.
#' @inheritParams tree_fix_brlen
#' @return A rooted `phylo` object.
#' @export
patristic_matrix_to_phylo <- function(patristic_matrix,
                                      clustering_method = "nj",
                                      fix_negative_brlen = TRUE,
                                      fixing_method = 0,
                                      ultrametric = TRUE,
                                      variance_matrix = NULL) {
  # patristic_matrix <- threebirds_result[[5]]
  if (!inherits(patristic_matrix, "matrix") & !inherits(patristic_matrix, "data.frame")) {
    message("patristic_matrix argument is not a matrix")
    return(NA)
  }
  # has to be matrix not data frame:
  if (inherits(patristic_matrix, "data.frame")) {
    patristic_matrix <- as.matrix(patristic_matrix)
  }
  clustering_method <- match.arg(arg = tolower(clustering_method),
                                 choices = c("nj", "upgma", "bionj", "triangle", "mvr"),
                                 several.ok = FALSE)
  if (anyNA(patristic_matrix)) {
    patristic_matrix <- patristic_matrix[rowSums(is.na(patristic_matrix)) != ncol(patristic_matrix), colSums(is.na(patristic_matrix)) != nrow(patristic_matrix)]
  } # Get rid of rows and columns with all NA or NaNs, leaves the ones with some NA or NaNs
  if (dim(patristic_matrix)[1] == 2) {
    phy <- ape::rtree(n = 2, rooted = TRUE, tip.label = rownames(patristic_matrix), br = 0.5 * patristic_matrix[1, 2])
    phy$clustering_method <- "ape::rtree"
    phy$citation <- names(patristic_matrix)
    return(phy)
  }
  phycluster <- cluster_patristicmatrix(patristic_matrix, variance_matrix)
  phy <- choose_cluster(phycluster, clustering_method)
  if (!inherits(phy, "phylo")) {
    message("Clustering patristic matrix to phylo failed.")
    phy$citation <- names(patristic_matrix)
    return(phy)
  }
  phy$negative_brlen <- NA
  mess1 <- "Converting from patristic distance matrix to a tree resulted in some negative branch lengths;"
  mess2 <- "the largest by magnitude is"
  # when original tree IS ultrametric and has negative edges:
  if (is.null(phy$edge.length.original) & any(phy$edge.length < 0)) {
    warning(paste(mess1, mess2, min(phy$edge.length)))
  }
  # when original tree is NOT ultrametric and has negative edges:
  if (!is.null(phy$edge.length.original) & any(phy$edge.length.original < 0)) {
    warning(paste(mess1, mess2, min(phy$edge.length.original)))
    # and when tree was forced ultrametric and there still are neg edges:
    if (any(phy$edge.length < 0)) {
      warning(paste("After phytools::forcing.ultrametric there still are negative branch lengths;", mess2, min(phy$edge.length)))
    }
  }
  if (any(phy$edge.length < 0)) {
    phy$negative_brlen <- list(edge_number = which(phy$edge.length < 0))
    # phy$edge.length[which(phy$edge.length < 0)] <- 0 #sometimes NJ returns tiny negative branch lengths. https://github.com/phylotastic/datelife/issues/11
    if (fix_negative_brlen) {
      phy$negative_brlen <- list(edge_number = which(phy$edge.length < 0))
      phy <- tree_fix_brlen(tree = phy, fixing_criterion = "negative", fixing_method = fixing_method)
      fixing_method_called <- as.list(environment())$fixing_method
      phy$negative_brlen <- c(phy$negative_brlen, list(fixing_method = fixing_method_called))
      warning(paste0("Negative branch lengths were fixed with tree_fix_brlen, fixing_method = ", fixing_method_called))
    }
  }
  # for cases when there are no neg branch lengths to fix (or we don't want them fixed)
  # and we still want the final tree to be ultrametric:
  if (ultrametric) {
    if (is.null(phy$edge.length.original)) {
      phy <- force_ultrametric(phy)
    }
  }
  phy$tip.label <- gsub(" ", "_", phy$tip.label)
  phy$citation <- names(patristic_matrix)
  class(phy) <- c(class(phy), "datelifeTree")
  return(phy)
}
#' Force a non-ultrametric `phylo` object to be ultrametric with [phytools::force.ultrametric()].
#' @inheritParams phylo_check
#' @return A `phylo` object.
#' @export
force_ultrametric <- function(phy) {
  # TODO: check if there is an edge.length.original already there
  # something like how many grepl("edge.length.original") in names(phy) and add an integer after it.
  if (!inherits(phy, "phylo")) {
    message("phy argument is not a phylo object.")
    return(NA)
  }
  if (!ape::is.ultrametric(phy)) {
    phy$edge.length.original <- phy$edge.length
    phy <- phytools::force.ultrametric(tree = phy, method = "extend")
    phy$force.ultrametric <- "extend"
  }
  return(phy)
}
#' Cluster a patristic matrix into a tree with various methods.
#'
#' @inheritParams patristic_matrix_to_phylo
#' @return A list of trees obtained with clustering methods detailed in [patristic_matrix_to_phylo()].
#' @details If clustering method fails, `NA` is returned.
#' @export
cluster_patristicmatrix <- function(patristic_matrix, variance_matrix = NULL) {
  if (!inherits(patristic_matrix, "matrix") & !inherits(patristic_matrix, "data.frame")) {
    message("patristic_matrix argument is not a matrix")
    return(NA)
  }
  # has to be matrix not data frame:
  if (inherits(patristic_matrix, "data.frame")) {
    patristic_matrix <- as.matrix(patristic_matrix)
  }
  if (dim(patristic_matrix)[1] < 2) {
    return(NA)
  }
  if (dim(patristic_matrix)[1] == 2) {
    message("patristic_matrix has two taxa only, you don't need to cluster.")
    return(NA)
  } else {
    phyclust <- vector(mode = "list", length = 9)
    names(phyclust) <- c("nj", "njs", "upgma", "upgma_daisy", "bionj", "bionjs", "triangMtd", "triangMtds", "mvrs")
    phyclust$nj <- tryCatch(ape::nj(patristic_matrix), error = function(e) NA)
    if (inherits(phyclust$nj, "phylo")) {
      phyclust$nj <- tryCatch(phytools::midpoint.root(phyclust$nj),
        error = function(e) NA
      )
    }
    # if (is.null(phyclust$nj)){ # case when we have missing data (NA) on patristic_matrix and regular nj does not work; e.g. clade thraupidae SDM.results have missing data, and upgma chokes
    # njs appears to be the only option for missing data with NJ
    # but see Criscuolo and Gascuel. 2008. Fast NJ-like algorithms to deal with incomplete distance matrices. BMC Bioinformatics 9:166
    phyclust$njs <- tryCatch(ape::njs(patristic_matrix), error = function(e) NA)
    if (inherits(phyclust$njs, "phylo")) {
      phyclust$njs <- tryCatch(phytools::midpoint.root(phyclust$njs),
        error = function(e) NA
      )
    }
    # } else {
    # root the tree on the midpoint (only for trees with ape::Ntip(phy) > 2):
    # phy <- tryCatch(phangorn::midpoint(phy), error = function(e) NULL)
    # using phytools::midpoint.root instead of phangorn::midpoint bc it's less error prone.
    # sometimes, nj and njs do not work if patristic matrices come from sdm. why? it was probably the midpoint function from phangorn. Using phytools one now.
    # use regular upgma (or implemented with daisy and hclust) when nj or midpoint.root fail
    phyclust$upgma <- tryCatch(phangorn::upgma(patristic_matrix), error = function(e) NA)
    # if (is.null(phyclust$upgma)){ # case when we have missing data (NA) on patristic_matrix and regular upgma does not work; e.g. clade thraupidae SDM.results have missing data, and upgma chokes
    phyclust$upgma_daisy <- tryCatch(
      {
        # using daisy to calculate dissimilarity matrix instead of as.dist (which is used in phangorn::upgma) when there are NAs in the matrix; agnes does not work with NAs either.
        patristic_matrix <- patristic_matrix * 0.5 # doing this because it's giving ages that are too old, so it must be taking total distance
        DD <- cluster::daisy(x = patristic_matrix, metric = "euclidean")
        hc <- stats::hclust(DD, method = "average") # original clustering method from phangorn::upgma. Using agnes() instead hclust() to cluster gives the same result.
        phy <- ape::as.phylo(hc)
        phy <- phylobase::reorder(phy, "postorder")
        phy
      },
      error = function(e) NA
    )
    # }
    phyclust$bionj <- tryCatch(ape::bionj(patristic_matrix), error = function(e) NA)
    # if (is.null(phyclust$bionj)){ # case when we have missing data (NA) on patristic_matrix and regular nj does not work; e.g. clade thraupidae SDM.results have missing data, and upgma chokes
    # njs appears to be the only option for missing data with NJ
    # but see Criscuolo and Gascuel. 2008. Fast NJ-like algorithms to deal with incomplete distance matrices. BMC Bioinformatics 9:166
    phyclust$bionjs <- tryCatch(ape::bionjs(patristic_matrix), error = function(e) NA)
    if (inherits(phyclust$bionjs, "phylo")) {
      phyclust$bionjs <- tryCatch(phytools::midpoint.root(phyclust$bionjs),
        error = function(e) NA
      )
    }
    # } else {
    if (inherits(phyclust$bionj, "phylo")) {
      phyclust$bionj <- tryCatch(phytools::midpoint.root(phyclust$bionj),
        error = function(e) NA
      )
    }
    phyclust$triangMtd <- tryCatch(ape::triangMtd(patristic_matrix), error = function(e) NA)
    # if (is.null(phyclust$triangMtd)){ # case when we have missing data (NA) on patristic_matrix and regular nj does not work; e.g. clade thraupidae SDM.results have missing data, and upgma chokes
    # njs appears to be the only option for missing data with NJ
    # but see Criscuolo and Gascuel. 2008. Fast NJ-like algorithms to deal with incomplete distance matrices. BMC Bioinformatics 9:166
    phyclust$triangMtds <- tryCatch(ape::triangMtds(patristic_matrix), error = function(e) NA)
    if (inherits(phyclust$triangMtds, "phylo")) {
      phyclust$triangMtds <- tryCatch(phytools::midpoint.root(phyclust$triangMtds),
        error = function(e) NA
      )
    }
    # } else {
    if (inherits(phyclust$triangMtd, "phylo")) {
      phyclust$triangMtd <- tryCatch(phytools::midpoint.root(phyclust$triangMtd),
        error = function(e) NA
      )
    }
    if (inherits(variance_matrix, "matrix")) {
      # not possible to use the version for complete matrices, how to fill a variance matrix with missing values? I tried filling it with 0s and it runs but the output trees are network like...
      phyclust$mvrs <- tryCatch(ape::mvrs(patristic_matrix, variance_matrix), error = function(e) NA)
      if (inherits(phyclust$mvrs, "phylo")) {
        if (any(is.na(phyclust$mvrs$edge.length))) {
          phyclust$mvrs <- NA
        }
      }
    }
    return(phyclust)
  }
}
#' Choose an ultrametric phylo object from [cluster_patristicmatrix()] obtained
#' with a particular clustering method, or the next best tree.
#' If there are no ultrametric trees, it does not force them to be ultrametric.
#'
#' @inheritParams patristic_matrix_to_phylo
#' @param phycluster An output from [cluster_patristicmatrix()]
#' @return A `phylo` object or `NA`.
#' @export
choose_cluster <- function(phycluster, clustering_method = "nj") {
  if (!mode(phycluster) %in% "list") {
    message("phycluster argument is not a list; check that out.")
    return(NA)
  }
  # Choose between nj, njs, upgma or upgma_daisy only for now.
  # keep <- match(c("nj", "njs", "upgma", "upgma_daisy"), names(phycluster))
  # phycluster <- phycluster[keep]
  phy_return <- NA
  if (length(phycluster) == 0) {
    message("phycluster argument is length 0")
    return(NA)
  }
  if (inherits(phycluster, "phylo")) { # it is a tree of two tips
    return(phycluster)
  } else { # it is a list of results from cluster_patristicmatrix
    fail <- sapply(phycluster, function(x) !inherits(x, "phylo"))
    if (all(fail)) {
      message("The patristic matrix could not be transformed into a tree with any of the default methods (NJ, UPGMA)")
      return(NA)
    }
    phycluster <- phycluster[!fail] # take out the fails or unattempted from cluster_patristicmatrix
    if (length(phycluster) == 1) {
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
      if (length(ultram) == 0 & length(ultram2) == 0) {
        message(paste("The patristic matrix could not be transformed into an ultrametric tree with any of the default methods:", toupper(names(phycluster))))
        # return(NA)
      }
      choice <- grepl(clustering_method, names(phycluster)) # choice can only be one
      ff <- which(choice & ultram) # if the chosen method gives an ultrametric tree
      if (length(ff) != 0) {
        ff <- ff[1]
        phy <- phycluster[[ff]]
        phy$clustering_method <- names(phycluster)[ff]
        return(phy)
      }
      ff <- which(!choice & ultram) # if not, take the not chosen but ultrametric
      if (length(ff) != 0) {
        ff <- ff[1]
        phy <- phycluster[[ff]]
        phy$clustering_method <- names(phycluster)[ff]
        return(phy)
      }
      ff <- which(choice & ultram2) # if not, take the chosen one but less ultrametric
      if (length(ff) != 0) {
        ff <- ff[1]
        phy <- phycluster[[ff]]
        # phy$edge.length.original <- phy$edge.length
        # phy <- phytools::force.ultrametric(tree = phy, method = "extend")
        # phy$force.ultrametric <- "extend"
        phy$clustering_method <- names(phycluster)[ff]
        return(phy)
      }
      ff <- which(!choice & ultram2) # if not, take the not chosen one but less ultrametric
      if (length(ff) != 1) {
        ff <- ff[1]
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

#' Go from a summary matrix to an ultrametric `phylo` object.
#' @param summ_matrix Any summary patristic distance matrix, such as the ones obtained with [datelife_result_sdm_matrix()] or [datelife_result_median_matrix()].
#' @inheritParams get_taxon_summary
#' @param total_distance Whether the input `summ_matrix` stores total age distance
#'  (from tip to tip) or distance from node to tip. Default to `TRUE`,
#'  divides the matrix in half, if `FALSE` it will take it as is.
#' @param use A character vector indicating what type of age to use for summary tree.
#'  One of the following:
#' \describe{
#' 	\item{"mean"}{It will use the [mean()] of the node ages in `summ_matrix`.}
#' 	\item{"median"}{It uses the [stats::median()] age of node ages in `summ_matrix`.}
#' 	\item{"min"}{It will use the [min()] age from node ages in `summ_matrix`.}
#' 	\item{"max"}{Choose this if you wanna be conservative; it will use the [max()]
#'        age from node ages in `summ_matrix`.}
#' 	\item{"midpoint"}{It will use the mean of minimum age and maximum age.}
#' }
#' @param target_tree A `phylo` object. Use this in case you want a specific
#'  backbone for the output tree.
#' @inheritDotParams summary_matrix_to_phylo_all
#' @return An ultrametric phylo object.
#' @details It can take a regular patristic distance matrix, but there are simpler
#'  methods for that implemented in [patristic_matrix_to_phylo()].
#' @export
# #' @examples
# #' summary_matrix_to_phylo(summ_matrix, use = "median")
summary_matrix_to_phylo <- function(summ_matrix,
                                    datelife_query = NULL,
                                    target_tree = NULL,
                                    total_distance = TRUE,
                                    use = "mean",
                                    ...) {
  # enhance: add other dating methods, besides bladj.
  use <- match.arg(use, c("midpoint", "mean", "median", "min", "max"))

  all_trees <- summary_matrix_to_phylo_all(summ_matrix = summ_matrix,
                                           datelife_query = datelife_query,
                                           target_tree = target_tree,
                                           total_distance = total_distance)
  ##############################################################################
  ##############################################################################
  # add info to return object, and return
  ############################################################################
  ############################################################################
  if ("min" %in% use) {
    new_phy <- all_trees$min
  }
  if ("max" %in% use) {
    new_phy <- all_trees$max
  }
  if ("mean" %in% use) {
    new_phy <- all_trees$mean
  }
  if ("median" %in% use) {
    new_phy <- all_trees$median
  }
  if ("midpoint" %in% use) {
    new_phy <- all_trees$midpoint
  }
  new_phy$calibration_distribution <- attributes(all_trees)$node_age_distributions
  if (is.null(new_phy$clustering_method)) {
    new_phy$clustering_method <- NULL
  }
  if (inherits(target_tree, "phylo")) {
    if (!is.null(target_tree$ott_ids) & is.null(new_phy$ott_ids)) {
      tt <- match(new_phy$tip.label, target_tree$tip.label)
      # match(c("a", "b", "c", "d"), c("c", "d", "a", "a", "a", "b"))
      new_phy$ott_ids <- target_tree$ott_ids[tt]
    }
  }
  return(new_phy)
}
