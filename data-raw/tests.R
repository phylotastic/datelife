if(dim(patristic_matrix)[1] == 2) {
    clum <- "ape::rtree"
    phy <- ape::rtree(n = 2, rooted = TRUE, tip.label = rownames(patristic_matrix), br = 0.5 * patristic_matrix[1,2])
} else {
    if (clustering_method == "nj"){
      tried <- c(tried, "nj")
      clum <- "nj"
      phy <- tryCatch(ape::nj(patristic_matrix), error = function (e) NULL)
      if (is.null(phy)){ # case when we have missing data (NA) on patristic_matrix and regular nj does not work; e.g. clade thraupidae SDM.results have missing data, and upgma chokes
          tried <- c(tried, "njs")
          clum <- "njs"
          # njs appears to be the only option for missing data with NJ
          # but see Criscuolo and Gascuel. 2008. Fast NJ-like algorithms to deal with incomplete distance matrices. BMC Bioinformatics 9:166
          phy <- tryCatch(ape::njs(patristic_matrix),
          error = function(e) {
            # message("unable to cluster patristic matrix with NJ, trying clustering_method = 'upgma'")
            return(NULL)
          })
      }
      if(ape::Ntip(phy) > 2) {
          # root the tree on the midpoint:
          # phy <- tryCatch(phangorn::midpoint(phy), error = function(e) NULL)
          # using phytools::midpoint.root instead of phangorn::midpoint bc it's less error prone.
          phy <- tryCatch(phytools::midpoint.root(phy),
          error = function(e) {
            # message("unable to root the tree after NJ, trying clustering_method = 'upgma'")
            return(NULL)
          })
      }
    }
    # use regular upgma (or implemented with daisy and hclust) when nj or midpoint.root fail
    # sometimes, nj and njs do not work if patristic matrices come from sdm. why? it was probably the midpoint function from phangorn. Using phytools one now.
    if (clustering_method == "upgma"){ # is.null(phy) |
        # tried <- c(tried, "upgma")
        clum <- "upgma"
        phy <- tryCatch(phangorn::upgma(patristic_matrix), error = function (e) NULL)
        if (is.null(phy)){ # case when we have missing data (NA) on patristic_matrix and regular upgma does not work; e.g. clade thraupidae SDM.results have missing data, and upgma chokes
            phy <- tryCatch({
              # tried <- c(tried, "upgma_daisy")
              clum <- "upgma_daisy"
              # using daisy to calculate dissimilarity matrix instead of as.dist (which is used in phangorn::upgma) when there are NAs in the matrix; agnes does not work with NAs either.
              DD <- cluster::daisy(x = patristic_matrix, metric = "euclidean")
              hc <- stats::hclust(DD, method = "average") # original clustering method from phangorn::upgma. Using agnes() instead hclust() to cluster gives the same result.
              phy <- ape::as.phylo(hc)
              phy <- phylobase::reorder(phy, "postorder")
              phy
            }, error = function(e) NULL)
        }
        # if(is.null(phy) & !"nj" %in% tried){ # Trying the exact same chunck for nj was not tried already
        #   clum <- "nj"
        #   phy <- tryCatch(ape::nj(patristic_matrix), error = function (e) NULL)
        #   if (is.null(phy)){
        #       clum <- "njs"
        #       phy <- tryCatch(ape::njs(patristic_matrix),
        #       error = function(e) {
        #         return(NULL)
        #       })
        #   }
        # }
    }
}
