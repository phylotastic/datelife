# uncertainty generation Functions

#' Generate uncertainty in branch lengths using a lognormal.
#' @inheritParams phylo_check
#' @inheritParams sample_trees
#' @param uncertainty_method A character vector specifying the method to generate uncertainty. mrbayes is default.
#' @inheritParams make_mrbayes_tree
#' @param age_sd The standard deviation around the age to use for generating the uncertainty. If not a numeric value, var will be used to calculate it.
#' @param age_var The variance to calculate age_sd and generate uncertainty.
#' @param age_scale How to scale sd by the depth of the node. If 0, same sd for all. If not, older nodes have more uncertainty
#' @param alpha The significance level on uncertainty to generate. By default 0.025
#' @param rescale Boolean. If true, observed age will be rescaled each round.
#' @inheritParams datelife_search
#' @return A phylo or multiPhylo object with the same topology as phy but different branch lengths
#' @export
#' @details If you want to change the size of sampled trees you do not need to run mrbayes again.
#' Just use sample_trees("mrbayes_trees_file_directory", size = new_size) and you will get a multiPhylo object with a new tree sample.
#' @examples
#' # Generate uncertainty over feline species SDM chronogram.
#' # Load the data:
#'
#' data(felid_sdm)
#'
#' # By default, generates a sample of 100 trees with var = 0.1:
#'
#' unc <- phylo_generate_uncertainty(felid_sdm$phy)
#' length(unc)
#'
#' # Make an LTT plot:
#'
#' max_age <- max(sapply(unc, ape::branching.times))
#' ape::ltt.plot(phy = unc[[1]], xlim = c(-max_age, 0), col = "#cce5ff50")
#' for (i in 2:100) {
#'   ape::ltt.lines(phy = unc[[i]], col = "#cce5ff50")
#' }
#' ape::ltt.lines(felid_sdm$phy, col = "red")
#' title(c("fake uncertainty", "in Felidae SDM chronogram"))
phylo_generate_uncertainty <- function(phy, size = 100, uncertainty_method = "other", age_distribution = "uniform", age_sd = NULL, age_var = 0.1, age_scale = 0, alpha = 0.025, rescale = TRUE) {
  phylo_check(phy)
  size <- round(as.numeric(size), digits = 0)
  uncertainty_method <- match.arg(uncertainty_method, c("mrbayes", "other"))

  # if(uncertainty_method == "mrbayes"){
  #   phy_name <- gsub("[[:punct:]]", "_", deparse(substitute(phy)))  # gets the name of the phy object to name mrbayes files
  #   mrbayes_file <- paste0(phy_name, "_mrbayes_uncertainty_", age_distribution, ".nexus")
  #   # print(mrbayes_file)
  #   phy <- tree_add_outgroup(tree = phy, outgroup = "fake_outgroup")
  #   phycontre <- make_mrbayes_tree(constraint = phy, ncalibration = phy, age_distribution = age_distribution, mrbayes_output_file = mrbayes_file)
  #   phycontre <- ape::drop.tip(phycontre, "fake_outgroup")
  #   phycloud <- sample_trees(trees_file = paste0(mrbayes_file, ".t"), burnin = 0.25, size = size)
  #   phycloud <- lapply(phycloud, ape::drop.tip, "fake_outgroup")  # ape::drop.tip does not have a multiPhylo option, so I used lapply
  #   if (length(phycloud) == 1){
  #     phycloud <- phycloud[[1]]
  #     class(phycloud) <- "phylo"
  #   } else {
  #     class(phycloud) <- "multiPhylo"
  #   }
  #   return(list(consensus_tree = phycontre, trees = phycloud))
  # }
  if (uncertainty_method == "other") {
    phy.new <- phy <- ape::reorder.phylo(phy, "postorder")
    phy_depths <- max(ape::branching.times(phy)) - phytools::nodeHeights(phy)
    # phy_depths.new <- phy_depths
    if (is.numeric(age_sd[1])) { # n is number of random deviates, age_mu is the observed age
      # my_rlnorm <- function(n, age_mu, variance = age_sd^2){
      #   res <- stats::rlnorm(n = n, meanlog = log(age_mu / sqrt(1 + (variance / age_mu^2))), sdlog = sqrt(log(1 + variance / age_mu^2)))
      #   return(res)
      # }
      my_rlnorm <- function(n, age_mu) {
        res <- stats::rlnorm(n = n, meanlog = log(age_mu / sqrt(1 + (age_sd^2 / age_mu^2))), sdlog = age_sd)
        return(res)
      }
    } else if (is.numeric(age_var[1])) {
      my_rlnorm <- function(n, age_mu) {
        res <- stats::rlnorm(n = n, meanlog = log(age_mu / sqrt(1 + (age_var / age_mu^2))), sdlog = sqrt(log(1 + age_var / age_mu^2)))
        return(res)
      }
    } else {
      stop("age_sd or age_var argument must be a numeric vector")
    }
    # consider a confidence interval alpha:
    # if (is.numeric(alpha[1])){
    #   age_sd <- ((log(sd_amount[1] * ape::branching.times(phy))) - ape::branching.times(phy)) / (-2*log(sqrt(2*pi) * alpha[1]))
    #   age_sd <- age_sd + age_sd * sd_depth * ape::branching.times(phy)
    # }
    nn <- sort(phylo_get_node_numbers(phy))
    res <- c()
    while (length(res) < size) {
      message("Uncertainty sample number ", length(res) + 1)

      phy_depths.final <- phy_depths.original <- phy_depths
      tot_age <- my_rlnorm(1, age_mu = ape::branching.times(phy)[as.character(nn[1])]) # variance is determined by the function above, depending on sd or var
      phy_depths.final[phy$edge == nn[1]] <- tot_age
      for (i in nn[-1]) {
        max_age <- nn_age <- phy_depths.final[phy$edge[, 2] == i, 1]
        tries <- 0 # by chance, sampled age can be older than that of parent node; resample 50 times to try to get a younger age.
        while (max_age - nn_age <= 0 & tries < 50) {
          if (rescale) {
            mu <- ape::branching.times(phy)[as.character(i)] * max_age / phy_depths.original[phy$edge[, 2] == i, 1]
            nn_age <- my_rlnorm(1, age_mu = mu)
            nn_age <- nn_age * ape::branching.times(phy)[as.character(i)] / mu
          } else {
            mu <- ape::branching.times(phy)[as.character(i)]
            nn_age <- my_rlnorm(1, age_mu = mu)
          }
          tries <- tries + 1
        }
        phy_depths.final[phy$edge == i] <- nn_age # if sampled age is younger than parent node age, use it
      }
      phy.new$edge.length <- phy_depths.final[, 1] - phy_depths.final[, 2]
      if (all(phy.new$edge.length > 0)) {
        res <- c(res, list(phy.new))
      }
    }
    # if (uncertainty_method=="poisson") {
    #   # internal.variances <- rep(NA,asd)
    # }
    if (length(res) == 1) {
      res <- res[[1]]
      class(res) <- "phylo"
    } else {
      class(res) <- "multiPhylo"
    }
    return(res)
  }
}

#' Sample trees from a file containing multiple trees. Usually from a bayesian analysis output trees file.
#' @param trees_file A character vector indicating the name and directory of file with trees to sample.
#' @param trees_object An R object containing a list of trees already read into R from a tree file from a bayesian analysis output.
#' @param burnin A numeric vector indicating the burnin fraction. It should be a number between 0 and 1. Default to 0.25
#' @param size A numeric vector indicating the number of samples to be generated.
#' @return A `multiPhylo` object with a random sample of trees.
#' @export

sample_trees <- function(trees_file, trees_object = NULL, burnin = 0.25, size = 100) {
  if (!inherits(trees_object, "list")) {
    phycloud <- ape::read.nexus(trees_file)
  } else {
    phycloud <- trees_object
  }
  # use ceiling instead of round to make sure 0 is never among the values to sample:
  phycloud <- phycloud[sample(ceiling(burnin * length(phycloud)):length(phycloud), size)] # sample size number of trees from the cloud of trees, after discarding 25% as burnin
  if (length(phycloud) == 1) {
    phycloud <- phycloud[[1]]
    class(phycloud) <- "phylo"
  } else {
    class(phycloud) <- "multiPhylo"
  }
  return(phycloud)
}
