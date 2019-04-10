test_that("felis/canidae divergence is accurate", {
    query <- make_datelife_query(input = c("felis", "canidae"), get_spp_from_taxon = TRUE)
    cats_and_dogs_results <- get_datelife_result(input = query)
    matrix_max_ages <- sapply(cats_and_dogs_results, max)
    taxa <- c("felis", "canidae")
    cats_and_dogs <- datelife_search(input = taxa, get_spp_from_taxon = TRUE,
      summary_format = "phylo_all")
    phylo_max_ages <- sapply(cats_and_dogs, function(x) max(ape::branching.times(x)))
    expect_true(all(names(matrix_max_ages) == names(phylo_max_ages)))
    # names(matrix_max_ages) <- names(phylo_max_ages)<- NULL
    ns <- 20
    format(round(sort(matrix_max_ages/2)), nsmall = ns) == format(round(sort(phylo_max_ages)), nsmall = ns)
    # ages from our cache range from 54.9 to 70.9, this includes upper limit confidence interval chronograms
    # timetree study-derived ages range from 39.7 to 67.1. This excludes confidence intervals
    median_phy <- summarize_datelife_result(datelife_result = cats_and_dogs_results, datelife_query = query,
      summary_format = "phylo_median")
    sdm_phy <- summarize_datelife_result(datelife_result = cats_and_dogs_results, datelife_query = query,
        summary_format = "phylo_sdm")
})

test_that("patristic_matrix_to_phylo works", {
  unpadded.matrices <- lapply(cetacea_result, patristic_matrix_unpad)
  good.matrix.indices <- get_goodmatrices(unpadded.matrices, verbose = TRUE)
  if(length(good.matrix.indices) > 1) {
    unpadded.matrices <- unpadded.matrices[good.matrix.indices]
    sdm_matrix <- get_sdm(unpadded.matrices, weighting = "flat", verbose = TRUE)
  }
  max(sdm_matrix, na.rm = TRUE)/2
  t0 <- summarize_datelife_result(datelife_result = datelife_result, summary_format = "phylo_sdm")
  max(ape::branching.times(t0))
  ape::is.ultrametric(t0, 1)
  ape::is.ultrametric(t0, 2)
  p1 <- patristic_matrix_to_phylo(subset2_sdm_matrix, clustering_method = "nj",
        fix_negative_brlen = TRUE, fixing_method = 0, ultrametric = TRUE)
  p2 <- patristic_matrix_to_phylo(subset2_sdm_matrix, clustering_method = "upgma",
        fix_negative_brlen = TRUE, fixing_method = 0, ultrametric = TRUE)
  # plot(p2, cex = 0.5)
})

test_that("cluster_patristicmatrix works", {
    c1 <- cluster_patristicmatrix(subset2_sdm_matrix)
})

test_that("threebirds sdm matrix to phylo is accurate", {
  threebirds_dq <- make_datelife_query(input=c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"),
get_spp_from_taxon = TRUE)
  utils::data(threebirds_dr)
  unpadded.matrices <- lapply(threebirds_dr, patristic_matrix_unpad)
  good.matrix.indices <- get_goodmatrices(unpadded.matrices, verbose = TRUE)
  if(length(good.matrix.indices) > 1) {
    unpadded.matrices <- unpadded.matrices[good.matrix.indices]
    sdm_matrix <- get_sdm(unpadded.matrices, weighting = "flat", verbose = TRUE)
  }
  max(sdm_matrix, na.rm = TRUE)/2
  sdm_phylo <- summarize_datelife_result(datelife_result = threebirds_dr, summary_format = "phylo_sdm")
  names(sdm_phylo)
})

test_that("summary_matrix_to_phylo works with and without target trees", {
  subset2_sdmphylo_mean1 <- summary_matrix_to_phylo(summ_matrix = subset2_sdm_matrix,
            use = "mean")
  # subset2_sdm_matrix[,"Lycopodium annotinum"]
  # plot(subset2_sdmphylo_mean1, cex = 0.5)
  subset2_otol <- get_otol_synthetic_tree(colnames(subset2_sdm_matrix))
  # plot(subset2_otol, cex = 0.5)

  subset2_sdmphylo_mean2 <- summary_matrix_to_phylo(summ_matrix = subset2_sdm_matrix,
          use = "mean", target_tree = subset2_otol)
  # plot(subset2_sdmphylo_mean2, cex = 0.5)
})
