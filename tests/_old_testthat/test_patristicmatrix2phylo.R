test_that("felis/canidae divergence is accurate", {
    skip_on_cran()
    skip_on_travis()
    matrix_max_ages <- sapply(catsanddogs_results, max)
    taxa <- c("felis", "canidae")
    cats_and_dogs <- datelife_search(input = taxa, get_spp_from_taxon = TRUE,
      summary_format = "phylo_all")
    phylo_max_ages <- sapply(cats_and_dogs, function(x) max(ape::branching.times(x)))
    expect_true(all(names(matrix_max_ages) == names(phylo_max_ages)))
    # names(matrix_max_ages) <- names(phylo_max_ages)<- NULL
    ns <- 20
    xx <- format(round(sort(matrix_max_ages/2)), nsmall = ns) == format(round(sort(phylo_max_ages)), nsmall = ns)
    # ages from our cache range from 54.9 to 70.9, this includes upper limit confidence interval chronograms
    # timetree study-derived ages range from 39.7 to 67.1. This excludes confidence intervals
    median_phy <- summarize_datelife_result(datelife_result = catsanddogs_results,
        datelife_query = catsanddogs_query, summary_format = "phylo_median")
    sdm_phy <- summarize_datelife_result(datelife_result = catsanddogs_results,
        datelife_query = catsanddogs_query, summary_format = "phylo_sdm")
})

test_that("patristic_matrix_to_phylo works: cetaceae", {
  unpadded.matrices <- lapply(cetacea_result, patristic_matrix_unpad)
  good.matrix.indices <- get_goodmatrices(unpadded.matrices, verbose = TRUE)
  if(length(good.matrix.indices) > 1) {
    unpadded.matrices <- unpadded.matrices[good.matrix.indices]
    sdm_matrix <- get_sdm(unpadded.matrices, weighting = "flat", verbose = TRUE)
  }
  # max(sdm_matrix, na.rm = TRUE)/2
  t0 <- summarize_datelife_result(datelife_result = cetacea_result, summary_format = "phylo_sdm")
  max(ape::branching.times(t0))
  expect_true(ape::is.ultrametric(t0, 1))
  expect_true(ape::is.ultrametric(t0, 2))
})
test_that("patristic_matrix_to_phylo works: subset2", {
  p1 <- patristic_matrix_to_phylo(subset2_sdm_matrix, clustering_method = "nj",
        fix_negative_brlen = TRUE, fixing_method = 0, ultrametric = TRUE)
  p2 <- patristic_matrix_to_phylo(subset2_sdm_matrix, clustering_method = "upgma",
        fix_negative_brlen = TRUE, fixing_method = 0, ultrametric = TRUE)
  # plot(p2, cex = 0.5)
  expect_true(inherits(p1, "phylo"))
  expect_true(inherits(p2, "phylo"))
})
test_that("patristic_matrix_to_phylo works: some ants", {
      xx <- patristic_matrix_to_phylo(patristic_matrix = some_ants_datelife_result[[1]])
      expect_s3_class(xx, "phylo")
      expect_true(ape::is.ultrametric(xx))
      #make sure it works with missing data:
      withNaN <- some_ants_datelife_result[[1]]
      withNaN[9, 8] <- NaN
      withNaN[8, 9] <- NaN
      xx <- patristic_matrix_to_phylo(patristic_matrix = withNaN)
      expect_s3_class(xx, "phylo")
      expect_true(ape::is.ultrametric(xx))# because NaN, it uses njs, and gives a tree
})

test_that("cluster_patristicmatrix works", {
    c1 <- cluster_patristicmatrix(subset2_sdm_matrix)
    expect_true(mode(c1) %in% "list")
})

test_that("summary_matrix_to_phylo works: threebirds", {
  unpadded.matrices <- lapply(threebirds_result, patristic_matrix_unpad)
  good.matrix.indices <- get_goodmatrices(unpadded.matrices, verbose = TRUE)
  if(length(good.matrix.indices) > 1) {
    unpadded.matrices <- unpadded.matrices[good.matrix.indices]
    sdm_matrix <- get_sdm(unpadded.matrices, weighting = "flat", verbose = TRUE)
  }
  # max(sdm_matrix, na.rm = TRUE)/2
  sdm_phylo <- summary_matrix_to_phylo(summ_matrix = sdm_matrix)
  # names(sdm_phylo)
  expect_true(inherits(sdm_phylo, "phylo"))
  expect_true(ape::is.ultrametric(sdm_phylo, 2))
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
