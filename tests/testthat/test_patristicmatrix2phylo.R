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

test_that("patristic_matrix_to_phylo gives ultrametric trees", {
  dq <- make_datelife_query(input = "cetacea", get_spp_from_taxon = TRUE)
  dr <- get_datelife_result(dq)
  t0 <- summarize_datelife_result(dq, dr, summary_format = "phylo_sdm")
  t1 <- datelife_result_sdm(dr, clustering_method = "nj")
  t2 <- datelife_result_sdm(dr, clustering_method = "upgma")
  tt <- list(t0, t1$phy, t2$phy)
  class(tt) <- "multiPhylo"
  ape::is.ultrametric(tt, 2)
  # plot(tt, cex = 0.5)
})

# test_that("We get species names back from subspecies and var names",{
#   # rotl::tnrs_match_names("Ceratopogon slossonae")
#   NA
# })
test_that("cetacea sdm matrix to phylo is accurate", {
  datelife_result <- get_datelife_result(input = "cetacea")
  unpadded.matrices <- lapply(datelife_result, patristic_matrix_unpad)
  good.matrix.indices <- get_goodmatrices(unpadded.matrices, verbose = TRUE)
  if(length(good.matrix.indices) > 1) {
    unpadded.matrices <- unpadded.matrices[good.matrix.indices]
    sdm_matrix <- get_sdm(unpadded.matrices, weighting = "flat", verbose = TRUE)
  }
  max(sdm_matrix, na.rm = TRUE)/2
  sdm_phylo <- summarize_datelife_result(datelife_result = datelife_result, summary_format = "phylo_sdm")
  max(ape::branching.times(sdm_phylo))
  phycluster <- cluster_patristicmatrix(patristic_matrix = sdm_matrix)
  sapply(phycluster, function(x) max(ape::branching.times(x)))
  plot(phycluster$njs, cex=0.2)
  ape::axisPhylo()
  plot(phycluster$upgma_daisy, cex=0.2)
  ape::axisPhylo()
  phy <- choose_cluster(phycluster, clustering_method = "nj")
  ape::ltt.plot(phycluster$upgma_daisy, col = "blue")
  ape::ltt.lines(force_ultrametric(phycluster$njs))
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
