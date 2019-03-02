test_that("SDM correctly returns three birds tree", {
  skip_on_cran()
  skip_on_os("linux") #b/c no pathd8 on travis linux
  utils::data(opentree_chronograms)
  taxa <- c("Rhea americana", "Pterocnemia pennata", "Struthio camelus")
  results_list <- lapply(opentree_chronograms$trees, get_subset_array_dispatch, taxa=taxa, phy=NULL)
  datelife_result <- results_list_process(results_list, taxa, TRUE)
  result.tree <- datelife_result_sdm(datelife_result)
  expect_true(inherits(result.tree, "phylo"))
  ape::is.ultrametric(result.tree, options = 2)
})

test_that("subset2_taxa gives ultrametric summary trees", {
    # class(subset2_result)
    xx <- get_best_grove(subset2_result)
    # class(xx$best_grove)
    xx_median <- datelife_result_median(xx$best_grove)
    ape::is.ultrametric(xx_median, option = 2)
    # plot(xx_median, cex = 0.5, label.offset = 50)
    # abline(v= max(ape::branching.times(xx_median)))
    xx_sdm <- datelife_result_sdm(xx$best_grove)
    ape::is.ultrametric(xx_sdm, option = 2)
    # nrow(sdm.matrix)
})

# for testing sdm with cetacea:
# datelife_result <- get_datelife_result(input = "cetacea")
# unpadded.matrices <- lapply(datelife_result, patristic_matrix_unpad)
# good.matrix.indices <- get_goodmatrices(unpadded.matrices, verbose = TRUE)
# if(length(good.matrix.indices) > 1) {
#   unpadded.matrices <- unpadded.matrices[good.matrix.indices]
#   sdm_matrix <- get_sdm(unpadded.matrices, weighting = "flat", verbose = TRUE)
# }
# which(sdm_matrix < 0)
# sdm_matrix[ceiling(7301/ncol(sdm_matrix)),] # Eubalaena japonica,
# sdm_matrix[,ceiling(261/nrow(sdm_matrix))]  # Eubalaena glacialis
# xx <- sdm_matrix #[1:5, 1:5]
# even removing negative values for small positive values gives back non ultrametric trees with njs
# sdm_matrix[which(sdm_matrix < 0)] <- 0.01
# test <- cluster_patristicmatrix(sdm_matrix)
# class(test) <- "multiPhylo"
# ape::is.ultrametric(test)
# plot(test$njs)
