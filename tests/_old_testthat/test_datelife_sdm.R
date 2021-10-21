test_that("SDM correctly returns three birds tree", {
  skip_on_cran()
  skip_on_os("linux") #b/c no pathd8 on travis linux
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
test_that("datelife_result_sdm works with all clustering methods", {
    t0 <- datelife_result_sdm(cetacea_result)
    # t1 <- datelife_result_sdm(cetacea_result, clustering_method = "nj")
    # t2 <- datelife_result_sdm(cetacea_result, clustering_method = "upgma")
})
