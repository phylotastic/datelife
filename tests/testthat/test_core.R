test_that("ProcessResultsList", {
data(opentree_chronograms)
taxa <- c("Rhea americana", "Pterocnemia pennata", "Struthio camelus")
results.list <- lapply(datelife.cache$trees,GetSubsetArrayDispatch, taxa=taxa, phy=NULL)
expect_gte(length(ProcessResultsList(results.list, taxa, TRUE)), 4)
})
