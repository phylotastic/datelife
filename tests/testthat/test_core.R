test_that("ProcessResultsList", {
data(opentree_chronograms)
datelife.cache <- results
taxa <- c("Rhea americana", "Pterocnemia pennata", "Struthio camelus")
results.list <- lapply(datelife.cache$trees,GetSubsetArrayDispatch, taxa=taxa, phy=NULL)

expect_gte(length(ProcessResultsList(results.list, taxa), 4)
})
