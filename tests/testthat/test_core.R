test_that("ProcessResultsList", {
data(opentree_chronograms)
taxa <- c("Rhea americana", "Pterocnemia pennata", "Struthio camelus")
results.list <- lapply(datelife.cache$trees,GetSubsetArrayDispatch, taxa=taxa, phy=NULL)
expect_gte(length(ProcessResultsList(results.list, taxa, TRUE)), 1)
})

test_that("Summarize as mrca works correctly", {
  data(opentree_chronograms)
  taxa <- c("Rhea americana", "Pterocnemia pennata", "Struthio camelus")
  results.list <- lapply(datelife.cache$trees,GetSubsetArrayDispatch, taxa=taxa, phy=NULL)
  filtered.results <- ProcessResultsList(results.list, taxa, TRUE)
  mrca.vector <- SummarizeResults(filtered.results, output.format="mrca", cache=datelife.cache)
  expect_equal(class(mrca.vector), "numeric")
  expect_gte(min(mrca.vector),5)
  expect_lte(max(mrca.vector),150)
})

test_that("Summarize as citations works correctly", {
  data(opentree_chronograms)
  taxa <- c("Rhea americana", "Pterocnemia pennata", "Struthio camelus")
  results.list <- lapply(datelife.cache$trees,GetSubsetArrayDispatch, taxa=taxa, phy=NULL)
  filtered.results <- ProcessResultsList(results.list, taxa, TRUE)
  citation.results <- SummarizeResults(filtered.results, output.format="cit", cache=datelife.cache)
  expect_equal(class(citation.results), "character")
  expect_gte(sum(grepl("Brown", citation.results)),1) #b/c any sort of list of bird phylogenies can't be complete without one authored by Joseph Brown
})

test_that("Summarize as newick.all works correctly", {
  data(opentree_chronograms)
  taxa <- c("Rhea americana", "Pterocnemia pennata", "Struthio camelus")
  results.list <- lapply(datelife.cache$trees,GetSubsetArrayDispatch, taxa=taxa, phy=NULL)
  filtered.results <- ProcessResultsList(results.list, taxa, TRUE)
  trees <- SummarizeResults(filtered.results, output.format="newick.all", cache=datelife.cache)
  expect_equal(class(trees), "character")
  expect_false(anyNA(trees))
  expect_equal(class(ape::read.tree(text=trees[1])), "phylo")
})

test_that("Summarize as newick.median works correctly", {
  data(opentree_chronograms)
  taxa <- c("Rhea americana", "Pterocnemia pennata", "Struthio camelus")
  results.list <- lapply(datelife.cache$trees,GetSubsetArrayDispatch, taxa=taxa, phy=NULL)
  filtered.results <- ProcessResultsList(results.list, taxa, partial=FALSE)
  tree <- SummarizeResults(filtered.results, output.format="newick.median", cache=datelife.cache)
  expect_equal(class(tree), "character")
  expect_false(anyNA(tree))
  expect_equal(class(ape::read.tree(text=tree)), "phylo")
})

test_that("Processing input newick", {
	skip_on_cran()
	skip_on_travis()
  data(opentree_chronograms)
  input.processed <- ProcessInput(ape::write.tree(rcoal(3, tip.label=c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"))), usetnrs=FALSE, approximatematch=TRUE)
  expect_equal(class(input.processed$phy)=="phylo")
})


test_that("Congruification works", {
	skip_on_cran()
	skip_on_travis()
  data(opentree_chronograms)
  filtered.results <- GetFilteredResults(ape::write.tree(rcoal(3, tip.label=c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"))), partial=TRUE, usetnrs=FALSE, approximatematch=TRUE, cache=datelife.cache)
  expect_gte(length(filtered.results, 4))
})

test_that("TNRS with approximate match works", {
	 taxa <- c("Rhea_americana", "Pterocnemia pennato", "Strutho camelus")
 	input.processed <- ProcessInput(taxa, usetnrs=TRUE, approximatematch=TRUE)
 	expect_true(all.equal(c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"), input.processed$cleaned.names))
})


test_that("TNRS with unmatchable taxa works", {
	taxa <- c("Rhea_americana", "Pterocnemia pennato", "Oscar the grouch", "Strutho camelus")
 	input.processed <- ProcessInput(taxa, usetnrs=TRUE, approximatematch=TRUE)
 	expect_true(all.equal(c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"), input.processed$cleaned.names))
})
