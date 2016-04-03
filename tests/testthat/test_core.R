test_that("ProcessResultsList", {
data(opentree_chronograms)
taxa <- c("Rhea americana", "Pterocnemia pennata", "Struthio camelus")
results.list <- lapply(datelife.cache$trees,GetSubsetArrayDispatch, taxa=taxa, phy=NULL)
expect_gte(length(ProcessResultsList(results.list, taxa, TRUE)), 4)
})

test_that("Summarize as mrca works correctly", {
  data(opentree_chronograms)
  taxa <- c("Rhea americana", "Pterocnemia pennata", "Struthio camelus")
  results.list <- lapply(datelife.cache$trees,GetSubsetArrayDispatch, taxa=taxa, phy=NULL)
  filtered.results <- ProcessResultsList(results.list, taxa, TRUE)
  mrca.vector <- SummarizeResults(filtered.results, output.format="mrca")
  expect_equal(class(mrca.vector), "numeric")
  expect_gte(min(mrca.vector),5)
  expect_lte(max(mrca.vector),150)
})

test_that("Summarize as citations works correctly", {
  data(opentree_chronograms)
  taxa <- c("Rhea americana", "Pterocnemia pennata", "Struthio camelus")
  results.list <- lapply(datelife.cache$trees,GetSubsetArrayDispatch, taxa=taxa, phy=NULL)
  filtered.results <- ProcessResultsList(results.list, taxa, TRUE)
  citation.results <- SummarizeResults(filtered.results, output.format="cit")
  expect_equal(class(citation.results), "character")
  expect_gte(sum(grepl("Brown", citation.results)),1) #b/c any sort of list of bird phylogenies can't be complete without one authored by Joseph Brown
  expect_gte(sum(grepl("Hedges", citation.results)),1) #the TimeTree phylogeny is an important one for this set, too
})

test_that("Summarize as newick works correctly", {
  data(opentree_chronograms)
  taxa <- c("Rhea americana", "Pterocnemia pennata", "Struthio camelus")
  results.list <- lapply(datelife.cache$trees,GetSubsetArrayDispatch, taxa=taxa, phy=NULL)
  filtered.results <- ProcessResultsList(results.list, taxa, TRUE)
  trees <- SummarizeResults(filtered.results, output.format="newick")
  expect_equal(class(trees), "character")
  expect_false(anyNA(trees))
  expect_equal(class(ape::read.tree(text=trees[1])), "phylo")
})

