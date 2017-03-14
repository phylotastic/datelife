test_that("ProcessResultsList", {
utils::data(opentree_chronograms)
taxa <- c("Rhea americana", "Pterocnemia pennata", "Struthio camelus")
results.list <- lapply(opentree_chronograms$trees,GetSubsetArrayDispatch, taxa=taxa, phy=NULL)
expect_gte(length(ProcessResultsList(results.list, taxa, TRUE)), 1)
})

test_that("Summarize as mrca works correctly", {
  utils::data(opentree_chronograms)
  taxa <- c("Rhea americana", "Pterocnemia pennata", "Struthio camelus")
  results.list <- lapply(opentree_chronograms$trees,GetSubsetArrayDispatch, taxa=taxa, phy=NULL)
  filtered.results <- ProcessResultsList(results.list, taxa, TRUE)
  mrca.vector <- SummarizeResults(filtered.results, output.format="mrca", cache=opentree_chronograms)
  expect_equal(class(mrca.vector), "numeric")
  expect_gte(min(mrca.vector),5)
  expect_lte(max(mrca.vector),150)
})

test_that("Summarize as citations works correctly", {
  utils::data(opentree_chronograms)
  taxa <- c("Rhea americana", "Pterocnemia pennata", "Struthio camelus")
  results.list <- lapply(opentree_chronograms$trees,GetSubsetArrayDispatch, taxa=taxa, phy=NULL)
  filtered.results <- ProcessResultsList(results.list, taxa, TRUE)
  citation.results <- SummarizeResults(filtered.results, output.format="cit", cache=opentree_chronograms)
  expect_equal(class(citation.results), "character")
  expect_gte(sum(grepl("Prum", citation.results)),1)
 })

test_that("Summarize as newick.all works correctly", {
  utils::data(opentree_chronograms)
  taxa <- c("Rhea americana", "Pterocnemia pennata", "Struthio camelus")
  results.list <- lapply(opentree_chronograms$trees,GetSubsetArrayDispatch, taxa=taxa, phy=NULL)
  filtered.results <- ProcessResultsList(results.list, taxa, TRUE)
  trees <- SummarizeResults(filtered.results, output.format="newick.all", cache=opentree_chronograms)
  expect_equal(class(trees), "character")
  expect_false(anyNA(trees))
  expect_equal(class(ape::read.tree(text=trees[1])), "phylo")
  expect_equal(class(ape::read.tree(text=trees[length(trees)])), "phylo")
  
})

test_that("Summarize as newick.median works correctly", {
  utils::data(opentree_chronograms)
  taxa <- c("Rhea americana", "Pterocnemia pennata", "Struthio camelus")
  results.list <- lapply(opentree_chronograms$trees,GetSubsetArrayDispatch, taxa=taxa, phy=NULL)
  filtered.results <- ProcessResultsList(results.list, taxa, partial=FALSE)
  tree <- SummarizeResults(filtered.results, output.format="newick.median", cache=opentree_chronograms)
  expect_equal(class(tree), "character")
  expect_false(anyNA(tree))
  expect_equal(class(ape::read.tree(text=tree)), "phylo")
})

test_that("Processing input newick", {
	skip_on_cran()
	skip_on_travis() #b/c no pathd8
  utils::data(opentree_chronograms)
  input.processed <- ProcessInput(ape::write.tree(ape::rcoal(3, tip.label=c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"))), usetnrs=FALSE, approximatematch=TRUE)
  expect_equal(class(input.processed$phy),"phylo")
})

# test_that("Processing complex newick works", {
# 	skip_on_cran()
# 	skip_on_travis()
#   utils::data(opentree_chronograms)
# 	expect_error(EstimateDates("((((((Typha latifolia,(Phragmites australis,(Sporobolus alterniflorus,Sporobolus pumilus)Sporobolus)PACMAD clade)Poales,(((Hydrilla verticillata,Vallisneria americana)Hydrocharitaceae,Potamogeton perfoliatus),Zostera marina,Ruppia maritima)Alismatales),(Lythrum salicaria,Myriophyllum spicatum)),(Ulva,Caulerpa taxifolia))Chloroplastida,((Skeletonema,(Gomphonema,Didymosphenia geminata)Bacillariophyceae)Bacillariophytina,Prorocentrum)SAR),Microcystis)Eukaryota;", output.format="phylo.all"), NA)
# })

test_that("EstimateDates returns phylo.all", {
  skip_on_cran()
  utils::data(opentree_chronograms)
  phylo.results <- EstimateDates(input=c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"), partial=TRUE, usetnrs=FALSE, approximatematch=TRUE, cache=get("opentree_chronograms"), output.format="phylo.all")
  expect_equal(class(phylo.results[[1]]),"phylo")
})

test_that("EstimateDates returns phylo.sdm", {
  skip_on_cran()
  skip_on_travis() #b/c no pathd8
  utils::data(opentree_chronograms)
phylo.results <- EstimateDates(input=c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"), partial=TRUE, usetnrs=FALSE, approximatematch=TRUE, cache=get("opentree_chronograms"), output.format="phylo.sdm")
expect_equal(class(phylo.results),"phylo")
expect_true(!is.null(phylo.results$edge.length))
})


test_that("EstimateDates returns mrca", {
  skip_on_cran()
  utils::data(opentree_chronograms)
  phylo.results <- EstimateDates(input=c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"), partial=TRUE, usetnrs=FALSE, approximatematch=TRUE, cache=get("opentree_chronograms"), output.format="mrca")
expect_equal(class(phylo.results),"numeric")
expect_gte(length(phylo.results), 2)
})

test_that("EstimateDates returns mrca", {
  skip_on_cran()
  utils::data(opentree_chronograms)
  phylo.results <- EstimateDates(input=c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"), partial=TRUE, usetnrs=FALSE, approximatematch=TRUE, cache=get("opentree_chronograms"), output.format="mrca")
expect_equal(class(phylo.results),"numeric")
expect_gte(length(phylo.results), 2)
})

test_that("GetFilteredResults works", {
  skip_on_cran()
  skip_on_travis() #b/c no pathd8
  utils::data(opentree_chronograms)
  filtered.results.in <- GetFilteredResults(input=c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"), partial=TRUE, usetnrs=FALSE, approximatematch=TRUE, cache=get("opentree_chronograms"))
  expect_equal(class(filtered.results.in), "list")
  expect_gte(length(filtered.results.in),4) #as of Nov 4, 2016, had length 8
  expect_equal(class(filtered.results.in[[1]]), "matrix")
})

test_that("Processing input string with spaces", {
	skip_on_cran()
  utils::data(opentree_chronograms)
  input.processed <- ProcessInput(c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"), usetnrs=FALSE, approximatematch=TRUE)
  expect_equal(length(input.processed$cleaned.names),3)
})

test_that("Processing input string with underscores", {
  skip_on_cran()
  utils::data(opentree_chronograms)
  input.processed <- ProcessInput(c("Rhea_americana", "Pterocnemia_pennata", "Struthio_camelus"), usetnrs=FALSE, approximatematch=TRUE)
  expect_equal(length(input.processed$cleaned.names),3)
})

test_that("Making OToL and BOLD tree works", {
  utils::data(opentree_chronograms)
  phy.pars <- GetBoldOToLTree(input=c("Rhea americana",  "Struthio camelus","Gallus gallus", "Pterocnemia pennata"), marker="COI", otol_version="v2", doML=FALSE)
  phy.ml <- GetBoldOToLTree(input=c("Rhea americana",  "Struthio camelus","Gallus gallus", "Pterocnemia pennata"), marker="COI", otol_version="v2", doML=TRUE)
  expect_equal(class(phy.pars), "phylo")
  expect_equal(class(phy.ml), "phylo")
  expect_gte(max(phy.pars$edge.length), 1)
  expect_lte(min(phy.ml$edge.length), 1)
})

test_that("Congruification works", {
	skip_on_cran()
	skip_on_travis() #b/c no pathd8
  utils::data(opentree_chronograms)
  filtered.results <- GetFilteredResults(ape::write.tree(ape::rcoal(3, tip.label=c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"))), partial=TRUE, usetnrs=FALSE, approximatematch=TRUE)
  expect_gte(length(filtered.results), 2)
})

test_that("SDM correctly returns tree", {
  skip_on_cran()
  skip_on_travis() #b/c no pathd8
  utils::data(opentree_chronograms)
  taxa <- c("Rhea americana", "Pterocnemia pennata", "Struthio camelus")
  results.list <- lapply(opentree_chronograms$trees,GetSubsetArrayDispatch, taxa=taxa, phy=NULL)
  filtered.results <- ProcessResultsList(results.list, taxa, TRUE)
  result.tree <- RunSDM(filtered.results)$phy
  expect_equal(class(result.tree), "phylo")
})

test_that("UseAllCalibrations actually works", {
  utils::data(opentree_chronograms)
  results <- UseAllCalibrations()
  expect_true(ape::is.ultrametric(results$phy, tol=0.00001))
  expect_equal(class(results$phy), "phylo")
})

test_that("Crop plant taxa work", {
  utils::data(opentree_chronograms)
  taxa <- c("Zea mays", "Oryza sativa", "Arabidopsis thaliana", "Glycine max", "Medicago sativa", "Solanum lycopersicum")
  results <- EstimateDates(input=taxa, output.format="phylo.all")
  expect_equal(class(results), "list")
  expect_equal(class(results[[1]]), "phylo")
  expect_gte(length(results), 2)
})

test_that("Processing newick input works", {
  processed <- ProcessInput("((Zea mays,Oryza sativa),((Arabidopsis thaliana,(Glycine max,Medicago sativa)),Solanum lycopersicum));")
  expect_equal(class(processed$phy), "phylo")
  expect_equal(ape::Ntip(processed$phy), 6)
  expect_equal(ape::Nnode(processed$phy), 5)
  expect_equal(length(processed$cleaned.names), 6)
})

test_that("Processing newick input works with underscores", {
  processed <- ProcessInput("((Zea_mays,Oryza_sativa),((Arabidopsis_thaliana,(Glycine_max,Medicago_sativa)),Solanum_lycopersicum));")
  expect_equal(class(processed$phy), "phylo")
  expect_equal(ape::Ntip(processed$phy), 6)
  expect_equal(ape::Nnode(processed$phy), 5)
  expect_equal(length(processed$cleaned.names), 6)
})

test_that("Crop plant newick works", {
skip_on_cran()
skip_on_travis() #b/c no pathd8

trees <- EstimateDates(input = "((Zea mays,Oryza sativa),((Arabidopsis thaliana,(Glycine max,Medicago sativa)),Solanum lycopersicum)Pentapetalae);", output.format = "phylo.all", partial = TRUE, usetnrs = FALSE, approximatematch = TRUE, cache = opentree_chronograms, method = "PATHd8")
expect_equal(class(trees[[1]]), "phylo")
})


# test_that("TNRS with approximate match works", {
	 # taxa <- c("Rhea_americana", "Pterocnemia pennato", "Strutho camelus")
 	# input.processed <- ProcessInput(taxa, usetnrs=TRUE, approximatematch=TRUE)
 	# expect_true(all.equal(c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"), input.processed$cleaned.names))
# })


# test_that("TNRS with unmatchable taxa works", {
	# taxa <- c("Rhea_americana", "Pterocnemia pennato", "Oscar the grouch", "Strutho camelus")
 	# input.processed <- ProcessInput(taxa, usetnrs=TRUE, approximatematch=TRUE)
 	# expect_true(all.equal(c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"), input.processed$cleaned.names))
# })
