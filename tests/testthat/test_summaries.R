# testing datelife functions to summarize source trees
testthat("get_taxon_summary works", {
  taxa <- c("Rhea americana", "Struthio camelus", "Gallus gallus")
  taxa <- "Pan"
	datelife_query <- make_datelife_query(taxa, get_spp_from_taxon = TRUE)
	datelife_result <- get_datelife_result(taxa)
	xx <- get_taxon_summary(datelife_query, datelife_result)
  # enhance: it should test if it's giving both matrix and summary
})

testthat("summarize_datelife_result works", {
  taxa <- "Pan"
  taxa <- c("Rhea americana", "Struthio camelus", "Gallus gallus")
	datelife_query <- make_datelife_query(taxa, get_spp_from_taxon = TRUE)
	datelife_result <- get_datelife_result(taxa)
  xx <- summarize_datelife_result(datelife_query, datelife_result, summary_format = "phylo_all", taxon_summary = "summary")
  xx <- summarize_datelife_result(datelife_query, datelife_result, summary_format = "phylo_median", taxon_summary = "summary")
  xx <- summarize_datelife_result(datelife_query, datelife_result, summary_format = "phylo_sdm", taxon_summary = "summary")

})

test_that("Summarize as newick_median works correctly", {
  utils::data(opentree_chronograms)
  taxa <- c("Rhea americana", "Pterocnemia pennata", "Struthio camelus")
  results_list <- lapply(opentree_chronograms$trees,get_subset_array_dispatch, taxa=taxa, phy=NULL)
  datelife_result <- results_list_process(results_list, taxa, partial=FALSE)
  tree <- summarize_datelife_result(datelife_result = datelife_result, summary_format="newick_median", cache=opentree_chronograms)
  expect_equal(class(tree), "character")
  expect_false(anyNA(tree))
  expect_equal(class(ape::read.tree(text=tree)), "phylo")
})
test_that("Summarize as mrca works correctly", {
  utils::data(opentree_chronograms)
  taxa <- c("Rhea americana", "Pterocnemia pennata", "Struthio camelus")
  results_list <- lapply(opentree_chronograms$trees,get_subset_array_dispatch, taxa=taxa, phy=NULL)
  datelife_result <- results_list_process(results_list, taxa, TRUE)
  mrca.vector <- summarize_datelife_result(datelife_result = datelife_result, summary_format="mrca", cache=opentree_chronograms)
  expect_equal(class(mrca.vector), "numeric")
  expect_gte(min(mrca.vector),5)
  expect_lte(max(mrca.vector),150)
})

test_that("Summarize as citations works correctly", {
  utils::data(opentree_chronograms)
  taxa <- c("Rhea americana", "Pterocnemia pennata", "Struthio camelus")
  results_list <- lapply(opentree_chronograms$trees,get_subset_array_dispatch, taxa=taxa, phy=NULL)
  datelife_result <- results_list_process(results_list, taxa, TRUE)
  citation.results <- summarize_datelife_result(datelife_result = datelife_result, summary_format="cit", cache=opentree_chronograms)
  expect_equal(class(citation.results), "character")
  expect_gte(sum(grepl("Prum", citation.results)),1)
 })

test_that("Summarize as newick_all works correctly", {
  utils::data(opentree_chronograms)
  taxa <- c("Rhea americana", "Pterocnemia pennata", "Struthio camelus")
  results_list <- lapply(opentree_chronograms$trees,get_subset_array_dispatch, taxa=taxa, phy=NULL)
  datelife_result <- results_list_process(results_list, taxa, TRUE)
  trees <- summarize_datelife_result(datelife_result = datelife_result, summary_format="newick_all", cache=opentree_chronograms)
  expect_equal(class(trees), "character")
  expect_false(anyNA(trees))
  expect_equal(class(ape::read.tree(text=trees[1])), "phylo")
})

test_that("taxon_summary argument from summarize_datelife_result() works", {
  utils::data(opentree_chronograms)
  taxa <- c("Rhea americana", "Pterocnemia pennata", "Struthio camelus")
  results_list <- lapply(opentree_chronograms$trees,get_subset_array_dispatch, taxa=taxa, phy=NULL)
  datelife_result <- results_list_process(results_list, taxa, TRUE)
  trees <- summarize_datelife_result(datelife_result = datelife_result, summary_format="newick_all", cache=opentree_chronograms, taxon_summary = "summary")
  # str(trees)
  # trees$taxon_distribution
  # trees$absent_taxa
  expect_gte(length(trees), 3)
  trees2 <- summarize_datelife_result(datelife_result = datelife_result, summary_format="newick_all", cache=opentree_chronograms, taxon_summary = "matrix")
  expect_gte(length(trees2), 3)
  # str(trees2)
  # trees2$taxon_distribution
  # trees2$absent_taxa
})
test_that("SDM correctly returns tree", {
  skip_on_cran()
  skip_on_os("linux") #b/c no pathd8 on travis linux
  utils::data(opentree_chronograms)
  taxa <- c("Rhea americana", "Pterocnemia pennata", "Struthio camelus")
  results_list <- lapply(opentree_chronograms$trees, get_subset_array_dispatch, taxa=taxa, phy=NULL)
  datelife_result <- results_list_process(results_list, taxa, TRUE)
  result.tree <- datelife_result_sdm(datelife_result)$phy
  expect_true(inherits(result.tree, "phylo"))
})

test_that("get_dated_otol_induced_subtree works", {
    utils::data(felid_sdm)
  xx <- get_dated_otol_induced_subtree()
  xx <- get_dated_otol_induced_subtree(input = felid_sdm$phy)
  xx <- get_dated_otol_induced_subtree(ott_id = c(563163, 770315)) # cat and human ott_ids
  # expect an NA result on the two following:
  xx <- get_dated_otol_induced_subtree(ott_id = c("563163", "770315", "mrcaott99")) # cat and human ott_ids
  xx <- get_dated_otol_induced_subtree(ott_id = c("563163", "mrcaott770315", "mrcaott99")) # cat and human ott_ids
  # "Hamamelidaceae", "Altingiaceae", "Zamiaceae", "Rutaceae", "Saxifragaceae", "Asparagaceae", "Cycadaceae", "Smilacaceae", "Boraginaceae"
  # with respective ott ids 737324, 853767, 614459, 43847, 1035588, 17704, 99242, 978709, 147029
  # are dropped from tree. I'm sure there are more.
  # next one returns NA:
  xx <- get_dated_otol_induced_subtree(input = c("Felis", "Canis", "Hamamelidaceae", "Altingiaceae"))
  # next one returns a tree of cats and wolves, no plants
  xx <- get_dated_otol_induced_subtree(input = c("Felis silvestris", "Homo sapiens", "Hamamelidaceae", "Altingiaceae"))
  # we should add an element to pylo object that contains input lineages that are excluded from tree:
  # expect_false(is.null(xx$dropped))
})
