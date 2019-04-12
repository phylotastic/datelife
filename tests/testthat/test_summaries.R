# testing datelife functions to summarize source trees
test_that("get_taxon_summary works", {
	xx <- get_taxon_summary(threebirds_result, threebirds_query)
	expect_true(all(names(xx) %in% c("matrix", "summary", "absent_taxa")))
  # enhance: it should test if it's giving both matrix and summary
})

test_that("summarize_datelife_result works", {
  # taxa <- "Pan"
  xx <- summarize_datelife_result(threebirds_result, threebirds_query, summary_format = "phylo_all", taxon_summary = "summary")
  xx <- summarize_datelife_result(threebirds_result, threebirds_query, summary_format = "phylo_median", taxon_summary = "summary")
  xx <- summarize_datelife_result(threebirds_result, threebirds_query, summary_format = "phylo_sdm", taxon_summary = "summary")
})

test_that("Summarize as newick_median works correctly", {
  tree <- summarize_datelife_result(datelife_result = threebirds_result, summary_format="newick_median", cache=opentree_chronograms)
  expect_true(inherits(tree, "character"))
  expect_false(anyNA(tree))
  expect_equal(class(ape::read.tree(text=tree)), "phylo")
})

test_that("Summarize as mrca works correctly", {
  mrca.vector <- summarize_datelife_result(datelife_result = threebirds_result, summary_format="mrca", cache=opentree_chronograms)
  expect_equal(class(mrca.vector), "numeric")
  expect_gte(min(mrca.vector),5)
  expect_lte(max(mrca.vector),150)
})


test_that("Summarize as citations works correctly", {
  citation.results <- summarize_datelife_result(datelife_result = threebirds_result, summary_format="cit", cache=opentree_chronograms)
  expect_equal(class(citation.results), "character")
  expect_gte(sum(grepl("Prum", citation.results)),1)
 })

test_that("Summarize as newick_all works correctly", {
  trees <- summarize_datelife_result(datelife_result = threebirds_result, summary_format="newick_all", cache=opentree_chronograms)
  trees <- sapply(datelife_result, patristic_matrix_to_newick)
  expect_equal(class(trees), "character")
  expect_false(anyNA(trees))
  expect_equal(class(ape::read.tree(text=trees[1])), "phylo")
})

test_that("taxon_summary argument from summarize_datelife_result() works", {
  trees <- summarize_datelife_result(datelife_result = threebirds_result, summary_format="newick_all", taxon_summary = "summary")
  expect_gte(length(trees), 3)
  trees2 <- summarize_datelife_result(datelife_result = threebirds_result, summary_format="newick_all", taxon_summary = "matrix")
  expect_gte(length(trees2), 3)
})

test_that("get_dated_otol_induced_subtree works", {
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
  # next one returns no plants
  xx <- get_dated_otol_induced_subtree(input = c("Felis silvestris", "Homo sapiens", "Hamamelidaceae", "Altingiaceae"))
  # we should add an element to pylo object that contains input lineages that are excluded from tree:
  # expect_false(is.null(xx$dropped))
})
