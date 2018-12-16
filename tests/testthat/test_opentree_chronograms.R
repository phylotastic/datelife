test_that("get_otol_chronograms works", {
	# skip_on_cran()
  #	skip_on_travis() #b/c no pathd8
  # skip_on_os("linux") #b/c no pathd8 on travis linux

  utils::data(opentree_chronograms)
  # xx <- get_otol_chronograms(verbose=TRUE, max_tree_count = 10)

})

test_that("is_good_chronogram works", {
  utils::data(felid_gdr_phylo_all)
  t1 <- felid_gdr_phylo_all$phylo_all[[1]]
  t1$tip.label[1] <-  "*tip_#1_not_mapped_to_OTT._Original_label_-_Hagensia_havilandi"
  t1$tip.label[2] <-  "*tip_#1_not_mapped_to_OTT._Original_label_-_Felis_silvestris"
  is_good_chronogram(t1)
  t1$tip.label[3] <-  gsub("_", " ", t1$tip.label[1])
  is_good_chronogram(t1)
  t1$tip.label[4] <-  "*tip_#1_not_mapped_to_OTT._Original_label_-_Hagensia_havilandi"
})

test_that("clean_ott_chronogram works as expected", {

})

test_that("opentree_chronograms object is ok", {
  utils::data(opentree_chronograms)
  # write(paste(names(opentree_chronograms), collapse = '", "'), file = "data-raw/names.txt")
  # test that all expected elements are in opentree_chronograms:
  expect_true(all(c("trees", "authors", "curators", "studies", "dois") %in% names(opentree_chronograms)))
  # test that all opentree_chronograms elements have the same length:
  expect_true(all(sapply(opentree_chronograms, length) == length(opentree_chronograms$trees)))
  # opentree_chronograms$trees
})
