# in this file we:
# debug by running the main functions over and over but with different set of taxa
# and register unexpected results from tnrs services to be hopefully corrected in the future

test_that("birds from wikipedia work", {
  # birds_wiki is defined in helper_test_data.R
  xx <- datelife_search(birds_wiki, summary_format = "phylo_median")
  expect_true(inherits(xx, "phylo"))
  xx <- datelife_search(birds_wiki, summary_format = "phylo_sdm")
  expect_true(inherits(xx, "phylo"))
})
