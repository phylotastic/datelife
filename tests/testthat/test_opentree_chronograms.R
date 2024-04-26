test_that("get_opentree_chronograms function runs", {
  xx <- update_datelife_cache(write = TRUE, max_tree_count = 2) # runs in 2 minutes
  expect_true(inherits(xx$trees, "multiPhylo"))
})

test_that("opentree_chronograms object is ok", {
  data("opentree_chronograms")
  # ls(opentree_chronograms)
  # opentree_chronograms$version
  # opentree_chronograms$trees[[1]]
  expect_true(
    all(sapply(opentree_chronograms$trees, inherits, "phylo"))
    )
})
