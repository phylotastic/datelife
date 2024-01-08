test_that("get_opentree_chronograms function runs", {
  xx <- update_datelife_cache(write = TRUE, max_tree_count = 2) # runs in 2 minutes
  # expect_true(inherits(xx$trees), "multiPhylo")
})
