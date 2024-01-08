test_that("extract_calibrations_phylo works", {
  class(felid_gdr_phylo_all$phylo_all) <- "multiPhylo"
  xx <- extract_calibrations_phylo(input = felid_gdr_phylo_all$phylo_all)
  expect_true(inherits(xx, "data.frame"))
})
