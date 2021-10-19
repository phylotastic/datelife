test_that("datelife_use", {
  datelife_use(input = "Rhea americana, Struthio camelus, Gallus gallus",
               each = FALSE,
               dating_method = "bladj")
  #
  expect_warning(datelife_use_datelifequery())
})

test_that("tree workflow"), {
  tree <- make_bold_otol_tree(c("Rhea americana",  "Struthio camelus", "Gallus gallus"))
  datelife_search(input = tree)
}
# test_that(),{
#   datelife_use_datelifequery(input = a datelifequery object with phy as NA)
#   use_all_calibration(each = TRUE)
#   use_calibrations(phy = a multiPhylo object, calibrations = NULL)
#   use_calibrations(phy = NULL)
#   use_calibrations(dating_method = "pathd8")
#   extract_calibrations_phylo(input = a multiPhylo object with no branch lengths)
#   extract_calibrations_phylo(input = a multiPhylo object with some trees with no branch lengths)
#   extract_calibrations_phylo(input = a Phylo object with no branch lengths)
#   extract_calibrations_dateliferesult(input = a datelifeResult object)
#   get_calibrations_vector()
#   get_calibrations_datelifequery(input = NULL)
#
# }
