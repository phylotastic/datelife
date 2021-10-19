test_that("datelife_use works", {
  datelife_use(input = "Rhea americana, Struthio camelus, Gallus gallus",
               each = FALSE,
               dating_method = "bladj")
  #
  expect_warning(datelife_use_datelifequery())
})

test_that("phylo workflows work", {
  # make a tee with branch lengths:
  tree <- make_bold_otol_tree(c("Rhea americana",  "Struthio camelus", "Gallus gallus"))
  datelife_search(input = tree)
  expect_na(match_all_calibrations(phy = tree, calibrations = NULL))
  tree$edge.length <- NULL
  expect_warning(extract_calibrations_phylo(input = tree))
  # test a datelifeQuery input with NA as phy:
  expect_warning(datelife_use_datelifequery(input = attributes(tree)$query))
})

# test_that("multiPhylo workflows work", {
#   extract_calibrations_phylo(input = a multiPhylo object with no branch lengths)
#   extract_calibrations_phylo(input = a multiPhylo object with some trees with no branch lengths)
#   use_all_calibration(each = TRUE)
#   use_calibrations(phy = a multiPhylo object, calibrations = NULL)
#   match_all_calibration(taxon names in phy are not in calibratiosn data frame)
#
# })

test_that("object checks work"),{
  match_all_calibrations(phy = NULL)
  use_calibrations(phy = NULL)
  get_calibrations_datelifequery(input = NULL)
})
#
# test_that("pathd8 workflow works"),{
#   use_calibrations(dating_method = "pathd8")
# })
#
# test_that("datelifeResult workflow works"),{
#   extract_calibrations_dateliferesult(input = a datelifeResult object)
#   get_calibrations_vector()
# })
