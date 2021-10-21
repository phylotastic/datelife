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
  expect_warning(match_all_calibrations(phy = tree, calibrations = NULL))
  tree$edge.length <- NULL
  expect_warning(extract_calibrations_phylo(input = tree))
  # test a datelifeQuery input with NA as phy:
  expect_warning(datelife_use_datelifequery(input = attributes(tree)$query))
  make_datelife_query(input = attributes(tree)$query)
})

test_that("input processing a newick string and multiPhylo workflows work", {
  newick <- "(Gallus_gallus,(Rhea_americana,Struthio_camelus)Palaeognathae)Aves;"
  phylo <- input_process(newick)
  newickBL <- "(Gallus_gallus:165.8333333,(Rhea_americana:128,Struthio_camelus:128)Palaeognathae:37.83333333)Aves;"
  phyloBL <- input_process(newickBL)
  multiphy <- structure(list(phylo, phylo), class = "multiPhylo")
  expect_warning(extract_calibrations_phylo(input = multiphy))
  multiphyBL <- structure(list(phylo, phyloBL), class = "multiPhylo")
  expect_warning(calibs <- extract_calibrations_phylo(input = multiphyBL))
  expect_error(use_calibrations(phy = multiphyBL, calibrations = NULL))
  use_calibrations(phy = multiphyBL, calibrations = calibs)
  use_all_calibrations(phy = multiphyBL, calibrations = calibs, each = TRUE)
  phylo$tip.label <- c("rooster", "nandu", "ostrich")
  expect_warning(match_all_calibrations(phy = phylo, calibrations  = calibs))

})

test_that("object checks work",{
  expect_warning(match_all_calibrations(phy = NULL))
  expect_error(use_calibrations(phy = NULL))
  expect_error(get_calibrations_datelifequery(input = NULL))
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
