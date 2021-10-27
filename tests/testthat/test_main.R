test_that("datelife_use workflows work", {
  du <- datelife_use(input = "Rhea americana, Struthio camelus, Gallus gallus",
               each = FALSE,
               dating_method = "bladj")
  class(du)
  attributes(du)
  expect_warning(datelife_use_datelifequery())
  # testing that phylo workflows work
  tree <- du$phy
  datelife_search(input = tree, summary_format = "citations")
  expect_warning(match_all_calibrations(phy = tree, calibrations = NULL))
  tree$edge.length <- NULL
  expect_warning(extract_calibrations_phylo(input = tree))
  # test a datelifeQuery input with NA as phy:
  expect_warning(datelife_use_datelifequery(input = attributes(tree)$query))
  # input is a datelifeQuery:
  make_datelife_query(input = attributes(du)$datelifeQuery)
  # test that pathd8 workflow works:
  use_calibrations(dating_method = "pathd8", phy = du$phy, calibrations = du$calibrations.df)
  # test get calibrations from a vector (calls datelife_search and extract_calibrations_phylo)
  gc <- get_calibrations_vector(input = c("Rhea americana", "Struthio camelus", "Gallus gallus"))
  # testing that you can get a datelife result from a tree work:
  dr <- get_datelife_result(input = du$phy)
  is_datelife_result_empty(dr)
  extract_calibrations_dateliferesult(input = dr)
  # test an empty datelifeResult object
  
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



# test match_all_calibrations when all(all_nodes < ape::Ntip(phy)) is not TRUE

# test make_bold_otol_tree that does not get a phylo object

