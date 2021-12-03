test_that("datelife_use workflows work", {
  du <- datelife_use(
    input = "Rhea americana, Struthio camelus, Gallus gallus",
    each = FALSE,
    dating_method = "bladj"
  )
  expect_true("phylo" %in% class(du))
  expect_error(datelife_use_datelifequery())
  # testing that phylo workflows work
  ds <- datelife_search(input = du, summary_format = "citations")
  expect_warning(match_all_calibrations(phy = du, calibrations = NULL))
  du00 <- du0 <- du
  du0$edge.length <- NULL
  expect_warning(extract_calibrations_phylo(input = du0))
  expect_error(extract_calibrations_phylo(input = attributes(du)$datelife_query))
  # test a datelifeQuery input with NA as phy:
  dq <- attributes(du)$datelife_query
  dq$phy <- NA
  expect_warning(datelife_use_datelifequery(input = dq))
  # when input is a datelifeQuery it returns it:
  make_datelife_query(input = attributes(du)$datelife_query)
  # test that pathd8 workflow works:
  uc <- use_calibrations(dating_method = "pathd8", phy = du, calibrations = attributes(du)$datelife_calibrations)
  # test that other bladj workflows work:
  sapply(c("mean", "min", "max"), function(x) {
    use_calibrations_bladj(
      phy = du0,
      calibrations = attributes(du)$datelife_calibrations,
      type = x
    )
  })
  # test that absence of phylogeny or no branch lengths returns warning
  expect_warning(use_calibrations_pathd8(phy = NA))
  expect_error(use_calibrations_bladj(phy = NA))
  expect_warning(use_calibrations_pathd8(phy = du0, calibrations = NULL))
  expect_error(use_calibrations_bladj(phy = du0, calibrations = NULL))
  du00$edge.length <- c(NaN, NaN, NaN, NaN)
  expect_warning(use_calibrations_pathd8(phy = du00))
  # test get calibrations from a vector (calls datelife_search and extract_calibrations_phylo)
  gc <- get_calibrations_vector(input = c("Rhea americana", "Struthio camelus", "Gallus gallus"))
  # testing that you can get a datelife result from a tree:
  dr <- get_datelife_result(input = du)
  is_datelife_result_empty(dr)
  extract_calibrations_dateliferesult(input = dr)
  # test an empty datelifeResult object
  input_process(input = du0)
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
  expect_warning(match_all_calibrations(phy = phylo, calibrations = calibs))
})

test_that("object checks work", {
  expect_warning(match_all_calibrations(phy = NULL))
  expect_error(use_calibrations(phy = NULL))
  expect_error(get_calibrations_datelifequery(input = NULL))
})


test_that("you can load opentree_chronograms",{
  data(opentree_chronograms, package = "datelife")
  expect_equal(ls(opentree_chronograms), c("authors","curators","dois","studies","trees"))
})


# test match_all_calibrations when all(all_nodes < ape::Ntip(phy)) is not TRUE

# test make_bold_otol_tree that does not get a phylo object

# test all function check_conflicting_calibrations in use_calibrations_bladj.R
