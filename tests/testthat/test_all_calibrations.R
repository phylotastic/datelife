test_that(" getting all calibrations to work", {
  # TEST make_bold_otol_tree and get_dated_otol_induced_subtree
  # from use_all_calibrations:
  r1 <- use_all_calibrations()
  expect_true(ape::is.ultrametric(r1$phy, option = 2))
  expect_s3_class(r1$phy, "phylo")

  # TEST with a tree with NO branch lengths
  r2 <- suppressWarnings(use_all_calibrations(input=threebirds_nbl))
  # initial branch lengths are required for pathd8 dating
  # the following will use bladj instead of pathd8 bc phy has no branch lengths:
  c2 <- use_calibrations(phy = threebirds_nbl, calibrations = r2$calibrations, dating_method = "pathd8")
  # test dating_method attribute is "bladj"
  # get all calibrations should not work with a tree with no branch lengths:
  expect_error(get_all_calibrations(input=threebirds_nbl))

  # TEST a tree WITH branch lengths
  r3 <- suppressWarnings(use_all_calibrations(input=threebirds_median))
  get_all_calibrations(input=threebirds_median)

  #test with a vector of Names
  # use_all_calibrations(input = c("Rhea americana", "Pterocnemia pennata",
                                 # "Struthio camelus"))

  # TEST with a datelife query object

  # TEST with a datelifeResult object
  # must throw error here:
  expect_error(suppressWarnings(use_all_calibrations(input=threebirds_result)))
  # but not here:
  get_all_calibrations(input= threebirds_result)

  # TEST with a multiPhylo object
  xx <- use_all_calibrations(input=threebirds_all)
  get_all_calibrations(input=threebirds_all, each = TRUE)

  threebirds_all0 <- threebirds_all
  threebirds_all0[[1]]$edge.length <- NULL
  get_all_calibrations(input=threebirds_all0)

  threebirds_all_nbl <- c(threebirds_nbl, threebirds_all0[[1]])
  class(threebirds_all_nbl) <- 'multiPhylo'
  get_all_calibrations(input=threebirds_all_nbl)

  #TEST match_all_calibrations
  match_all_calibrations(phy = NULL, calibrations= NULL)
  match_all_calibrations(phy = xx$phy, calibrations= NULL)


  skip_on_cran()
  use_calibrations(phy = threebirds_median, calibrations = r2$calibrations, dating_method = "pathd8")
  use_calibrations(phy = threebirds_median, calibrations = r2$calibrations, dating_method = "PATHd8")

})



test_that("use_calibrations fails when necessary",{
  # when phy is missing
  expect_error(use_calibrations())
  # when calibrations are missing
  expect_error(use_calibrations(phy = threebirds_nbl))
})


test_that("get_otol_synthetic_tree works", {
  otol_tree <- get_otol_synthetic_tree(input = c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"))
  expect_s3_class(otol_tree, "phylo") # output is phylo
  expect_gte(length(otol_tree), 4)  # it has no branch lengths
  # otol_tree <- get_otol_synthetic_tree(input = c("Struthio camelus"))
  # should it return all names always?
  # it should return the same number of names matched by tnrs_match
  # enhance: add such a test
  # child <- get_ott_children("Felis")
  # get_otol_synthetic_tree(input = child$Felis)
})

# TODO
# test the following data set:
# load a file from url
# load(url("https://github.com/LunaSare/datelife_benchmark/tree/master/data/0_global/
# aves_targets/aves_tree_7000.rda"))
# load(url("https://github.com/LunaSare/datelife_benchmark/tree/master/data/0_global/
# aves_mat_samples/samp25_mat17_47.rda"))
# data in url is not yet available :(
# the following gives a tree with NaN as edge.length, why?
# use_all_calibrations(phy = aves_tree_7000, all_calibrations = samp25_mat17_47)
