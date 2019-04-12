
test_that("use_all_calibrations actually works", {
  skip_on_cran()
  skip_on_travis() #b/c no pathd8 
  # skip_on_os("linux") #b/c no pathd8 on travis linux
  results <- suppressWarnings(use_all_calibrations())
  # expect_true(ape::is.ultrametric(results$phy, tol=0.0000))
  expect_true(ape::is.ultrametric(results$phy, option = 2))
  expect_s3_class(results$phy, "phylo")
  # skip("data in url is not yet available")
  # # enhance:
  # # the following loads a file from an url:
  # load(url("https://github.com/LunaSare/datelife_benchmark/tree/master/data/0_global/
  # aves_targets/aves_tree_7000.rda"))
  # load(url("https://github.com/LunaSare/datelife_benchmark/tree/master/data/0_global/
  # aves_mat_samples/samp25_mat17_47.rda"))
  # # the following gave a tree with NaN as edge.length, why?
  # use_all_calibrations(phy = aves_tree_7000, all_calibrations = samp25_mat17_47)

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


test_that("get_all_calibrations works", {
  # this function does not use pathd8
  xx <- get_all_calibrations()
  expect_s3_class(xx, "data.frame") # should be a data.frame
})
