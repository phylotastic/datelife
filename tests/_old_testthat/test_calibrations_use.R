test_that("use_calibrations errors when necessary", {
  # error when phy is not provided by user:
  expect_error(use_calibrations())
  phy <- "((Pterocnemia_pennata,Rhea_americana),Struthio_camelus);"
  # error when tree is not a phylo object:
  expect_error(use_calibrations(phy = phy))
  # error when calibrations are not provided by user:
  phy <- ape::read.tree(text = phy)
  expect_error(use_calibrations(phy = phy))
})

# test_that("use_calibrations works",{
#
# })
