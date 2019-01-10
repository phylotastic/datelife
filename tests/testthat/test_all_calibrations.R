test_that("get_otol_synthetic_tree works", {
  otol_tree <- get_otol_synthetic_tree(input = c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"))
  expect_s3_class(otol_tree, "phylo") # output is phylo
  expect_gte(length(otol_tree), 4)  # it has no branch lengths
  # otol_tree <- get_otol_synthetic_tree(input = c("Struthio camelus"))
  # should it return all names always?
  # it should return the same number of names matched by tnrs_match
  # enhance: add such a test
  child <- get_ott_children("Felis")
  input <- child$Felis
  names(input)
  get_otol_synthetic_tree(input = child$Felis)

})
