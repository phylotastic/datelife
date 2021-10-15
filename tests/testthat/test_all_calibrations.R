test_that("getting, mapping and using all calibrations work", {
  # 1) TEST datelife with a vector of names
  # tests make_bold_otol_tree and get_dated_otol_induced_subtree from datelife
  u1 <- datelife_use(input = c("Rhea americana", "Struthio camelus", "Gallus gallus"))
  expect_true(ape::is.ultrametric(u1$phy, option = 2))
  expect_s3_class(u1$phy, "phylo")
  # TODO:
  # u1.1 <- datelife_use(input = c("Felis catus", "Canis canis", "Elephas maximus"))
  # Giving error:
  # Error in convertAlnRows(result$msa, type) : There is an invalid aln file!
    # In addition: Warning message:
    # In scan(file = file, what = what, sep = sep, quote = quote, dec = dec,  :
  # EOF within quoted string

  # TODO
  # these are taking too long to get bold tree I think? Check it out:
  # u1.2 <- datelife_use(input = c("Felis catus", "Homo sapiens", "Elephas maximus"))
  # u1.2 <- datelife_use(input = c("Delphinus delphus", "Homo sapiens", "Elephas maximus"))
  u1.5 <- datelife_use(input = c("Chen caerulescens", "Cygnus columbianus", "Anas acuta"))

  #TEST match_all_calibrations
  expect_equal(match_all_calibrations(phy = NULL, calibrations= NULL), NA)
  # match_all_calibrations(phy = u1.2$phy, calibrations= NULL)
  expect_equal(match_all_calibrations(phy = u1.5$phy, calibrations = u1$calibrations.df), NA)
  expect_equal(match_all_calibrations(phy = u1$phy, calibrations = u1.5$calibrations.df), NA)

  # 2) TEST datelife with a tree with NO branch lengths
  u2 <- datelife_use(input=threebirds_nbl)
  # initial branch lengths are required for pathd8 dating
  # the following will use bladj instead of pathd8 bc phy has no branch lengths:
  c1 <- use_calibrations(phy = threebirds_nbl, calibrations = u2$calibrations, dating_method = "pathd8")
  # test dating_method attribute is "bladj"
  # get all calibrations works with a tree with branch lengths too!
  c2 <- get_all_calibrations(input=threebirds_nbl)

  # 3) TEST datelife with a tree WITH branch lengths
  u3 <- suppressWarnings(datelife_use(input=threebirds_median))
  g3 <- get_all_calibrations(input=threebirds_median)
  g3.2 <- get_all_calibrations(input=threebirds_median, each = TRUE)


  # 4) TEST datelife with a multiPhylo object
  u4.1 <- datelife_use(input=threebirds_all) # tests argument each = FALSE
  u4.2 <- datelife_use(input=threebirds_all, each = TRUE)
  g4 <- get_all_calibrations(input=threebirds_all, each = TRUE)

  threebirds_all0 <- threebirds_all
  threebirds_all0[[1]]$edge.length <- NULL
  get_all_calibrations(input=threebirds_all0)

  threebirds_all_nbl <- c(threebirds_nbl, threebirds_all0[[1]])
  class(threebirds_all_nbl) <- 'multiPhylo'
  get_all_calibrations(input=threebirds_all_nbl)


  # 5) TEST datelife with a datelifeResult object
  # must throw error here:
  expect_error(suppressWarnings(datelife_use(input=threebirds_result)))
  # but not here:
  get_all_calibrations(input= threebirds_result)

  # TODO: 6) TEST datelife with a datelife query object

  # TEST datelife with pathd8

  skip_on_cran()
  use_calibrations(phy = threebirds_median, calibrations = u2$calibrations, dating_method = "pathd8")
  use_calibrations(phy = threebirds_median, calibrations = u2$calibrations, dating_method = "PATHd8")

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
# datelife_use(phy = aves_tree_7000, all_calibrations = samp25_mat17_47)
