# test that groves functions are doing what they're supposed to do

test_that("test groves", {
  grove_taxa <- filter_for_grove(subset2_result, criterion = "taxa")
  grove_tree <- filter_for_grove(subset2_result, criterion = "tree")
})

# test that best grove obtained for median method in datelife_summary function
# is the best one that can be used for sdm (maybe groves with smaller n work for
# sdm and not for median and visceversa)


test_that("median and sdm give a tree when source trees have different degrees of overlapping names", {
  x1 <- summarize_datelife_result(
    datelife_query = subset2_query,
    datelife_result = subset2_result, summary_format = "phylo_median"
  )
  expect_true(inherits(summarize_datelife_result(
    datelife_query = subset2_query,
    datelife_result = subset2_result, summary_format = "phylo_median"
  ), "phylo"))
  expect_true(inherits(summarize_datelife_result(
    datelife_query = subset2_query,
    datelife_result = subset2_result, summary_format = "phylo_sdm"
  ), "phylo"))
  expect_equal(class(summarize_datelife_result(
    datelife_query = subset2_query,
    datelife_result = subset2_result, summary_format = "newick_median"
  )), "character")
  expect_equal(class(summarize_datelife_result(
    datelife_query = subset2_query,
    datelife_result = subset2_result, summary_format = "newick_sdm"
  )), "character")
})
