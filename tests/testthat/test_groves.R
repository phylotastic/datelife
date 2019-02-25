# test that groves functions are doing what they're supposed to do

test_that("test groves", {
    utils::data(names_subset2)
    subset2_query <- make_datelife_query(names_subset2)
    subset2_dl <- get_datelife_result(subset2_query)
    grove_taxa <- filter_for_grove(subset2_dl, criterion = "taxa")
    grove_tree <- filter_for_grove(subset2_dl, criterion = "tree")
})

# test that best grove obtained for median method in datelife_summary function
# is the best one that can be used for sdm (maybe groves with smaller n work for
# sdm and not for median and visceversa)


test_that("median and sdm give a tree when source trees have different degrees of overlapping names", {
    utils::data(names_subset2)
    spp_query <- make_datelife_query(names_subset2)
    spp_dl_result <- get_datelife_result(spp_query)
    # length(spp_dl_result)  #24
    xx <- get_best_grove(spp_dl_result)
    
    xxphylo <- sdm_matrix_to_phylo(xx)
    # sdm_matrix_to_phylo(spp_dl_result)
    x1 <- summarize_datelife_result(datelife_query = spp_query,
      datelife_result = spp_dl_result, summary_format = "phylo_median")
    expect_true(inherits(summarize_datelife_result(datelife_query = spp_query,
      datelife_result = spp_dl_result, summary_format = "phylo_median"), "phylo"))
    expect_true(inherits(summarize_datelife_result(datelife_query = spp_query,
      datelife_result = spp_dl_result, summary_format = "phylo_sdm"), "phylo"))
    expect_equal(class(summarize_datelife_result(datelife_query = spp_query,
      datelife_result = spp_dl_result, summary_format = "newick_median")), "character")
    expect_equal(class(summarize_datelife_result(datelife_query = spp_query,
      datelife_result = spp_dl_result, summary_format = "newick_sdm")), "character")
})
