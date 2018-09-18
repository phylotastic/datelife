# test that groves functions are doing what they're supposed to do

test_that("test "){
    utils::data(names_subset2)
    subset2_query <- make_datelife_query(names_subset2)
    subset2_dl <- get_datelife_result(subset2_query)
    grove_taxa <- filter_for_grove(subset2_dl, criterion = "taxa")
    grove_tree <- filter_for_grove(subset2_dl, criterion = "tree")
}
