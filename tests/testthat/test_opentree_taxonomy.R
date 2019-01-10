test_that("get_valid_children works", {
    xx <- get_valid_children("Felis")
    yy <- tnrs_match(rownames(xx$Felis$children))
    expect_false("EXTINCT" %in%  yy$flags) # we need a function to clean taxonomy_taxon_info objects
})

test_that("clean_taxon_info_children works", {
    # Felis ott_id is 563165
    taxon_info <- rotl::taxonomy_taxon_info(563165, include_children = TRUE)
    clean_taxon_info_children(taxon_info)
})
