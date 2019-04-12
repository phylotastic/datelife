test_that("clean_taxon_info_children works", {
    # Felis ott_id is 563165
    taxon_info <- rotl::taxonomy_taxon_info(563165, include_children = TRUE)
    # clean_taxon_info_children(taxon_info)
})
test_that("get_valid_children works", {
    xx <- get_valid_children("Felis")
    yy <- tnrs_match(rownames(xx$Felis$children))
    expect_false("EXTINCT" %in%  yy$flags)
})
test_that("get_ott_lineage works", {
    xx <- get_ott_lineage(input = c("Homo"))
    xx <- get_ott_lineage(input = c("random", "Homo")) # fix: this does not work well yet
    xx <- get_ott_lineage(input = c("perro", "canis", "Homo"))
    xx <- get_ott_lineage(input = c("Lamiaceae", "Campanulaceae", "Fabaceae"))
    xx <- get_ott_lineage(input = c("Lamiaceae", "Campanulaceae", "Salvia"))
})
test_that("get_ott_clade works", {
    xx <- get_ott_clade(input = c("random", "Homo"))
    xx <- get_ott_clade(input = c("random", "Homo"), ott_rank = c("family", "order", "class"))
    xx <- get_ott_clade(input = c("perro", "canis", "Homo"), ott_rank = c("family", "order", "class", "genus"))
    xx <- get_ott_clade(input = c("Lamiaceae", "Campanulaceae", "Fabaceae"), ott_rank = "family")
    expect_false(all(is.na(xx$family))) # this should return the same lineages since they are all family level already
    xx <- get_ott_clade(input = c("Lamiaceae", "Campanulaceae", "Salvia"), ott_rank = "family")
    expect_false(all(is.na(xx$family))) # this should return the same lineages since they are all family level already
})
