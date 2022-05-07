test_that("datelife queries work", {
  make_datelife_query(input = "anatidae")
  # get species of penguins
  xx <- make_datelife_query(input = "eudyptes", get_spp_from_taxon = TRUE)
}
