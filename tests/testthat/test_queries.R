test_that("datelife queries work", {
  # does not retrieve all species:
  xx <- make_datelife_query(input = "anatidae")
  expect_true(inherits(xx, "datelifeQuery"))
  expect_true(xx$cleaned_names == "Anatidae")
  expect_true(is.numeric(xx$ott_ids))
  expect_true(is.na(xx$phy))
  # get species of penguins
  xx <- make_datelife_query(input = "eudyptes", get_spp_from_taxon = TRUE)
  expect_true(inherits(xx, "datelifeQuery"))
  expect_true(is.character(xx$cleaned_names))
  expect_true(length(xx$cleaned_names) > 1)
  expect_true(is.numeric(xx$ott_ids))
  expect_true(length(xx$ott_ids) > 1)
  expect_true(is.na(xx$phy))
})
