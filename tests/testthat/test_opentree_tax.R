test_that("get_ott_clade works", {
    xx <- get_ott_clade(input = c("random", "Homo"))
    xx <- get_ott_clade(input = c("random", "Homo"), ott_rank = c("family", "order", "class"))
    xx <- get_ott_clade(input = c("perro", "canis", "Homo"), ott_rank = c("family", "order", "class", "genus"))
    xx <- get_ott_clade(input = c("Lamiaceae", "Campanulaceae", "Fabaceae"), ott_rank = "family")
    expect_false(all(is.na(xx$family))) # this should return the same lineages since they are all family level already
    xx <- get_ott_clade(input = c("Lamiaceae", "Campanulaceae", "Salvia"), ott_rank = "family")
    expect_false(all(is.na(xx$family))) # this should return the same lineages since they are all family level already
})
