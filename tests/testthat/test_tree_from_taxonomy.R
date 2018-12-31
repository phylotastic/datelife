test_that("classification_paths_from_taxonomy works", {
    taxa <- c("Homo sapiens", "twilight sparkle", "Equus quagga", "Archaeopteryx")
    results <- classification_paths_from_taxonomy(taxa)
    expect_true(inherits(results$resolved, "data.frame"))
    expect_true(length(results$unresolved)>0)
})

test_that("tree_from_taxonomy works", {
    taxa <- c("Homo sapiens", "twilight sparkle", "Equus quagga", "Archaeopteryx")
    results <- tree_from_taxonomy(taxa)
    expect_equal(inherits(results$phy, "phylo"))
})
