test_that("classification_paths_from_taxonomy works", {
    taxa <- c("Homo sapiens", "twilight sparkle", "Equus quagga", "Archaeopteryx")
    results <- classification_paths_from_taxonomy(taxa)
    expect_true(inherits(results$resolved, "data.frame"))
    expect_true(length(results$unresolved)>0)
})

test_that("tree_from_taxonomy works", {
    taxa <- c("Homo sapiens", "twilight sparkle", "Equus quagga", "Archaeopteryx")
    results <- tree_from_taxonomy(taxa)
    expect_true(inherits(results$phy, "phylo"))
})

test_that("tree_from_taxonomy works with weird inputs and PBDB", {
    taxa <- c("Homo sapiens", "twilight sparkle", "Equus quagga", "Archaeopteryx", "Marchantiophyta", "Polypodiopsida")
    classifications <- classification_paths_from_taxonomy(taxa, sources="The Paleobiology Database")
    expect_true(inherits(classifications$resolved, "data.frame"))
    results <- tree_from_taxonomy(taxa, sources="The Paleobiology Database")
    expect_true(inherits(results$phy, "phylo"))


    taxa <- c("Gorilla", "Panthera", "Tyto", "Dromaius", "Aedes", "Solenopsis", "Caretta", "Crocodylus", "Brassica", "Solanum", "Zea", "Prunus", "Rosa", "Climacograptus")
    classifications <- classification_paths_from_taxonomy(taxa, sources="The Paleobiology Database")
    expect_true(inherits(classifications$resolved, "data.frame"))
    results <- tree_from_taxonomy(taxa, sources="The Paleobiology Database")
    expect_true(inherits(results$phy, "phylo"))
})
