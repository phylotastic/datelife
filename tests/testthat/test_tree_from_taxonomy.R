test_that("tree_from_taxonomy works with Phyllostomidae (reticulations)", {
    taxa_dq <- make_datelife_query("Phyllostomidae", get_spp_from_taxon = TRUE)
    taxa <- unname(taxa_dq$cleaned_names)
    taxtree_ncbi <- tree_from_taxonomy(taxa, source = "NCBI")
    # taxtree$phy$tip.label[is.na(match(taxtree$phy$tip.label, taxa))]
    taxtree_col <- tree_from_taxonomy(taxa, source = "Catalogue of Life")
    # this runs but gives a bad tree
    # Artibeus aztecus gives trouble in Catalogue of Life because it is assigned
    # a subspecies, so ranks do not coincide anymore, ending up with nodes that
    # are parents of themselves. Compare to Dermanura azteca wich has a coherent taxonomy.
    # we still need to develop the test, this is just the example:
    classification_results <- classification_paths_from_taxonomy(taxa=c("Dermanura azteca", "Artibeus aztecus"), sources="Catalogue of Life")
    paths <- classification_results$resolved$classification_path
})

test_that("classification_paths_from_taxonomy works", {
    taxa <- c("Homo sapiens", "twilight sparkle", "Equus quagga", "Archaeopteryx")
    results <- classification_paths_from_taxonomy(taxa)
    expect_true(inherits(results$resolved, "data.frame"))
    expect_true(length(results$unresolved)>0)
    # the following happens when using cached opentree_chronograms from load(data-raw/opentree_chronograms_oct2018.rda)
    # classification_paths_from_taxonomy(opentree_chronograms$trees[[50]], sources = "Open Tree of Life Reference Taxonomy")
    # gives: Error: Request-URI Too Long (HTTP 414)
    # not sure yet why it happens
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

test_that("get_fossil_range works", {
    dates <- get_fossil_range("Tyrannosaurus rax") # yep, with mispellings
    expect_true(all(dates$min_ma > 64))

    dates <- get_fossil_range("Tyrannosaurus rex", recent=TRUE) # I'm not allowed to say how I know this, but T. rex is still alive.
    expect_true(min(dates$min_ma) == 0)

})

test_that("summarize_fossil_range works", {
    dates <- summarize_fossil_range("Tyrannosaurus rax") # yep, with mispellings
    expect_true(dates$min_ma > 64)
    expect_true(nrow(dates)==1)
    expect_true(rownames(dates)=="Tyrannosaurus rax") # we get the name we put in (which may match our tree), not the GNR'ed name
})
