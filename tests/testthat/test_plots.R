test_that("wrap_string_to_plot works", {
    wrap_string_to_plot(string = "Tree", max_cex = 0.75, whole = FALSE)
    string <- "Cetaceae Species Presence across chronograms in DateLife Data Base"
    # when calling plotSimmap from RStudio the nest line
    # graphics::par("din")[2]-graphics::par("pin")[2]- graphics::par("omi")[1]-graphics::par("mai")[1] - 0.2
    # it gives a negative number! so check this before plotting
})

test_that("plot_phylo_all works", {
    utils::data(felid_gdr_phylo_all)
    plot_phylo_all(felid_gdr_phylo_all$phylo_all)
    # phyloch::axisGeo(GTS = NULL, unit = c("period","epoch"),
    #   col = c("gray80", "white"), gridcol = c("gray80", "white"), cex = 0.5,
    #   gridty = "twodash")
})

# getting an error when phangorn::densitree plotting datelife_result chronograms from the following taxa
test_that("plot_densitree works", {
    taxa <- c("Rhea americana", "Pterocnemia pennata", "Struthio camelus", "Gallus gallus")
    four_birds <- datelife_search(input = taxa, summary_format = "phylo_all")
    plot_densitree(trees = four_birds, include_all = FALSE)
    plot_densitree(trees = four_birds, include_all = TRUE)
    skip("example generating an error from package phangorn")
    taxa <- c("Rhea americana", "Pterocnemia pennata", "Struthio camelus", "Gallus", "Felis")
    birds_and_cats <- datelife_search(input = taxa, summary_format = "phylo_all", get_spp_from_taxon = TRUE)
    plot_densitree(trees = birds_and_cats, include_all = FALSE)
    plot_densitree(trees = birds_and_cats, include_all = TRUE)
  })
