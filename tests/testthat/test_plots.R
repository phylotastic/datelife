test_that("wrap_string_to_plot works", {
    wrap_string_to_plot(string = "Tree", max_cex = 0.75, whole = FALSE)
    string2 <- "Cetaceae Species Presence across chronograms in DateLife Data Base"
    wrap_string_to_plot(string = string2, max_cex = 0.75, whole = FALSE)
    # when calling plotSimmap from RStudio the next line
    # graphics::par("din")[2]-graphics::par("pin")[2]- graphics::par("omi")[1]-graphics::par("mai")[1] - 0.2
    # it gives a negative number! so check this before plotting
})

test_that("plot_phylo_all works", {
    skip_on_cran()
    skip_on_travis()
    skip("plotting phylo_all")
    plot_phylo_all(felid_gdr_phylo_all$phylo_all)
})

# getting an error when phangorn::densitree plotting datelife_result chronograms from the following taxa
test_that("plot_densitree works", {
    skip_on_cran()
    skip_on_travis()
    skip("plotting densitree")
    taxa <- c("Rhea americana", "Pterocnemia pennata", "Struthio camelus", "Gallus gallus")
    four_birds <- datelife_search(input = taxa, summary_format = "phylo_all")
    plot_densitree(trees = four_birds, include_all = FALSE)
    plot_densitree(trees = four_birds, include_all = TRUE)
    taxa <- c("Rhea americana", "Pterocnemia pennata", "Struthio camelus", "Gallus", "Felis")
    birds_and_cats <- datelife_search(input = taxa, summary_format = "phylo_all", get_spp_from_taxon = TRUE)
    plot_densitree(trees = birds_and_cats, include_all = FALSE)
    plot_densitree(trees = birds_and_cats, include_all = TRUE)
  })

# test_that("plot_ltt_phyloall works", {
#     plot_ltt_phyloall(taxon = "Three birds", phy = threebirds_all, ltt_colors = NULL, tax_datedotol = NULL,
#         file_name = NULL, file_dir = file.path(getwd(), "//data-raw"), height = 3.5, width = 7, add_legend = FALSE, add_title = FALSE)
# })
