test_that("plot_phylo_all works", {
    utils::data(felid_gdr_phylo_all)
    plot_phylo_all(felid_gdr_phylo_all$phylo_all)
    phyloch::axisGeo(GTS = NULL, unit = c("period","epoch"),
      col = c("gray80", "white"), gridcol = c("gray80", "white"), cex = 0.5,
      gridty = "twodash")
})
