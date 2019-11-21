
#' Lineage through time plots of  all chronograms in multiphylo object and its summaries.
#'
#' @param phy_summ A multiphylo object with chronograms from a summary.
#' @param phy_summ_type A character value.
#' @param phy_summ_col A character value.
#' @param max_tips A numeric value
#' @param length_arrowhead A numeric value
#' @param ... parameters passed to ape::ltt.lines
ltt_summary <- function(phy_summ, phy_summ_type = NULL, phy_summ_col = NULL, max_tips,
  length_arrowhead = 0.075, ...){
    if(!inherits(phy_summ_type, "character")){
      phy_summ_type <- "Summary"
    }
    if(inherits(phy_summ, "phylo")){
      phy_summ <- list(phy_summ)
    }
    foo <- function(phy, color_here, labels_here, length_arrowhead, max_tips){
        ape::ltt.lines(phy = phy, col = paste0(color_here, "90"), lty = 1, lwd = 2, ...)
        # points(x = -max(ape::branching.times(tax_phycluster[[i]])), y = 2, pch = 25, col = paste0(col_here, "60"), lwd = 0.75)
        x0 <- x1 <- -max(ape::node.depth.edgelength(phy))
        graphics::arrows(x0, y0 = 2.5+max_tips*0.1, x1, y1 = 2.5, length = length_arrowhead,
            col = paste0(color_here, "99"), lwd = 2.5, lty = 1)
        graphics::text(x = x0, y = 2.5+max_tips*0.14, labels = labels_here, srt = 45,
            adj = 0, cex = 0.85, col = color_here, font = 2)
    }
    lab <- paste(gsub("_tree", " ", names(phy_summ)), ifelse("median" %in% tolower(phy_summ_type),
      tolower(phy_summ_type), phy_summ_type))
    for (i in seq(phy_summ)){
      foo(phy = phy_summ[[i]], color_here = phy_summ_col, labels_here = lab[i], length_arrowhead, max_tips)
    }
}
#' Lineage through time plots of all chronograms in multiphylo object and its summaries.
#'
#' @param taxon Character vector indicating the name of the taxon or lineage that the chronograms in phy belong to.
#' @param phy A phylo or multiphylo object with chronograms (trees with branch lengths proportional to geologic time), ideally.
#' @param phy_sdm A multiphylo object with chronograms from a SDM summary.
#' @param phy_median A multiphylo object with chronograms from a median summary.
#' @param tax_datedotol A chronogram to compare other chronograms to.
#' @param file_name A character string giving the name of the pdf file.
#' @param file_dir A character string giving the path to write the file to.
#' @param height Height of the plot
#' @param width Width of the plot
#' @inheritParams ape::plot.phylo
#' @inheritDotParams ape::ltt.plot
#' @param add_legend Boolean
#' @param add_title Boolean
#' @param col_sdm Color of the SDM tree
#' @param col_median Color of the median tree
#' @export
# modified from make_lttplot_summchrono2 function in datelife_examples
plot_ltt_summary <- function(taxon, phy, phy_sdm, phy_median,
        file_name = NULL, file_dir = NULL, height = 3.5, width = 7,
        add_legend = TRUE, add_title = FALSE, col_sdm = "#00AFBB", col_median = "#CC79A7",
        tax_datedotol = NULL, ...){

    if(!inherits(taxon, "character")){
      taxon <- "Some species"
    }
    if(!inherits(file_dir, "character")){
      file_dir <- "~//datelife//data-raw//"
    }
    if(!inherits(file_name, "character")){
      file_name <- paste0(gsub(" ", "_", taxon), "_LTTplot_summary.pdf")
    }
    file_out <- paste0(file_dir, file_name)
    phy_mrca <- sapply(phy, function(x) max(ape::branching.times(x)))
    leg <- "source chronograms"
    leg_color <- "#77889980"
    trees <- phy
    if(inherits(tax_datedotol, "phylo")){
        tax_datedotol <- ape::collapse.singles(tax_datedotol)
        tax_datedotol <- phytools::force.ultrametric(tax_datedotol)
        trees <- c(trees, tax_datedotol)
        leg <- c(leg, "dated OToL tree")
        leg_color <- c(leg_color, "#80808080")
    }
    max_ages <- sapply(trees, function(x) max(ape::branching.times(x)))
    xlim0 <- round(max(max_ages)+5, digits = -1)
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # general variables for source chronogram plotting:
    col_sample <- paste0("#778899", sample(20:90, length(unique(names(phy))))) #color is lightslategrey
    treesall <- c(trees, phy_sdm, phy_median)
    max_tips <- max(sapply(treesall, function(x) max(ape::Ntip(x))))
    length_arrowhead <- 0.075
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # start the plot:
    grDevices::pdf(file = file_name, height = height, width = width)
    graphics::par(mai = c(1.02, 0.82, 0.2, 0.2))
    ltt_phyloall(phy, trees = treesall, max_tips, max_ages, xlim0, taxon, phy_mrca,
      col_sample, length_arrowhead, lwd_phyloall = 1.5, ...)
    ltt_summary(phy_summ = phy_median, phy_summ_type = "Median",
        phy_summ_col = col_median, max_tips, length_arrowhead) # default color pinkish
    ltt_summary(phy_summ = phy_sdm, phy_summ_type = "SDM",
        phy_summ_col = col_sdm, max_tips, length_arrowhead) # default color teal
    if(add_legend){
        leg <- paste(taxon, c("source chronograms", "summary chronograms"))
        graphics::legend(x = -xlim0, y = max_tips*1.1, legend = leg, cex = 0.75, pch = 19,
            bty = "n", xpd = NA, col = c("#77889980", "#00AFBB80"), inset = -1)
    }
    if(add_title){
      graphics::mtext(text = "SDM and median summary chronograms", side = 3, cex = 1, font = 2, line = -1.5)
    }
    grDevices::dev.off()
}
