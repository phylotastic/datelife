# plots all chronograms from a phyloall datelife summary, using a sample of colors from grDevices::rainbow()
# can plot dated otol tree too
#' Easy visualization of lineage through time plots of all chronograms in a multiphylo object.
#'
#' @param taxon Character vector indicating the name of the taxon or lineage that the chronograms in phy belong to.
#' @param phy A phylo or multiphylo object with chronograms (trees with branch lengths proportional to geologic time), ideally.
#' @param ltt_colors A character vector indicating the colors to be used for plotting ltt
#' @param tax_datedotol A chronogram to compare trees in phy to.
#' @param file_name A character string giving the name of the pdf file.
#' @param file_dir A character string giving the path to write the file to.
#' @param height Height of the plot output
#' @param width width of the plot output
#' @inheritParams ape::plot.phylo
#' @inheritDotParams ape::ltt.plot -phy
#' @param add_legend Boolean
#' @param add_title Boolean
#' @param title_text Character vector
#' @param study_number_cex A numeric value
#' @param lwd_arrows A numeric value
#' @export
plot_ltt_phyloall <- function(taxon = NULL, phy, ltt_colors = NULL, tax_datedotol = NULL,
    file_name = NULL, file_dir = NULL, height = 3.5, width = 7, add_legend = FALSE,
    add_title = FALSE, title_text = NULL,  study_number_cex = 0.75, lwd_arrows = 2, ...){

  if(!inherits(taxon, "character")){
    taxon <- "Some species"
  }
  if(inherits(phy, "phylo")){
    phy <- list(phy)
    class(phy) <- "multiPhylo"
    names(phy) <- taxon
  }
  if(!inherits(file_dir, "character")){
    file_dir <- getwd()
  }
  if(!inherits(file_name, "character")){
    file_name <- "_LTTplot_phyloall.pdf"
  }
  file_out <- paste0(file_dir, "//", paste0(gsub(" ", "_", taxon), file_name))
  message(file_out)
  phy_mrca <- sapply(phy, function(x) max(ape::branching.times(x)))
  leg <- "source chronograms"
  leg_col <- "#8B008B" # "darkmagenta"
  trees <- phy
  if(inherits(tax_datedotol, "phylo")){
      # ape::is.ultrametric(tax_datedotol)
      # ape::is.binary(tax_datedotol)
      tax_datedotol <- ape::collapse.singles(tax_datedotol)
      tax_datedotol <- phytools::force.ultrametric(tax_datedotol)
      trees <- c(trees, tax_datedotol)
      leg_col <- c(leg_col, "#808080") #gray
      leg <- c(leg, "dated OToL tree")

  }
  class(trees) <- "multiPhylo"
  # class(phy)
  # ape::is.ultrametric(phy)
  max_ages <- sapply(trees, function(x) max(ape::branching.times(x)))
  xlim0 <- round(max(max_ages)+5, digits = -1)
  max_tipsall <- sapply(trees, function(x) max(ape::Ntip(x)))
  max_tips <- max(max_tipsall)
  # col_phyloall <- "#cce5ff" # blue
  y1 <- -max_tips*0.015
  y0 <- -max_tips*0.075
  length_arrowhead <- 0.075
  nn <- unique(names(phy))[order(unique(names(phy)))] # get ordered names
  # get a sample of colors the size of the number of studies (one color for each study):
  if(!inherits(ltt_colors, "character")){
    # col_sample <- sample(gray.colors(n = length(nn)), length(nn))
    col_sample <- sample(grDevices::rainbow(n = length(nn)), length(nn))
    write(paste0('"', paste(col_sample, collapse = '", "'), '"'),
      file = paste0(file_dir, "//", taxon, "_lttplot_phyloall_colors.txt"))
  } else {
    col_sample <- ltt_colors
  }
  # match the color to each study:
  col_phyloall_sample <- col_sample[match(names(phy), nn)]
  study_number <- seq(length(nn))[match(names(phy), nn)]
  ss <- which(table(study_number)>1)
  for(ii in ss){ # case when a study has multiple chronograms and we need to adjust x position of number
      tt <- which(ii==study_number)
      dd <- abs(diff(phy_mrca[tt]))
      eq <- which(dd < 0.02*xlim0)
      if(length(eq) == 0) next
      for(j in eq){ # uses the mean age for those chronograms that are closer by less than 0.5 myrs
          max_ages[tt[c(j, j+1)]] <- mean(phy_mrca[ii==study_number][c(j, j+1)])
      }
  }
  y_numbers <- rep(-max_tips*0.14, length(max_ages))
  cond1 <- duplicated(round(max_ages)) & !duplicated(study_number)
  y_numbers[cond1] <- -max_tips*0.23

  grDevices::pdf(file = file_out, height = height, width = width)
  graphics::par(mai = c(1.02, 0.82, 0.2, 0.2))
  ape::ltt.plot(trees[[which.max(max_tipsall)]], xlim = c(-xlim0, 0),
        ylim = c(-max_tips*0.30, max_tips),
        col = paste0("#ffffff", "80"), ylab = paste(taxon, "Species N"),
        xlab = "", ...)
  graphics::mtext("Time (MYA)", side = 1, cex = 1, font = 1, line = 2)
  cond2 <- (!duplicated(study_number) | !duplicated(round(max_ages)))
  for (i in order(phy_mrca, decreasing = TRUE)){ # plot the oldest chronogrm first, it looks better in graph
    col_phyloall <- col_phyloall_sample[i]
    ape::ltt.lines(phy = phy[[i]], col = paste0(col_phyloall), lwd = 1.5)
    x0 <- x1 <- -phy_mrca[i]
    graphics::arrows(x0, y0, x1, y1, length = length_arrowhead, col = paste0(col_phyloall), lwd = lwd_arrows)
    graphics::text(x = -max_ages[i], y = y_numbers[i], labels = ifelse(cond2[i], study_number[i], ""),
        font = 4, col = col_phyloall, cex = study_number_cex)
  }
  if(add_title){
    # text(labels = paste(taxon, "source chronograms"), x = -xlim0, y = max_tips*0.925, cex = 0.75, adj = 0, font = 1)
    if(!inherits(title_text, "character")){
      title_text <- "source chronograms"
    }
    graphics::mtext(text = paste(taxon, title_text), side = 3, cex = 1, font = 2, line = -1.5)
  }
  if(add_legend){
      leg <- paste(taxon, leg)
      graphics::legend(x = "topleft", #round(-max_age, digits = -1), y = round(max_tips, digits = -2),
        legend = leg, col = leg_col, cex = 0.5, pch = 19, bty = "n")
  }
  grDevices::dev.off()
}

# enhance: use the following plot inside other plots requiring ltt phylo all plotting
ltt_phyloall <- function(phy, trees = NULL, max_tips, max_ages, xlim0, taxon, phy_mrca,
  col_sample, length_arrowhead = 0.075, lwd_phyloall = 1.5,  study_number_cex = 0.75, ...){

  if(!inherits(trees, "multiPhylo")){
    trees <- phy
  }
  nn <- unique(names(phy))[order(unique(names(phy)))] # get ordered names
  col_phyloall_sample <- col_sample[match(names(phy), nn)]
  study_number <- seq(length(nn))[match(names(phy), nn)]
  ss <- which(table(study_number)>1)
  for(ii in ss){ # case when a study has multiple chronograms and we need to adjust x position of number
      tt <- which(ii==study_number)
      dd <- abs(diff(phy_mrca[tt]))
      eq <- which(dd < 0.02*xlim0)
      if(length(eq) == 0) next
      for(j in eq){ # uses the mean age for those chronograms that are closer by less than 0.5 myrs
          phy_mrca[tt[c(j, j+1)]] <- mean(phy_mrca[ii==study_number][c(j, j+1)])
      }
  }
  y_numbers <- rep(-max_tips*0.14, length(max_ages))
  cond1 <- duplicated(round(max_ages)) & !duplicated(study_number)
  y_numbers[cond1] <- -max_tips*0.23
  ape::ltt.plot(trees[[which.max(max_tips)]], xlim = c(-xlim0, 0),
        ylim = c(-max_tips*0.30, max_tips),
        col = paste0("#ffffff", "10"), ylab = paste(taxon, "Species N"),
        xlab = "", ...) # we need to plot it white because argument plot = FALSE is not working with ltt.plot
  graphics::mtext("Time (MYA)", side = 1, cex = 1, font = 1, line = 2)
  cond2 <- (!duplicated(study_number) | !duplicated(round(max_ages)))
  for (i in order(phy_mrca, decreasing = TRUE)){
    col_phyloall <- col_phyloall_sample[i]
    ape::ltt.lines(phy = phy[[i]], col = paste0(col_phyloall), lwd = lwd_phyloall)
    x0 <- x1 <- -phy_mrca[i]
    graphics::arrows(x0, y0 = -max_tips*0.075, x1, y1 = 0, length = length_arrowhead,
      col = paste0(col_phyloall), lwd = 2)
    graphics::text(x = -max_ages[i], y = y_numbers[i], labels = ifelse(cond2[i],
      study_number[i], ""), font = 4, col = col_phyloall, cex = study_number_cex)
  }
}
