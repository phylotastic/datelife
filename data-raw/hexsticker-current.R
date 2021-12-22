### `datelife`'s hexsticker latest version (date)

#### The chronogram

The chronogram as a `phylo` object is generated at random with:

```r
set.seed(1525)
phy <- ape::rcoal(n = 15)
```

We modify branch lengths to a max age (including the root) of 530 Mya, a bit after the start of the Cambrian (approximately at 541 Mya):

```
phy$edge.length <- (phy$edge.length * 510) / max(ape::branching.times(phy))
phy$root.edge <- 10
ape::plot.phylo(phy, root.edge = TRUE)
```

The chronogram plot was generated with `strap`:

```{r eval = FALSE}
#devtools::install_github("phylotastic/datelifeplot")
#datelifeplot::plot_phylo(phylo_sdm, title = "", plot_type = "strap")

library(strap) # required to load the geochronostratigraphic scale
phylo_length <- max(ape::branching.times(phy)) + phy$root.edge
# time_depth <- round(phylo_length*1.2, digits = -1)
phy$root.time <- phylo_length
unit = c("Era", "Period", "Epoch")
pdf("inst/figures/phy-strap-rightwards.pdf",
    width = 3,
    height = 1.8,
    bg = "transparent")
strap::geoscalePhylo(tree = phy,
                     direction = "rightwards",
                     x.lim = c(0, phylo_length),
                     cex.tip = 0.7,
                     show.tip.label = FALSE,
                     cex.ts = 0.001, # 0.2
                     width = 3,
                     tick.scale = "no",
                     boxes = unit[length(unit)],
                     quat.rm = TRUE,
                     units = unit,
                     root.edge = TRUE,
                     edge.color = "#a0a0a0") # "#e5e5e5" https://www.color-hex.com/color/e5e5e5
dev.off()
```

#### The hexsticker

The hexsticker is generated with `hexSticker` and a custom function `datelife_hexsticker()` that lives in the `data-raw` folder:

```{r eval = FALSE}
source("data-raw/datelife_hexsticker.R")

# Generate the sticker:
imgurl <- "~/pj_datelife/datelife/inst/figures/phy-strap-rightwards.pdf"
# background <- "#e5e5e5" #gray, default to "white"
datelife_hexsticker(subplot = imgurl,
                    package = " ",
                    p_x = c(0.69, 1.27),
                    p_y = c(1.28, 1.245),
                    p_size = c(135, 160),
                    s_x = 1.026,
                    s_y = 0.96,
                    s_width = 0.975,
                    s_height= 0.25,
                    filename = "man/figures/datelife-hexsticker-latest.png")
```
