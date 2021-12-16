# hexsticker latest version
#

The chronogram was generated with:

```{r eval = FALSE}
phylo_sdm <- datelife::datelife_search(input = c("Delphinus_delphis", "Gallus gallus", "elephas Maximus", "felis_catus", "homo-sapiens"),
                                       use_tnrs = TRUE,
                                       summary_format = "phylo_sdm")
class(phylo_sdm) <- "phylo"


```

The chronogram plot was generated with:

```{r eval = FALSE}
#devtools::install_github("phylotastic/datelifeplot")
#datelifeplot::plot_phylo(phylo_sdm, title = "", plot_type = "strap")

source("data-raw/datelife-hexsticker.R")
tree <- phylo_sdm
phylo_length <- max(ape::branching.times(tree))
time_depth <- round(phylo_length*1.2, digits = -1)
tree$root.time <- phylo_length
unit = c("Era", "Period")
pdf("inst/figures/phylo-sdm-strap-rightwards.pdf",
    width = 3,
    height = 2,
    bg = "transparent")
strap::geoscalePhylo(tree = tree,
                     direction = "rightwards",
                     x.lim = c(0, phylo_length),
                     cex.tip = 0.7,
                     show.tip.label = FALSE,
                     cex.ts = 0.7,
                     width = 4,
                     tick.scale = "no",
                     boxes = unit[length(unit)],
                     quat.rm = TRUE,
                     units = unit)
dev.off()
```

The hexsticker was generated with:

```{r eval = FALSE}
source("data-raw/datelife_hexsticker.R")

# Generate the sticker:
imgurl <- "~/pj_datelife/datelife/inst/figures/phylo-sdm-strap-upwards-gray.pdf"
background <- "#e5e5e5" #gray, default to "white"
datelife_hexsticker(subplot = imgurl,
                    p_x = c(1, 1.5),
                    p_y = c(1, 0.75),
                    p_size = c(135, 140),
                    s_x = 0.98,
                    s_y = 1,
                    s_width = 0.96,
                    s_height= 0.4,
                    h_fill = background,
                    filename = "inst/figures/datelife-hex-upwards-gray-test.png")
```
