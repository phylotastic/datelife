# hexsticker current version
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

install.packages("strap")
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
# Get absolute path of chronogram plot:
imgurl <- system.file("inst/figures/phylo-sdm-emoji0.pdf", package = "datelife")

# https://www.pngwing.com/en/search?q=phylogenetic+Tree
# https://www.pngwing.com/en/free-png-nndso
# Load text types
install.packages(showtext)
## Loading Google fonts (http://www.google.com/fonts):
sysfonts::font_add_google("Gochi Hand", "gochi")
## Automatically use showtext to render text for future devices:
showtext::showtext_auto()

# Generate the sticker:
install.packages("hexSticker")
green <- "#79c843" # datelife green color obtained with https://imagecolorpicker.com/
imgurl <- "~/pj_datelife/datelife/inst/figures/phylo-sdm-strap-upwards-gray.pdf"
background <- "#e5e5e5" #gray "white" # 

s <- hexSticker::sticker(imgurl,
          package = c("date", "life"), 
          p_color = c("black", green),
          p_x = c(1, 1.5),
          p_family = c("Aller_Rg", "gochi"),
          p_y = c(1,0.75),
          p_size = c(135,140), 
          s_x = 0.98, 
          s_y = 1, 
          s_width = 0.96, 
          s_height= 0.4,
          h_fill = background, 
          h_color = green, 
          h_size = 2.5,
          filename = "inst/figures/datelife-hex-upwards-gray.png",
          dpi = 1600)
          
```