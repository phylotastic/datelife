# `datelife`'s hexsticker

We used https://imagecolorpicker.com/ to find hex color codes from an image

We set up package requirements in the file `data-raw/hexsicker_env.R`.
`sourc`ing it checks that required packages are installed, and installs them if not:

```
source("data-raw/hexsicker_env.R")
```

## An `emojifont` tree with `ggtree` and `ggimage`

```
phylo_sdm <- datelife::datelife_search(input = c("Delphinus_delphis", "Gallus gallus", "elephas Maximus", "felis_catus", "homo-sapiens"),
                                       use_tnrs = TRUE,
                                       summary_format = "phylo_sdm")
class(phylo_sdm) <- "phylo"

phylo_sdm$tip.label_original <- phylo_sdm$tip.label
phylo_sdm$tip.label <- c("elephant",
                         "woman_cartwheeling",
                         "dolphin",
                         "cat",
                         "rooster")

p <- ggtree::ggtree(phylo_sdm, size=0.5, color='black') +
  ggtree::geom_tiplab(parse = "emoji", size=15, vjust=0.01) +
  ggtree::theme_tree(bgcolor = "#00640010") +
  coord_cartesian(clip = 'off') +
  ggtree::xlim_tree(c(0, 350))

pdf("inst/figures/phylo-sdm-emoji0.pdf",
    width = 10,
    height = 5,
    bg = "transparent")
print(p)
dev.off()
```

## A `strap` tree

### UPWARDS phlyogeny

```
tree <- phylo_sdm
phylo_length <- max(ape::branching.times(tree))
time_depth <- round(phylo_length*1.2, digits = -1)
tree$root.time <- phylo_length
library(strap)

pdf("inst/figures/phylo-sdm-strap-upwards-gray.pdf",
    width = 3,
    height = 1.8,
    bg = "transparent")
unit = c("Era", "Period", "Epoch")
strap::geoscalePhylo(tree = tree,
                     direction = "upwards",
                     x.lim = c(0, phylo_length),
                     cex.tip = 0.7,
                     show.tip.label = FALSE,
                     cex.ts = 0.5,
                     cex.age = 0.7,
                     width = 4,
                     tick.scale = "no",
                     boxes = unit[length(unit)],
                     quat.rm = TRUE,
                     units = unit, arotate = 45,
                     edge.color = "#708090")
dev.off()
```

### RIGHTWARDS phlyogeny

```
unit = c("Era", "Period")

pdf("inst/figures/phylo-sdm-strap-rightwards-black.pdf",
    width = 3,
    height = 2,
    bg = "transparent")
strap::geoscalePhylo(tree =  ape::rotate(tree, node=6),
                     direction = "rightwards",
                     x.lim = c(0, phylo_length),
                     cex.tip = 0.7,
                     show.tip.label = FALSE,
                     cex.ts = 0.9,
                     width = 4,
                     tick.scale = "no",
                     boxes = unit[length(unit)],
                     quat.rm = TRUE,
                     units = unit,
                     edge.color = "black") #"#708090" # gray
dev.off()
```

## CREATING THE HEXSTICKER

```
imgurl <- "~/pj_datelife/datelife/inst/figures/phylo-sdm-emoji0.pdf"
imgurl <- "~/pj_datelife/datelife/inst/figures/darwin-i-think.png"
imgurl <- "~/pj_datelife/datelife/inst/figures/phylo-sdm-strap-upwards.pdf"


imgurl <- "~/pj_datelife/datelife/inst/figures/phylo-sdm-strap-rightwards-black.pdf"
green <- "#79c843" # datelife green color obtained with https://imagecolorpicker.com/
background <- "white" # "#e5e5e5" #gray
s <- hexSticker::sticker(imgurl,
          package = c("date", "life"),
          p_color = c("black", green),
          p_x = c(0.6, 1.2),
          p_y = c(1,0.975),
          p_family = c("Aller_Rg", "gochi"),
          p_size = c(135, 170),
          s_x = 1,
          s_y = 0.95,
          s_width = 0.93,
          s_height= 0.55,
          h_fill = background,
          h_color = green,
          h_size = 2.5,
          filename = "inst/figures/datelife-hex-rightward-black.png",
          dpi = 1600)
```

## A hexsticker with an empty background, using an empty plot

### 1) Choose a background and rim color

```
green <- "#79c843" # datelife green color obtained with https://imagecolorpicker.com/
background <- "white" # "#e5e5e5" #gray
```

### 2) Create an empty image

```
pdf("inst/figures/empty.pdf",
    width = 2,
    height = 2,
    bg = "transparent")
plot.new() # creates an empty plot
dev.off()
```

### 3) Create the hexSticker

```
imgurl <- "inst/figures/empty.pdf"
s <- hexSticker::sticker(subplot = imgurl,
          package = "",
          h_fill = background,
          h_color = green,
          h_size = 2.5,
          filename = "inst/figures/datelife-hex-empty.png",
          dpi = 1600)
```
