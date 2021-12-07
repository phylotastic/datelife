######### HEXSTICKER DEVELOPMENT

BiocManager::install("ggtree")
install.packages("emojifont")
library(emojifont)
library(ggimage)

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
############# STRAP: 
########## UPWARDS phlyogeny:
tree <- phylo_sdm
phylo_length <- max(ape::branching.times(tree))
time_depth <- round(phylo_length*1.2, digits = -1)
tree$root.time <- phylo_length

pdf("inst/figures/phylo-sdm-strap-upwards.pdf",
    width = 3, 
    height = 1.8, 
    units = "in"
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
                     units = unit, arotate = 45)
dev.off()

########## RIGHTWARDS phlyogeny:
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
                     cex.ts = 0.9,
                     width = 4,
                     tick.scale = "no",
                     boxes = unit[length(unit)],
                     quat.rm = TRUE,
                     units = unit)
dev.off()

#################### CREATING THE HEXSTICKER
imgurl <- "~/pj_datelife/datelife/inst/figures/phylo-sdm-emoji0.pdf"
imgurl <- "~/pj_datelife/datelife/inst/figures/darwin-i-think.png"
imgurl <- "~/pj_datelife/datelife/inst/figures/phylo-sdm-strap-upwards.pdf"
imgurl <- "~/pj_datelife/datelife/inst/figures/phylo-sdm-strap-rightwards.pdf"
s <- hexSticker::sticker(imgurl,
          package = c("date", "life"), 
          p_color = c("black", green),
          p_x = c(0.8, 1.23),
          p_family = c("Aller_Rg", "gochi"),
          p_y = c(0.5,0.46),
          p_size = 107, 
          s_x = 1, 
          s_y = 1.1, 
          s_width = 0.75, 
          s_height= 0.5,
          h_fill = "white", 
          h_color = green, 
          h_size = 2.5,
          filename = "inst/figures/datelife-hex.png",
          dpi = 1600)

