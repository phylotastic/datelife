if (!requireNamespace("hexSticker")) {
  install.packages("hexSticker")
}
# https://www.pngwing.com/en/search?q=phylogenetic+Tree
# https://www.pngwing.com/en/free-png-nndso
# Load text types
if(!requireNamespace("showtext")){
  install.packages(showtext)
}
## Loading Google fonts (http://www.google.com/fonts):
sysfonts::font_add_google("Gochi Hand", "gochi")
## Automatically use showtext to render text for future devices:
showtext::showtext_auto()

if (!requireNamespace("strap")) {
  install.packages("strap")
}

if (!requireNamespace("BiocManager")) {
  install.packages("BiocManager")
}

if (!requireNamespace("ggtree")) {
  BiocManager::install("ggtree")
}

if (!requireNamespace("emojifont")) {
  install.packages("emojifont")
}

library(emojifont)
library(ggimage)

if (!requireNamespace("datelife")) {
  devtools::install_github("phylotastic/datelife")
}
