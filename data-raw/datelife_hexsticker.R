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

datelife_hexsticker <- function(subplot,
                                p_x = c(0.6, 1.2),
                                p_y = c(1, 0.975),
                                p_size = c(135, 170),
                                s_x = 1,
                                s_y = 0.95,
                                s_width = 0.93,
                                s_height= 0.55,
                                h_fill = "white",
                                filename = "inst/figures/datelife-hex.png",
                                ...) {
  s <- hexSticker::sticker(subplot = subplot,
            package = c("date", "life"),
            p_color = c("black", green),
            p_family = c("Aller_Rg", "gochi"),
            p_x = p_x,
            p_y = p_y,
            p_size = p_size,
            s_x = s_x,
            s_y = s_y,
            s_width = s_width,
            s_height = s_height,
            h_fill = h_fill,
            h_color = green,
            h_size = 2.5,
            filename = filename,
            dpi = 1600,
            white_around_sticker = TRUE
          )
    return(s)
}
