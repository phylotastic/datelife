source("data-raw/hextsicker_env.R")

datelife_hexsticker <- function(subplot,
                                package = c("date", "life"),
                                p_x = c(0.6, 1.2),
                                p_y = c(1, 0.975),
                                p_size = c(135, 170),
                                s_x = 1,
                                s_y = 0.95,
                                s_width = 0.93,
                                s_height= 0.55,
                                filename = "inst/figures/datelife-hex.png",
                                ...) {
#
  green <- "#79c843" # datelife green color obtained with https://imagecolorpicker.com/

  s <- hexSticker::sticker(subplot = subplot,
            package = package,
            p_color = c("black", green),
            p_family = c("Aller_Rg", "gochi"),
            p_x = p_x,
            p_y = p_y,
            p_size = p_size,
            s_x = s_x,
            s_y = s_y,
            s_width = s_width,
            s_height = s_height,
            h_fill = "white",
            h_color = green,
            h_size = 0.5,
            filename = filename,
            dpi = 1600,
            white_around_sticker = TRUE
          )
    return(s)
}
