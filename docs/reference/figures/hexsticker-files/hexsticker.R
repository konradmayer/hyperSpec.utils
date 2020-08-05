library(hyperSpec)
library(hexSticker)
library(cowplot)

dat <- read.txt.Witec.Graph('man/figures/hexsticker-files/hexbg (Header).txt', encoding = 'latin1', type = 'map')
spc <- read.jdx("man/figures/hexsticker-files/coniferylalkohol.jdx", encoding = 'latin1')

spcplt <- ~{plot(rev(spc[ , , 1000~1800]$spc[1, ]),
                  type = 'l', bty = 'n', axes = F, ylab = "", xlab = "")}

pltbg <- plotmap_viridis(dat[, , 1600], spc ~ -x * -y)

plt <- ggdraw() +
  draw_plot(pltbg) +
  draw_plot(spcplt, width = 0.8, height = 0.8,
            hjust = 0.09, vjust = -0.15)

hex_sticker <- sticker(plt,  s_width = 3.5, s_height = 3.8,
        s_x = 1.3,  s_y = 0.9, white_around_sticker = TRUE,
        package = "hyperSpec.utils", p_size = 18, p_x = 1.1, p_y = 1.3,
        h_color = "#a7da34",
        filename = "man/figures/hexsticker.png",
        url = "konradmayer.github.io/hyperSpec.utils", u_size = 5.3)

saveRDS(hex_sticker,'man/figures/hexsticker-files/hexsticker.rds')
