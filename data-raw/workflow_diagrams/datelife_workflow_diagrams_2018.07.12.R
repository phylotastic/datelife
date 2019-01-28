#  July, 12 2018
# datelife congruify workflow
# 

# foodweb plot:
# install.packages("mvbutils")
library(mvbutils)
pdf(file="~/Google Drive/datelife/workflow_diagrams/use_all_calibrations_foodweb1.pdf", height=10, width=25)
foodweb(where = "package:datelife", prune = "use_all_calibrations",
        border = TRUE,
        boxcolor = "aliceblue",
        textcolor = "black", cex = 1.0, lwd=2)
mtext("use_all_calibrations function foodweb", line = 2, font = 2, cex = 1.5)
dev.off()


# run some examples with diagram package:
# install.packages("diagram")
require(diagram)
openplotmat(main = "textbox shapes")
rx <- 0.1
ry <- 0.05
pos <- coordinates(c(1, 1, 1, 1, 1, 1, 1,1 ), mx = -0.2)
textdiamond(mid = pos[1,], radx = rx, rady = ry, lab = LETTERS[1], cex = 2, shadow.col = "lightblue")

# actual diagrams:
pdf(file = "congruify_workflow_use_all_calibrations.pdf", height = 14)
elpos <- coordinates (c(1, rep(2, 11)))
elpos2 <- coordinates (c(1, 2, 2, rep(1, 9)))
from <- c(1, 1, 2:5, 14:19)
to <-   c(2, 3, 4:7, 16:21)
nr <- length(from)
arrpos <- matrix(ncol = 2, nrow = nr)
par(mar = c(1, 1, 1, 1))
openplotmat()
for (i in 1:nr){
  arrpos[i, ] <- straightarrow (to = elpos[to[i], ],
                                from = elpos[from[i], ],
                                lwd = 2, arr.pos = 0.6, arr.length = 0.5)
}
# textellipse(elpos[1,], 0.1, 0.05, lab = c("datelife_search", "function"), box.col = "white",
#               shadow.col = "black", shadow.size = 0.005, cex = 0.8)
textrect (elpos[1,], 0.15, 0.025,lab = c("use_all_calibrations", "function"), box.col = "white",
          shadow.col = "black", shadow.size = 0.005, cex = 0.9)
texthexa(elpos[2,], 0.15, 0.025, lab = c("input =","taxon names"), box.col = "orange",
         shadow.col = "red", shadow.size = 0.005, cex = 0.8)

texthexa (elpos[3,], 0.15, 0.025,lab = c("input =", "a tree"), box.col = "#AFEEEE",
          shadow.col = "darkgreen", shadow.size = 0.005, cex = 0.9)

textrect(elpos[4,], 0.15, 0.025, lab = c("make_bold_otol_tree", "function"), box.col = "orange",
         shadow.col = "red", shadow.size = 0.005, cex = 0.8)
textrect (elpos2[6,], 0.3, 0.015,lab = c("get_all_calibrations function"), box.col = "white",
          shadow.col = "black", shadow.size = 0.005, cex = 0.9)
texthexa (elpos2[10,], 0.3, 0.015, lab = c("congruify with scale = NA", "output is a table of calibrations"), box.col = "green",
          shadow.col = "darkgreen", shadow.size = 0.005, cex = 0.8)
textrect (elpos2[11,], 0.3, 0.015,lab = c("geiger::PATHd8.phylo function"), box.col = "white",
          shadow.col = "black", shadow.size = 0.005, cex = 0.9)
textrect (elpos2[12,], 0.3, 0.015,lab = c("expand calibrations until !is.null(chronogram)"), box.col = "white",
          shadow.col = "black", shadow.size = 0.005, cex = 0.9)
textrect (elpos2[13,], 0.3, 0.025,lab = c("output is ONE chronogram scaled", "with secondary calibrations"), box.col = "white",
          shadow.col = "black", shadow.size = 0.005, cex = 0.9)
dev.off()

pdf(file = "congruify_workflow_datelife_search.pdf")
elpos <- coordinates (c(1, rep(2, 6)))
elpos2 <- coordinates (c(1, 2, 1, 1, 2, 2, 1))
from <- c(2, 2, 3:12) - 1
to <-   c(3, 4, 5:14) - 1 
nr <- length(from)
arrposA <- matrix(ncol = 2, nrow = nr)
par(mar = c(1, 1, 1, 1))
openplotmat()
for (i in 1:nr){
  arrposA[i, ] <- straightarrow (to = elpos[to[i], ],
                                 from = elpos[from[i], ],
                                 lwd = 2, arr.pos = 0.6, arr.length = 0.5)
}
# textellipse(elpos[1,], 0.1, 0.05, lab = c("datelife_search", "function"), box.col = "white",
#               shadow.col = "black", shadow.size = 0.005, cex = 0.8)
textrect (elpos[1,], 0.15, 0.05,lab = c("datelife_search", "function"), box.col = "white",
          shadow.col = "black", shadow.size = 0.005, cex = 0.9)
texthexa(elpos[2,], 0.15, 0.05, lab = c("input =","taxon names"), box.col = "pink",
         shadow.col = "red", shadow.size = 0.005, cex = 0.8)

texthexa (elpos[3,], 0.15, 0.05,lab = c("input =", "a tree"), box.col = "lightblue",
          shadow.col = "darkgreen", shadow.size = 0.005, cex = 0.9)
# textrect (elpos2[4,], 0.3, 0.03,lab = "make_datelife_query function: process names and/or tree", box.col = "white",
#           shadow.col = "black", shadow.size = 0.005, cex = 0.9)
textrect (elpos2[4,], 0.3, 0.03,lab = "get_datelife_results and inside functions", box.col = "white",
          shadow.col = "black", shadow.size = 0.005, cex = 0.9)
texthexa (elpos2[5,], 0.3, 0.05,lab = c("subset chronograms from otol"), box.col = "white",
          shadow.col = "black", shadow.size = 0.005, cex = 0.9)
texthexa (elpos[9,], 0.15, 0.05, lab = c("congruify", "scaling = PATHd8"), box.col = "green",
          shadow.col = "darkgreen", shadow.size = 0.005, cex = 0.8)
texthexa(elpos[10,], 0.15, 0.05, lab = c("output:", "chronograms from", "primary studies ***"), box.col = "pink",
         shadow.col = "red", shadow.size = 0.005, cex = 0.8) # chronograms with tips found in otol only

texthexa (elpos[11,], 0.15, 0.05,lab = c("output:", "chronograms scaled","with PATHd8"), box.col = "lightblue",
          shadow.col = "darkgreen", shadow.size = 0.005, cex = 0.8)
textrect (elpos2[10,], 0.3, 0.05,lab = c("summarize_datelife_result function:", "(a) all chronograms, (b) one summary chronogram,", " (c) a table of ages, (d) citations of primary studies"), box.col = "white",
          shadow.col = "black", shadow.size = 0.005, cex = 0.9)

dev.off()


pdf(file = "congruify_workflow_tree_add_dates.pdf")
elpos <- coordinates (c(1, rep(3, 6)))
elpos2 <- coordinates (c(1, 3, 3, 3, rep(1, 3)))
from <- c(1, 1, 1, 2:10)
to <-   c(2, 3, 4, 5:13)  
nr <- length(from)
arrposA <- matrix(ncol = 2, nrow = nr)
par(mar = c(1, 1, 1, 1))
openplotmat()
for (i in 1:nr){
  arrposA[i, ] <- straightarrow (to = elpos[to[i], ],
                                 from = elpos[from[i], ],
                                 lwd = 2, arr.pos = 0.6, arr.length = 0.5)
}
# textellipse(elpos[1,], 0.1, 0.05, lab = c("datelife_search", "function"), box.col = "white",
#               shadow.col = "black", shadow.size = 0.005, cex = 0.8)
textrect (elpos[1,], 0.15, 0.05,lab = c("***tree_add_dates", "function"), box.col = "white",
          shadow.col = "black", shadow.size = 0.005, cex = 0.9)
texthexa(elpos[2,], 0.15, 0.05, lab = c("adding_criterion =","random"), box.col = "plum",
         shadow.col = "violet", shadow.size = 0.005, cex = 0.8)
texthexa(elpos[3,], 0.15, 0.05, lab = c("adding_criterion =","taxonomy"), box.col = "mistyrose",
         shadow.col = "red", shadow.size = 0.005, cex = 0.8)
texthexa(elpos[4,], 0.15, 0.05, lab = c("adding_criterion =","tree"), box.col = "azure",
         shadow.col = "darkgreen", shadow.size = 0.005, cex = 0.8)
textrect(elpos[7,], 0.15, 0.05, lab = c("make_bold_otol_tree","function"), box.col = "azure",
         shadow.col = "darkgreen", shadow.size = 0.005, cex = 0.8)
texthexa(elpos[10,], 0.15, 0.05, lab = c("congruify", "scaling = NA"), box.col = "green",
         shadow.col = "darkgreen", shadow.size = 0.005, cex = 0.8)
textrect (elpos2[11,], 0.4, 0.05,lab = c("date with: mrbayes, bladj, PATHd8"), box.col = "white",
          shadow.col = "black", shadow.size = 0.005, cex = 0.9)
dev.off()