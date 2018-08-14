# # a variation of plot.phylo. It allows changing root edge color
# I need to fix this befor euncommenting:
# plot_chronogram.phylo: no visible global function definition for ‘par’
# plot_chronogram.phylo : .nodeHeight: no visible binding for global
#   variable ‘node_height’
# plot_chronogram.phylo : .nodeDepth: no visible binding for global
#   variable ‘node_depth’
# plot_chronogram.phylo : .nodeDepthEdgelength: no visible binding for
#   global variable ‘node_depth_edgelength’
# plot_chronogram.phylo: no visible global function definition for
#   ‘is.ultrametric’
# ... 29 lines ...
#   ‘.PlotPhyloEnv’
# Undefined global functions or variables:
#   .PlotPhyloEnv circular.plot cladogram.plot is.ultrametric node_depth
#   node_depth_edgelength node_height node_height_clado par
#   phylogram.plot plot.default polar2rect rect2polar reorder segments
#   strheight strwidth text unrooted.xy
# Consider adding
#   importFrom("graphics", "par", "plot.default", "segments", "strheight",
#              "strwidth", "text")
#   importFrom("stats", "reorder")
# to your NAMESPACE file.
# #
#
# plot_chronogram.phylo <- function (x, type = "phylogram", use.edge.length = TRUE, node.pos = NULL,
#       show.tip.label = TRUE, show.node.label = FALSE, edge.color = "black",
#       edge.width = 1, edge.lty = 1,
#       root.edge.color = "black", root.edge.width = 1, root.edge.lty = 1,
#       font = 3, cex = par("cex"),
#       adj = NULL, srt = 0, no.margin = FALSE, root.edge = FALSE,
#       label.offset = 0, underscore = FALSE, x.lim = NULL, y.lim = NULL,
#       direction = "rightwards", lab4ut = NULL, tip.color = "black",
#       plot = TRUE, rotate.tree = 0, open.angle = 0, node.depth = 1,
#       align.tip.label = FALSE, ...){
#
#       Ntip <- length(x$tip.label)
#       if (Ntip < 2) {
#           warning("found less than 2 tips in the tree")
#           return(NULL)
#       }
#       .nodeHeight <- function(edge, Nedge, yy) .C(node_height,
#           as.integer(edge[, 1]), as.integer(edge[, 2]), as.integer(Nedge),
#           as.double(yy))[[4]]
#       .nodeDepth <- function(Ntip, Nnode, edge, Nedge, node.depth) .C(node_depth,
#           as.integer(Ntip), as.integer(edge[, 1]), as.integer(edge[,
#               2]), as.integer(Nedge), double(Ntip + Nnode), as.integer(node.depth))[[5]]
#       .nodeDepthEdgelength <- function(Ntip, Nnode, edge, Nedge,
#           edge.length) .C(node_depth_edgelength, as.integer(edge[,
#           1]), as.integer(edge[, 2]), as.integer(Nedge), as.double(edge.length),
#           double(Ntip + Nnode))[[5]]
#       Nedge <- dim(x$edge)[1]
#       Nnode <- x$Nnode
#       if (any(x$edge < 1) || any(x$edge > Ntip + Nnode))
#           stop("tree badly conformed; cannot plot. Check the edge matrix.")
#       ROOT <- Ntip + 1
#       type <- match.arg(type, c("phylogram", "cladogram", "fan",
#           "unrooted", "radial"))
#       direction <- match.arg(direction, c("rightwards", "leftwards",
#           "upwards", "downwards"))
#       if (is.null(x$edge.length)) {
#           use.edge.length <- FALSE
#       }
#       else {
#           if (use.edge.length && type != "radial") {
#               tmp <- sum(is.na(x$edge.length))
#               if (tmp) {
#                   warning(paste(tmp, "branch length(s) NA(s): branch lengths ignored in the plot"))
#                   use.edge.length <- FALSE
#               }
#           }
#       }
#       if (is.numeric(align.tip.label)) {
#           align.tip.label.lty <- align.tip.label
#           align.tip.label <- TRUE
#       }
#       else {
#           if (align.tip.label)
#               align.tip.label.lty <- 3
#       }
#       if (align.tip.label) {
#           if (type %in% c("unrooted", "radial") || !use.edge.length ||
#               is.ultrametric(x))
#               align.tip.label <- FALSE
#       }
#       if (type %in% c("unrooted", "radial") || !use.edge.length ||
#           is.null(x$root.edge) || !x$root.edge)
#           root.edge <- FALSE
#       phyloORclado <- type %in% c("phylogram", "cladogram")
#       horizontal <- direction %in% c("rightwards", "leftwards")
#       xe <- x$edge
#       if (phyloORclado) {
#           phyOrder <- attr(x, "order")
#           if (is.null(phyOrder) || phyOrder != "cladewise") {
#               x <- reorder(x)
#               if (!identical(x$edge, xe)) {
#                   ereorder <- match(x$edge[, 2], xe[, 2])
#                   if (length(edge.color) > 1) {
#                     edge.color <- rep(edge.color, length.out = Nedge)
#                     edge.color <- edge.color[ereorder]
#                   }
#                   if (length(edge.width) > 1) {
#                     edge.width <- rep(edge.width, length.out = Nedge)
#                     edge.width <- edge.width[ereorder]
#                   }
#                   if (length(edge.lty) > 1) {
#                     edge.lty <- rep(edge.lty, length.out = Nedge)
#                     edge.lty <- edge.lty[ereorder]
#                   }
#               }
#           }
#           yy <- numeric(Ntip + Nnode)
#           TIPS <- x$edge[x$edge[, 2] <= Ntip, 2]
#           yy[TIPS] <- 1:Ntip
#       }
#       z <- reorder(x, order = "postorder")
#       if (phyloORclado) {
#           if (is.null(node.pos))
#               node.pos <- if (type == "cladogram" && !use.edge.length)
#                   2
#               else 1
#           if (node.pos == 1)
#               yy <- .nodeHeight(z$edge, Nedge, yy)
#           else {
#               ans <- .C(node_height_clado, as.integer(Ntip), as.integer(z$edge[,
#                   1]), as.integer(z$edge[, 2]), as.integer(Nedge),
#                   double(Ntip + Nnode), as.double(yy))
#               xx <- ans[[5]] - 1
#               yy <- ans[[6]]
#           }
#           if (!use.edge.length) {
#               if (node.pos != 2)
#                   xx <- .nodeDepth(Ntip, Nnode, z$edge, Nedge,
#                     node.depth) - 1
#               xx <- max(xx) - xx
#           }
#           else {
#               xx <- .nodeDepthEdgelength(Ntip, Nnode, z$edge, Nedge,
#                   z$edge.length)
#           }
#       }
#       else {
#           twopi <- 2 * pi
#           rotate.tree <- twopi * rotate.tree/360
#           if (type != "unrooted") {
#               TIPS <- x$edge[which(x$edge[, 2] <= Ntip), 2]
#               xx <- seq(0, twopi * (1 - 1/Ntip) - twopi * open.angle/360,
#                   length.out = Ntip)
#               theta <- double(Ntip)
#               theta[TIPS] <- xx
#               theta <- c(theta, numeric(Nnode))
#           }
#           switch(type, fan = {
#               theta <- .nodeHeight(z$edge, Nedge, theta)
#               if (use.edge.length) {
#                   r <- .nodeDepthEdgelength(Ntip, Nnode, z$edge,
#                     Nedge, z$edge.length)
#               } else {
#                   r <- .nodeDepth(Ntip, Nnode, z$edge, Nedge, node.depth)
#                   r <- 1/r
#               }
#               theta <- theta + rotate.tree
#               if (root.edge) r <- r + x$root.edge
#               xx <- r * cos(theta)
#               yy <- r * sin(theta)
#           }, unrooted = {
#               nb.sp <- .nodeDepth(Ntip, Nnode, z$edge, Nedge, node.depth)
#               XY <- if (use.edge.length) unrooted.xy(Ntip, Nnode,
#                   z$edge, z$edge.length, nb.sp, rotate.tree) else unrooted.xy(Ntip,
#                   Nnode, z$edge, rep(1, Nedge), nb.sp, rotate.tree)
#               xx <- XY$M[, 1] - min(XY$M[, 1])
#               yy <- XY$M[, 2] - min(XY$M[, 2])
#           }, radial = {
#               r <- .nodeDepth(Ntip, Nnode, z$edge, Nedge, node.depth)
#               r[r == 1] <- 0
#               r <- 1 - r/Ntip
#               theta <- .nodeHeight(z$edge, Nedge, theta) + rotate.tree
#               xx <- r * cos(theta)
#               yy <- r * sin(theta)
#           })
#       }
#       if (phyloORclado) {
#           if (!horizontal) {
#               tmp <- yy
#               yy <- xx
#               xx <- tmp - min(tmp) + 1
#           }
#           if (root.edge) {
#               if (direction == "rightwards")
#                   xx <- xx + x$root.edge
#               if (direction == "upwards")
#                   yy <- yy + x$root.edge
#           }
#       }
#       if (no.margin)
#           par(mai = rep(0, 4))
#       if (show.tip.label)
#           nchar.tip.label <- nchar(x$tip.label)
#       max.yy <- max(yy)
#       getLimit <- function(x, lab, sin, cex) {
#           s <- strwidth(lab, "inches", cex = cex)
#           if (any(s > sin))
#               return(1.5 * max(x))
#           Limit <- 0
#           while (any(x > Limit)) {
#               i <- which.max(x)
#               alp <- x[i]/(sin - s[i])
#               Limit <- x[i] + alp * s[i]
#               x <- x + alp * s
#           }
#           Limit
#       }
#       if (is.null(x.lim)) {
#           if (phyloORclado) {
#               if (horizontal) {
#                   xx.tips <- xx[1:Ntip]
#                   if (show.tip.label) {
#                     pin1 <- par("pin")[1]
#                     tmp <- getLimit(xx.tips, x$tip.label, pin1,
#                       cex)
#                     tmp <- tmp + label.offset
#                   }
#                   else tmp <- max(xx.tips)
#                   x.lim <- c(0, tmp)
#               }
#               else x.lim <- c(1, Ntip)
#           }
#           else switch(type, fan = {
#               if (show.tip.label) {
#                   offset <- max(nchar.tip.label * 0.018 * max.yy *
#                     cex)
#                   x.lim <- range(xx) + c(-offset, offset)
#               } else x.lim <- range(xx)
#           }, unrooted = {
#               if (show.tip.label) {
#                   offset <- max(nchar.tip.label * 0.018 * max.yy *
#                     cex)
#                   x.lim <- c(0 - offset, max(xx) + offset)
#               } else x.lim <- c(0, max(xx))
#           }, radial = {
#               if (show.tip.label) {
#                   offset <- max(nchar.tip.label * 0.03 * cex)
#                   x.lim <- c(-1 - offset, 1 + offset)
#               } else x.lim <- c(-1, 1)
#           })
#       }
#       else if (length(x.lim) == 1) {
#           x.lim <- c(0, x.lim)
#           if (phyloORclado && !horizontal)
#               x.lim[1] <- 1
#           if (type %in% c("fan", "unrooted") && show.tip.label)
#               x.lim[1] <- -max(nchar.tip.label * 0.018 * max.yy *
#                   cex)
#           if (type == "radial")
#               x.lim[1] <- if (show.tip.label)
#                   -1 - max(nchar.tip.label * 0.03 * cex)
#               else -1
#       }
#       if (phyloORclado && direction == "leftwards")
#           xx <- x.lim[2] - xx
#       if (is.null(y.lim)) {
#           if (phyloORclado) {
#               if (horizontal)
#                   y.lim <- c(1, Ntip)
#               else {
#                   pin2 <- par("pin")[2]
#                   yy.tips <- yy[1:Ntip]
#                   if (show.tip.label) {
#                     tmp <- getLimit(yy.tips, x$tip.label, pin2,
#                       cex)
#                     tmp <- tmp + label.offset
#                   }
#                   else tmp <- max(yy.tips)
#                   y.lim <- c(0, tmp)
#               }
#           }
#           else switch(type, fan = {
#               if (show.tip.label) {
#                   offset <- max(nchar.tip.label * 0.018 * max.yy *
#                     cex)
#                   y.lim <- c(min(yy) - offset, max.yy + offset)
#               } else y.lim <- c(min(yy), max.yy)
#           }, unrooted = {
#               if (show.tip.label) {
#                   offset <- max(nchar.tip.label * 0.018 * max.yy *
#                     cex)
#                   y.lim <- c(0 - offset, max.yy + offset)
#               } else y.lim <- c(0, max.yy)
#           }, radial = {
#               if (show.tip.label) {
#                   offset <- max(nchar.tip.label * 0.03 * cex)
#                   y.lim <- c(-1 - offset, 1 + offset)
#               } else y.lim <- c(-1, 1)
#           })
#       }
#       else if (length(y.lim) == 1) {
#           y.lim <- c(0, y.lim)
#           if (phyloORclado && horizontal)
#               y.lim[1] <- 1
#           if (type %in% c("fan", "unrooted") && show.tip.label)
#               y.lim[1] <- -max(nchar.tip.label * 0.018 * max.yy *
#                   cex)
#           if (type == "radial")
#               y.lim[1] <- if (show.tip.label)
#                   -1 - max(nchar.tip.label * 0.018 * max.yy * cex)
#               else -1
#       }
#       if (phyloORclado && direction == "downwards")
#           yy <- y.lim[2] - yy
#       if (phyloORclado && root.edge) {
#           if (direction == "leftwards")
#               x.lim[2] <- x.lim[2] + x$root.edge
#           if (direction == "downwards")
#               y.lim[2] <- y.lim[2] + x$root.edge
#       }
#       asp <- if (type %in% c("fan", "radial", "unrooted"))
#           1
#       else NA
#       plot.default(0, type = "n", xlim = x.lim, ylim = y.lim, xlab = "",
#           ylab = "", axes = FALSE, asp = asp, ...)
#       if (plot) {
#           if (is.null(adj))
#               adj <- if (phyloORclado && direction == "leftwards")
#                   1
#               else 0
#           if (phyloORclado && show.tip.label) {
#               MAXSTRING <- max(strwidth(x$tip.label, cex = cex))
#               loy <- 0
#               if (direction == "rightwards") {
#                   lox <- label.offset + MAXSTRING * 1.05 * adj
#               }
#               if (direction == "leftwards") {
#                   lox <- -label.offset - MAXSTRING * 1.05 * (1 -
#                     adj)
#               }
#               if (!horizontal) {
#                   psr <- par("usr")
#                   MAXSTRING <- MAXSTRING * 1.09 * (psr[4] - psr[3])/(psr[2] -
#                     psr[1])
#                   loy <- label.offset + MAXSTRING * 1.05 * adj
#                   lox <- 0
#                   srt <- 90 + srt
#                   if (direction == "downwards") {
#                     loy <- -loy
#                     srt <- 180 + srt
#                   }
#               }
#           }
#           if (type == "phylogram") {
#               phylogram.plot(x$edge, Ntip, Nnode, xx, yy, horizontal,
#                   edge.color, edge.width, edge.lty)
#           }
#           else {
#               if (type == "fan") {
#                   ereorder <- match(z$edge[, 2], x$edge[, 2])
#                   if (length(edge.color) > 1) {
#                     edge.color <- rep(edge.color, length.out = Nedge)
#                     edge.color <- edge.color[ereorder]
#                   }
#                   if (length(edge.width) > 1) {
#                     edge.width <- rep(edge.width, length.out = Nedge)
#                     edge.width <- edge.width[ereorder]
#                   }
#                   if (length(edge.lty) > 1) {
#                     edge.lty <- rep(edge.lty, length.out = Nedge)
#                     edge.lty <- edge.lty[ereorder]
#                   }
#                   circular.plot(z$edge, Ntip, Nnode, xx, yy, theta,
#                     r, edge.color, edge.width, edge.lty)
#               }
#               else cladogram.plot(x$edge, xx, yy, edge.color, edge.width,
#                   edge.lty)
#           }
#           if (root.edge) {
#               # rootcol <- if (length(edge.color) == 1)
#               #     edge.color
#               # else "black"
#               rootcol <- root.edge.color
#               # rootw <- if (length(edge.width) == 1)
#               #     edge.width
#               # else 1
#               rootw <- root.edge.width
#               # rootlty <- if (length(edge.lty) == 1)
#               #     edge.lty
#               # else 1
#               rootlty <- root.edge.lty
#
#               if (type == "fan") {
#                   tmp <- polar2rect(x$root.edge, theta[ROOT])
#                   segments(0, 0, tmp$x, tmp$y, col = rootcol, lwd = rootw,
#                     lty = rootlty)
#               }
#               else {
#                   switch(direction, rightwards = segments(0, yy[ROOT],
#                     x$root.edge, yy[ROOT], col = rootcol, lwd = rootw,
#                     lty = rootlty), leftwards = segments(xx[ROOT],
#                     yy[ROOT], xx[ROOT] + x$root.edge, yy[ROOT],
#                     col = rootcol, lwd = rootw, lty = rootlty),
#                     upwards = segments(xx[ROOT], 0, xx[ROOT], x$root.edge,
#                       col = rootcol, lwd = rootw, lty = rootlty),
#                     downwards = segments(xx[ROOT], yy[ROOT], xx[ROOT],
#                       yy[ROOT] + x$root.edge, col = rootcol, lwd = rootw,
#                       lty = rootlty))
#               }
#           }
#           if (show.tip.label) {
#               if (is.expression(x$tip.label))
#                   underscore <- TRUE
#               if (!underscore)
#                   x$tip.label <- gsub("_", " ", x$tip.label)
#               if (phyloORclado) {
#                   if (align.tip.label) {
#                     xx.tmp <- switch(direction, rightwards = max(xx[1:Ntip]),
#                       leftwards = min(xx[1:Ntip]), upwards = xx[1:Ntip],
#                       downwards = xx[1:Ntip])
#                     yy.tmp <- switch(direction, rightwards = yy[1:Ntip],
#                       leftwards = yy[1:Ntip], upwards = max(yy[1:Ntip]),
#                       downwards = min(yy[1:Ntip]))
#                     segments(xx[1:Ntip], yy[1:Ntip], xx.tmp, yy.tmp,
#                       lty = align.tip.label.lty)
#                   }
#                   else {
#                     xx.tmp <- xx[1:Ntip]
#                     yy.tmp <- yy[1:Ntip]
#                   }
#                   text(xx.tmp + lox, yy.tmp + loy, x$tip.label,
#                     adj = adj, font = font, srt = srt, cex = cex,
#                     col = tip.color)
#               }
#               else {
#                   angle <- if (type == "unrooted")
#                     XY$axe
#                   else atan2(yy[1:Ntip], xx[1:Ntip])
#                   lab4ut <- if (is.null(lab4ut)) {
#                     if (type == "unrooted")
#                       "horizontal"
#                     else "axial"
#                   }
#                   else match.arg(lab4ut, c("horizontal", "axial"))
#                   xx.tips <- xx[1:Ntip]
#                   yy.tips <- yy[1:Ntip]
#                   if (label.offset) {
#                     xx.tips <- xx.tips + label.offset * cos(angle)
#                     yy.tips <- yy.tips + label.offset * sin(angle)
#                   }
#                   if (lab4ut == "horizontal") {
#                     y.adj <- x.adj <- numeric(Ntip)
#                     sel <- abs(angle) > 0.75 * pi
#                     x.adj[sel] <- -strwidth(x$tip.label)[sel] *
#                       1.05
#                     sel <- abs(angle) > pi/4 & abs(angle) < 0.75 *
#                       pi
#                     x.adj[sel] <- -strwidth(x$tip.label)[sel] *
#                       (2 * abs(angle)[sel]/pi - 0.5)
#                     sel <- angle > pi/4 & angle < 0.75 * pi
#                     y.adj[sel] <- strheight(x$tip.label)[sel]/2
#                     sel <- angle < -pi/4 & angle > -0.75 * pi
#                     y.adj[sel] <- -strheight(x$tip.label)[sel] *
#                       0.75
#                     text(xx.tips + x.adj * cex, yy.tips + y.adj *
#                       cex, x$tip.label, adj = c(adj, 0), font = font,
#                       srt = srt, cex = cex, col = tip.color)
#                   }
#                   else {
#                     if (align.tip.label) {
#                       POL <- rect2polar(xx.tips, yy.tips)
#                       POL$r[] <- max(POL$r)
#                       REC <- polar2rect(POL$r, POL$angle)
#                       xx.tips <- REC$x
#                       yy.tips <- REC$y
#                       segments(xx[1:Ntip], yy[1:Ntip], xx.tips,
#                         yy.tips, lty = align.tip.label.lty)
#                     }
#                     if (type == "unrooted") {
#                       adj <- abs(angle) > pi/2
#                       angle <- angle * 180/pi
#                       angle[adj] <- angle[adj] - 180
#                       adj <- as.numeric(adj)
#                     }
#                     else {
#                       s <- xx.tips < 0
#                       angle <- angle * 180/pi
#                       angle[s] <- angle[s] + 180
#                       adj <- as.numeric(s)
#                     }
#                     font <- rep(font, length.out = Ntip)
#                     tip.color <- rep(tip.color, length.out = Ntip)
#                     cex <- rep(cex, length.out = Ntip)
#                     for (i in 1:Ntip) text(xx.tips[i], yy.tips[i],
#                       x$tip.label[i], font = font[i], cex = cex[i],
#                       srt = angle[i], adj = adj[i], col = tip.color[i])
#                   }
#               }
#           }
#           if (show.node.label)
#               text(xx[ROOT:length(xx)] + label.offset, yy[ROOT:length(yy)],
#                   x$node.label, adj = adj, font = font, srt = srt,
#                   cex = cex)
#       }
#       L <- list(type = type, use.edge.length = use.edge.length,
#           node.pos = node.pos, node.depth = node.depth, show.tip.label = show.tip.label,
#           show.node.label = show.node.label, font = font, cex = cex,
#           adj = adj, srt = srt, no.margin = no.margin, label.offset = label.offset,
#           x.lim = x.lim, y.lim = y.lim, direction = direction,
#           tip.color = tip.color, Ntip = Ntip, Nnode = Nnode, root.time = x$root.time,
#           align.tip.label = align.tip.label)
#       assign("last_plot.phylo", c(L, list(edge = xe, xx = xx, yy = yy)),
#           envir = .PlotPhyloEnv)
#       invisible(L)
# }

#' get a densiTree plot from a set of opentree_chronograms
#' if densiTree plot function throws an error, it chooses the tree with the most tips as consensus (using get_biggest_phylo)
#' we found that densiTree errors commonly from failing to do a consensus tree.
#' @inheritParams get_biggest_phylo
#' @param include_all Boolean. Default to FALSE, exclude chronograms with only two species. If TRUE, it plots all chronograms.
#' @inheritDotParams phangorn::densiTree -x -consensus
#' @export
plot_densitree <- function(trees, include_all = FALSE, ...){
  if (!include_all){
    trees <- trees[sapply(trees, function (x) ape::Ntip(x)) > 2]
  }
  class(trees) <- "multiPhylo"
  tryCatch(phangorn::densiTree(x = trees, ...), error = function(e) {
    biggest_phylo <- get_biggest_phylo(trees = trees)
    phangorn::densiTree(x = trees, consensus = biggest_phylo, ...)
  })
}
