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

#' A data frame containing the stratigraphic chart by issued in 2012 by the International Commission on Stratigraphy.
#'
#' @name strat2012
#' @docType data
#' @format A data frame
#' \describe{
#'   \item{eon}{Eon name.}
#'   \item{era}{Era name.}
#'   \item{period}{Period name.}
#'   \item{epoch}{Epoch name.}
#'   \item{stage}{Stage name.}
#'   \item{MA}{Estimated boundary age for the associated interval.}
#'   \item{error}{Estimated errors associated with each age estimate.}
#'   \item{GSSP}{Binary response denoting whether the age estimate is defined by a basal Global Standard Section and Point}
#' }
#' @source \url{https://github.com/fmichonneau/phyloch}
#' @keywords plot gradstein geochronology
#' @details
#' See phyloch package data description.
#' Enhance: There are a couple errors in the original strat2012 object that are not too meaningful but could be improved.
#' Generated with utils::data("strat2012", package = "phyloch")
#' use_data(strat2012)
"strat2012"

#' subset trees for plotting densitree plots and phylo_all plots
#' @inheritParams get_biggest_phylo
#' @param include Boolean or numeric vector. Default to TRUE, keep all chronograms
#' in trees. If FALSE, exclude chronograms with only two tips. If numeric, it is used
#' as indices to subset trees object.
subset_trees <- function(trees, include = TRUE){
  if(is.numeric(include)){  #if it is numeric
    include <- round(include)
    include <- include[which(include <= length(trees))]
    include <- sort(unique(include))
    trees <- trees[include]
  }
  if(is.logical(include) & !is.na(include)){
      if (!include){
      trees <- trees[sapply(trees, function (x) ape::Ntip(x)) > 2]
    }
  }
  trees
}

#' get a densiTree plot from a set of opentree_chronograms
#' if densiTree plot function throws an error, it chooses the tree with the most tips as consensus (using get_biggest_phylo)
#' we found that densiTree errors commonly from failing to do a consensus tree.
#' @inheritParams get_biggest_phylo
#' @inheritParams subset_trees
#' @inheritDotParams phangorn::densiTree -x -consensus
#' @export
plot_densitree <- function(trees, include = TRUE, ...){
  trees <- subset_trees(trees, include = include)
  # for densitree plot it does not really matter if trees have the same length
  # max.depth <- round(max(sapply(trees, function(x) max(ape::branching.times(x)))) + 5, digits = -1)
  # # max.depth <- round(max(summarize_datelife_result(datelife_result = get.filtered.results(),
  # #             summary_format = "mrca")) + 5)  # in case we used it for a datelife object
  #
  # # plot all trees from the same depth:
  # trees <- lapply(trees, function(x) {
  #  x$root.edge <- max.depth - max(ape::branching.times(x))
  #  x
  # })
  class(trees) <- "multiPhylo"
  # if we use biggest phylo as consensus for all, some data set are not plotted correctly
  # biggest_phylo <- get_biggest_phylo(trees = trees)
  # try(phangorn::densiTree(x = trees, consensus = biggest_phylo, ...))
  tryCatch(phangorn::densiTree(x = trees, ...),
  error = function(e) {
    biggest_phylo <- get_biggest_phylo(trees = trees)
    try(phangorn::densiTree(x = trees, consensus = biggest_phylo, ...))
  })
}

#' get the outer margin of a graphics device from a number of tips
#' @param tree phylo object to be plotted
phylo_height_omi <- function(tree){
  tipnum <- ape::Ntip(tree)
  if(tipnum > 10){
    hei <- 50 + (30 * tipnum)
  } else {
    hei <- 300
  }
  if(tipnum == 2){
    omi1 <- 2
  } else if(tipnum == 3){
    omi1 <- 1.5
  } else if(tipnum == 4){
    omi1 <- 1.2
  } else if (tipnum >= 5 & tipnum <= 7){
    # omi1 <- 3
    omi1 <- 1 # this works good for pdf images
  } else if (tipnum >= 8 & tipnum <= 10){
    omi1 <- 2.5
  } else {
    omi1 <- 2
  }
  return(list(height = hei, omi1 = omi1))

}

#' wrap a character string to a plotting area. It gives the optimal cex and width. It works once plot device has been called
#' idea to use strwrap from https://stackoverflow.com/questions/7367138/text-wrap-for-plot-titles
#'
#' @param string A character vector with the text to be plotted
#' @param max_cex A real number, giving the maximum *cex* (**c**haracter **ex**pansion) for the string to be plotted
#' @param min_cex minimum character expansion to be used on the title
#' @param string_font font type to be used on the title
#' @param max_height A real number, giving the maximum height to be covered by the text
#' @param init_strwrap_width A real number indicating the minimum number of characters to be plotted by line
#' @param min_width A real number, giving the minimum width to be occupied by the string
#' @param max_width A real number, giving the maximum width to be occupied by the string
#' @param whole Boolean, indicating if the whole string should be plotted even if it surpasses the limits established in previous arguments
#' @export
wrap_string_to_plot <- function(string, max_cex = 1, min_cex = 0.5, string_font = 2,
                           max_height = graphics::par("din")[2]-graphics::par("pin")[2]- graphics::par("omi")[1]-graphics::par("mai")[1] - 0.2,
                           init_strwrap_width = 50,  min_width = 0.9*graphics::par("din")[1],
                           max_width = 0.9*graphics::par("din")[1], whole = TRUE){
  # collapse string to a vetor of one element in case it has more elements
  if(max_height <= 0){
    message("string cannot be adjusted because there is not enough space on upper margin of plotting device")
    return(NA)
  }
  wraps <- strwrap(string, width = init_strwrap_width)
  sw <- graphics::strwidth(s = wraps, units = "inches", cex = max_cex, font = string_font)
  sh <- graphics::strheight(s = paste(wraps, collapse = "\n"), units = "inches", cex = max_cex, font = string_font)
  wi <- init_strwrap_width - 1
  string_cex <- max_cex
  # next if when title is too short and fits in one line with regular cex:
  if(sh < max_height & max(sw) < max_width){ #if(length(wraps) ==1)
      return(list(wrapped = paste(wraps, collapse = "\n"), wraps = wraps,
      string_cex = string_cex, strwrap_width = wi, string_font = string_font))
  }
  # next while to find the appropriate cex to fit a long title with a max and min width:
  while(sh > max_height | max(sw) > max_width | max(sw) < min_width) {  #max(sw) < min_width | max(sw) > max_width
    while(max(sw) < min_width) {
      wi <- wi + 1
      wraps <- strwrap(string, width = wi)
      sw <- graphics::strwidth(s = wraps, units = "inches", cex = string_cex, font = string_font)
    }
    sh <- graphics::strheight(s = paste(wraps, collapse = "\n"), units = "inches", cex = string_cex, font = string_font)
    if(sh > max_height){ #length(wraps) > n_lines |
      if(string_cex <= min_cex){
        break
      }
      string_cex <- string_cex - 0.01
      sw <- graphics::strwidth(s = wraps, units = "inches", cex = string_cex, font = string_font)
      sh <- graphics::strheight(s = paste(wraps, collapse = "\n"), units = "inches", cex = string_cex, font = string_font)
    } else {
      break
    }
  }
  if(!whole){
    while(sh > max_height){
      wraps <- wraps[-length(wraps)]
      sh <- graphics::strheight(s = paste(wraps, collapse = "\n"), units = "inches", cex = string_cex, font = string_font)
    }
  }
  return(list(wrapped = paste(wraps, collapse = "\n"), wraps = wraps,
  string_cex = string_cex, strwrap_width = wi, string_font = string_font))
}

#' plot all trees with study titles and geochronological axis
#'
#' @inheritParams get_biggest_phylo
#' @inheritParams subset_trees
#' @inheritParams plot_phylo
#' @param individually Boolean indicating if trees should be plotted one by one or all on the same file
#' @inheritParams ape::plot.phylo
#' @param file A character string giving the name and path to write the files to.
#' @export
plot_phylo_all <- function(trees, cex = graphics::par("cex"), include = TRUE, individually = TRUE, write = "no", file = "phylo_all"){
  trees <- subset_trees(trees, include = include)
  if(any("tip.label" %in% names(trees))){ # in case there is just one tree in trees
    trees <- list(trees)
  }
  # if(isTRUE(all.equal(round(sapply(trees, function(x) max(ape::branching.times(x))), digits = 3))))
  max.depth <- round(max(sapply(trees, function(x) max(ape::branching.times(x)))) + 5, digits = -1)
  # max.depth <- round(max(summarize_datelife_result(datelife_result = get.filtered.results(),
  #             summary_format = "mrca")) + 5)  # in case we used it for a datelife object

  # plot all trees from the same depth:
  trees <- lapply(trees, function(x) {
   x$root.edge <- max.depth - max(ape::branching.times(x))
   x
  })
  mai4 <- unique(unlist(sapply(trees, "[", "tip.label")))
  ind <- which.max(nchar(mai4))
  mai4 <- graphics::strwidth(s = mai4[ind], units = "inches", cex = cex, font = 3)
  # if(any(lapply(trees, ape::Ntip) > 3))
  # png("~/tmp/axisgeo.png", units = "in")
  if(!any(c("png", "pdf") %in% write)){
  # if (!grDevices::devAskNewPage() && !names(grDevices::dev.cur()) %in% c("pdf", "postscript")) {
      grDevices::devAskNewPage(TRUE)
      # dev.size("px")
      on.exit(grDevices::devAskNewPage(FALSE))
  } else {
      dir.create(path = gsub("\\.png$|\\.pdf$", "", file))
  }
  for (i in 1:length(trees)){
    file_name <- paste0(gsub("\\.png$|\\.pdf$", "", file), "/", gsub("\\.png$|\\.pdf$", "", file), "_", i, ".", write)
    plot_phylo(trees[[i]], names(trees)[i], time_depth = max.depth, axis_type = 1, cex, mai4, write, file_name, GTS = NULL)
  }
  # getting an "unrecoverable" error with merge PDF:
  # if(!individually){
  #   # install_github("trinker/plotflow")
  # if we decide to use this, we should add plotflow functions in datelife package so we don't have to add it to description...
  #   plotflow:::mergePDF(
  #       in.file= paste(file.path(gsub("\\.png$|\\.pdf$", "", file), dir(gsub("\\.png$|\\.pdf$", "", file)))),
  #       file= paste0(gsub("\\.png$|\\.pdf$", "", file), ".", write)
  #       # file= "merged.pdf"
  #   )
  # }
}
#' plot one tree with study title and geochronological axis
#'
#' @inheritParams get_biggest_phylo
#' @inheritParams tree_fix_brlen
#' @param title A character string giving the name and path to write the files to.
#' @param time_depth A numeric vector indicating the upper limit on the time x axis scale.
#' @param axis_type A numeric vector indicating the type of geochronological axis to plot. See examples.
#' @param mai4 A numeric vector of length one indicating the space needed for plotting whole tip labels (right margin of the plot).
#' @param write A character vector of length 1. Use pdf or png to write a file on those formats respectively. Anything else will not write any image file.
#' @inheritParams ape::plot.phylo
#' @param file_name A character string giving the name and path to write the files to.
#' @param GTS A dataframe of geochronological limits.
#' @export
# enhance: examples of axis_types!
plot_phylo <- function(tree, title = "Tree", time_depth = NULL, axis_type = 1,
cex = graphics::par("cex"), mai4 = NULL, write = "nothing", file_name = NULL, GTS = getAnywhere("strat2012")){
  if(is.null(GTS)){
    # utils::data(strat2012)
    GTS <- getAnywhere("strat2012")
  }
  if(is.null(time_depth) & !is.null(tree$edge.length)){
    if(is.null(tree$root.edge)){
      time_depth <- round(max(ape::branching.times(tree)) + 5, digits = -1)
    } else {
      time_depth <- max(ape::branching.times(tree)) + tree$root.edge
    }
  }
  if(is.null(mai4)){
    ind <- which.max(nchar(tree$tip.label))
    mai4 <- graphics::strwidth(s = tree$tip.label[ind], units = "inches", cex = cex, font = 3)
  }
  ho <- phylo_height_omi(tree = tree)
  if(any(c("png", "pdf") %in% write)){
    if("png" %in% write){
      grDevices::png(file = file_name, height = ho$height)
    } else {
      grDevices::pdf(file = file_name, height = ho$height/72)
    }
  }
  graphics::par(xpd = NA, mai = c(0,0,0,mai4), omi = c(ho$omi1, 0, 1, 0))
  # plot_chronogram.phylo(trees[[i]], cex = 1.5, edge.width = 2, label.offset = 0.5,
    # x.lim = c(0, max.depth), root.edge = TRUE, root.edge.color = "white")
  if(is.null(tree$edge.length)){
      ape::plot.phylo(tree, cex = cex, #edge.width = 2,
        label.offset = 0.5, plot = TRUE)  #, ...
  } else {
      ape::plot.phylo(tree, cex = cex, #edge.width = 2,
        label.offset = 0.5, x.lim = c(0, time_depth), root.edge = TRUE, plot = TRUE)  #, ...
  }
  graphics::par(xpd = FALSE)
  if(!is.null(tree$edge.length)){
      if(axis_type == 1){
        axisGeo(GTS = GTS, unit = c("period"),
            col = c("gray80", "white"), gridcol = c("gray80", "white"), cex = 0.5,
            gridty = "twodash")
      }
      if(axis_type == 2){
        axisGeo(GTS = GTS, unit = c("period","epoch"),
            col = c("gray80", "white"), gridcol = c("gray80", "white"), cex = 0.5,
            gridty = "twodash")
      }
      if(axis_type == 3){
        strap::geoscalePhylo(tree=tree, ages=tree$ranges.used, cex.tip=0.7,
          cex.ts=0.7,cex.age=0.7, width=4, tick.scale = 15, boxes = "Epoch", erotate = 90,
          quat.rm=TRUE, units=c("Period","Epoch"), x.lim=c(65,-10))
      }
      graphics::mtext("Time (MYA)", cex = cex, side = 1, font = 2, line = (ho$omi1-0.2)/0.2,
      outer = TRUE, at = 0.4)
  } else (
      message("tree has no edge.length, so time axis will not be plotted")
  )
  if(!is.null(title)){
    titlei <- wrap_string_to_plot(string = title, max_cex = 1, whole = FALSE)
    graphics::mtext(text = titlei$wrapped, outer = TRUE,
      cex = titlei$string_cex, font = titlei$string_font, line = 1)
  }
  if(any(c("png", "pdf") %in% write)){
    grDevices::dev.off()
  }
}
# tree <- plant_bold_otol_tree
# plot_phylo_gg <- function(tree, title = "Tree", time_depth = NULL, axis_type = 1,
# cex = graphics::par("cex"), mai4 = NULL, write = "nothing", file_name = NULL, GTS = getAnywhere("strat2012")){
#   max_age <- max(ape::branching.times(tree))
#   age_lim <- max_age*1.2
#   grDevices::pdf("test.pdf")
#   p <- ggtree::ggtree(tree) + ggtree::geom_tiplab()  + #ggplot2::xlim(age_lim*0.1,-age_lim) +
#   ggplot2::coord_cartesian(xlim = c(age_lim*0.5,-age_lim), ylim = c(-1, ape::Ntip(tree)), expand = FALSE) +
#   ggplot2::scale_x_continuous(breaks=seq(-age_lim,0,100), labels=abs(seq(-age_lim,0,100))) +
#   ggtree::theme_tree2()
#   p <- ggtree::revts(p)
#   deeptime::gggeo_scale(p, neg = TRUE)
#   print(p)
#   grDevices::dev.off()
# }
# .PlotPhyloEnv <- new.env()
