# Functions from `phyloch` package used in datelife
# Hosting phyloch functions here for now so we can submit to CRAN
# phyloch on github at https://github.com/fmichonneau/phyloch


#' Function from phyloch package
#'  Add interval bars to nodes on a phylo plot.
#' @param phy A phylo object with two extra elements specifying the MAX and MIN limits of intervals to plot.
#' @param label A character vector corresponding to the name used to label intervals stored in phy
#' @param tab Something I'm not sure of, should ask Christoph Heibl.
#' @param nodes Numeric. Use this if you wanna plot bars on certain nodes only.
#' @param col Color for the bars, default to "grey75"
#' @param lwd Line width for the bars, default to 5.
#' @param broken Something I'm not sure of, should ask Christoph Heibl.
#' @inheritDotParams graphics::segments
#' @importFrom ape .PlotPhyloEnv
#' @return Adds bars to nodes on a phylo plot.
#' @export
#' @details Make sure ape is loaded otherwise it won't find .PlotPhyloEnv
#' @author Christoph Heibl
HPDbars <- function(phy, label = "height_95%_HPD", tab = NULL, nodes, col, lwd, broken = FALSE, ...){

	# check input data
	# ----------------
	if (!inherits(phy, "phylo"))
        stop("object 'phy' is not of class 'phylo'")
    if (!label %in% gsub("_MIN|_MAX", "", names(phy)))
        stop("there is no element '", label, "' in object 'phy'")
    if (!paste(label, "MIN", sep = "_") %in% names(phy))
        stop("there is no lower bound for '", label,
            "' in 'phy': the corresponding element must be named '", 			label, "_MIN'")
    if (!paste(label, "MAX", sep = "_") %in% names(phy))
        stop("there is no upper bound for '", label,
            "' in 'phy': the corresponding element must be named '", 			label, "_MAX'")
	lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
	if (!lastPP$use.edge.length)
	    stop("function needs edge length information")
	if (lastPP$type != "phylogram")
	    stop("currently only 'type == phylogram' supported")

	# set default parameters
	# ----------------------
	if (missing(lwd)) lwd <- 5
	if (missing(col)) col <- "grey75"

	if (missing(nodes)) nodes <- lastPP$Ntip + 1:lastPP$Nnode
	nodes <- sort(nodes)

	# get x and y values:
	# -------------------
	label <- paste(label, c("MIN", "MAX"), sep = "_")
	if (is.null(tab)){
		# label %in% names(phy)
		# phy <- subset2_sdmphylo_mean
		# label <- "height_95%_HPD"
		id <- which(names(phy) %in% label)
		# previous line was not working anymore so modified it:
		x <- cbind(phy[[id[1]]],
							 phy[[id[2]]])
		# length(nodes - lastPP$Ntip)
		# enhance: add a tryCatch here, it will tell when there are less nodes in calibrations than actual nodes in phy
		if(nrow(x) != length(nodes)){
			x <- x[nodes - lastPP$Ntip, ]
		}
	}	else {
		id <- which(colnames(tab) %in% label)
		x <- tab[tab[, 1] %in% nodes, id]
	}
	phy$node.label <- NULL # important: if phy already has node labels the next two lines won't work
	bt <- ape::branching.times(phy)
	bt <- bt[names(bt) %in% nodes]
	x <- cbind(x, bt)
	x <- cbind(x[, 3] - x[, 1], x[, 2] - x[, 3], x)
	colnames(x) <- c("min", "max", "MIN", "MAX", "mean")

	if (lastPP$dir %in% c("upwards", "downwards")){
		y <- x
		x <- lastPP$xx[nodes]
	}
	else y <- lastPP$yy[nodes]

	# plot bars:
	# -------------
	if (lastPP$dir == "rightwards"){
		xx <- lastPP$xx[nodes]
		xx <- cbind(xx - x[, 2], xx + x[, 1])
		id <- which(xx[, 1] < lastPP$x.lim[1])
		if (length(id) > 0) {
			ns <- paste(nodes[id], collapse = ", ")
			minx <- min(xx[id, 1])
			new.xlim <- c(minx, lastPP$x.lim[2] - minx)
			new.xlim <- round(new.xlim, digits = 5)
			new.xlim <- paste(new.xlim, collapse = ", ")
			message("HPD bar(s) for nodes ", ns,
			  " exceed(s) plot region.",
			  "\n  Try setting \"x.lim\" to c(", new.xlim, ")")
			graphics::segments(xx[id, 1], y[id],
			    lastPP$xx[lastPP$Ntip + id], y[id],
		        col = col, lwd = lwd, lend = 1, lty = "11")
		    graphics::segments(lastPP$xx[lastPP$Ntip + id], y[id],
			    xx[id, 2], y[id],
		        col = col, lwd = lwd, lend = 0)
		    xx <- xx[-id, ]; y <- y[-id]
		}
		graphics::segments(xx[, 1], y, xx[, 2], y,
		    col = col, lwd = lwd, ...)
    }
    if (lastPP$dir == "leftwards"){
    	xx <- lastPP$xx[nodes]
    	xx <- cbind(xx - x[, 1], xx + x[, 2])
		id <- which(xx[, 2] > lastPP$x.lim[2])
		if (length(id) > 0) {
			ns <- paste(nodes[id], collapse = ", ")
			warning("HPD bar(s) for nodes ", ns,
			  " exceed(s) plot region (issue currently unsolved)")
		}
	    graphics::segments(xx[, 1], y, xx[, 2], y,
		    col = col, lwd = lwd, ...)
	}
	if (lastPP$dir == "upwards"){
		yy <- lastPP$yy[nodes]
		yy <- cbind(yy - y[, 2], yy + y[, 1])
		id <- which(yy[, 1] < lastPP$y.lim[1])
		if (length(id) > 0) {
			ns <- paste(nodes[id], collapse = ", ")
			miny <- min(yy[id, 1])
			new.ylim <- c(miny, lastPP$y.lim[2] - miny)
			new.ylim <- round(new.ylim, digits = 5)
			new.ylim <- paste(new.ylim, collapse = ", ")
			warning("HPD bar(s) for nodes ", ns,
			  " exceed(s) plot region.",
			  "\n  Try setting \"y.lim\" to c(", new.ylim, ")")
		    graphics::segments(x[id], yy[id, 1],
			    x[id], lastPP$yy[lastPP$Ntip + id],
		        col = col, lwd = lwd, lend = 1, lty = "11")
		    graphics::segments(x[id], lastPP$yy[lastPP$Ntip + id],
			    x[id], yy[id, 2],
		        col = col, lwd = lwd, lend = 0)
		    yy <- yy[-id, ]; x <- x[-id]
		}
	    graphics::segments(x, yy[, 1], x, yy[, 2],
	        col = col, lwd = lwd,  ...)
	}
	if (lastPP$dir == "downwards"){
		yy <- lastPP$yy[nodes]
		yy <- cbind(yy - y[, 1], yy + y[, 2])
		id <- which(yy[, 2] > lastPP$y.lim[2])
		if (length(id) > 0){
			ns <- paste(nodes[id], collapse = ", ")
			warning("HPD bar(s) for node(s) ", ns,
			  " exceed(s) plot region (issue currently unsolved)")

		}
        graphics::segments(x, yy[, 1], x, yy[, 2],
            col = col, lwd = lwd,  ...)
	}
}

#' Function from Christoph Heibl's phyloch package
#' Add a Geological Time Axis for Time-Calibrated Phylogenies.
#' @param GTS a data frame containing the geological time scale such as provided by \code{\link{strat2012}}.
#' @param tip.time a real number giving time units back from present were the youngest taxa were sampled.
#' @param unit a vector of mode character used to select geological time units that shall be displayed. When using \code{gradstein04}, \code{"eon"}, \code{"era"}, \code{"period"}, \code{"epoch"}, and \code{"stage"} are available.
#' @param ages logical: if \code{TRUE}, a real-numbered axis of unit boundaries is added to the plot.
#' @param cex a real number giving the \bold{c}haracter \bold{ex}pansion.
#' @param col a vector of mode character used to select background color(s). Will be recycled if necessary.
#' @param texcol a vector of mode character used to select font color(s). Will be recycled if necessary.
#' @param gridty an integer representing the line type of the grid exactly in the same way as \code{lty} argument in \code{\link{par}}: 0 or \code{blank}, 1 or \code{"solid"}, 2 or \code{"dashed"}, 3 or \code{"dotted"}, 4 or \code{"dotdash"}, 5 or \code{"longdash"}, and 6 or \code{"twodash"}.
#' @param gridcol a vector of mode character used to select grid color(s). Will be recycled if necessary.
#' @importFrom ape .PlotPhyloEnv
#' @return Adds geological time scale to a phylo plot.
#' @export
#' @details
#' Make sure ape is loaded otherwise it won't find .PlotPhyloEnv
#' If the name of a time interval does not fit the space on the plotting device provided by that time frame, \code{axisGeo} tries to fit the name string by lowering the original character expansion. In order to avoid undiscernable small font sizes only character expansion of 50 percent or more of the original \code{cex} is allowed; in all other cases strings are abbreviated and plotted with the original character expansion.
#' @author Christoph Heibl
#' @references
#' Gradstein F.M., J.G. Ogg & A.G. Smith. 2004. A Geologic Time Scale. Cambridge University Press, Cambridge, UK. \url{www.cambridge.org/uk/catalogue/catalogue.asp?isbn=0521786738}.
#' International Commission on Stratigraphy. 2012. International stratigraphic chart: \url{http://www.stratigraphy.org/upload/ISChart2009.pdf}.

axisGeo <- function(GTS, tip.time = 0, unit = c("epoch", "period"), ages =
  TRUE, cex = 1, col = "white", texcol = "black", gridty = 0,
                    gridcol = "black"){

  adjustCex <- function(space, string, cex){
    while (graphics::strwidth(string, cex = cex) >= space & cex > .001)
      cex <- cex - .001
    cex
  }

  lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  ntips <- lastPP$Ntip
  root <- ntips + 1

  if ( lastPP$direction == "rightwards" )
    maxage <- max(lastPP$xx) + tip.time
  if ( lastPP$direction == "upwards" )
    maxage <- max(lastPP$yy) + tip.time

  # Geological Time Scale:
  # ----------------------
  gts <- GTS
  maid <- grep("MA", names(gts))
  gts <- cbind(gts[, 1:maid], c(0, utils::head(gts[, maid], -1)), gts[, ((maid + 1):dim(gts)[2])])
  names(gts)[maid:(maid + 1)] <- c("fromMA", "toMA")
  if ( sum(gts[1, maid:(maid + 1)]) == 0 )
    gts <- gts[-1, ]

  # Select corresponding window
  # -----------------------------------------------------------
  ind <- which(gts$fromMA <= maxage & gts$toMA >= tip.time)
  gts <- gts[c(ind, max(ind) + 1), ]
  # cut ends of scale:
  gts$toMA[1] <- tip.time
  gts$fromMA[dim(gts)[1]] <- maxage

  # allow plot to extent into outer margin of device region
  # -------------------------------------------------------
  graphics::par(xpd = NA)

  #
  # ----------------------------
  plotGeo <- function(gts, unit, yy){

    id <- which(names(gts) %in% c(unit, "fromMA", "toMA"))
    gts <- gts[id]
    names(gts) <- c("unit", "from", "to")

    stages <- unique(gts$unit)

    # manage colors
    # --------------
    if (length(col) == 1) col <- rep(col, 2)
    col1 <- rep(col, length(stages))
    col1 <- utils::head(col1, length(col1) / 2)
    if (length(texcol) == 1) texcol <- rep(texcol, 2)
    col2 <- rep(texcol, length(stages))
    col2 <- utils::head(col2, length(col2) / 2)

    xgrid <- NULL
    for (i in seq(along = stages)){
      message(paste0("\nStage ", i, ": ", stages[i], sep = ""))

      # plot boxes of stages:
      # ---------------------
      from <- maxage - max(gts[gts$unit == stages[i], 2])
      to <- maxage - min(gts[gts$unit == stages[i], 3])
      graphics::rect(from, yy[1], to, yy[2], col = col1[i],
           border = "black", lwd = 0.5)
      xgrid <- c(xgrid, from, to)

      # plot names of stages:
      # ---------------------
      en <- as.character(stages[i])
      if ((to - from) > graphics::strwidth(en, cex = cex))
        graphics::text(mean(c(from, to)), mean(yy), en, 					cex = cex, col = col2[4])				else {
          thiscex <- adjustCex(to - from, en, cex) * 0.95
          if (2 * thiscex >= cex)
            graphics::text(mean(c(from, to)), mean(yy), en, 						cex = thiscex, col = col2[i])				else {
              while (nchar(en) > 0 & 								graphics::strwidth(en, cex = cex) >= (to - from))
                en <- paste(utils::head(unlist(strsplit(en, "")), 					-1), collapse = "")
              if (nchar(en) > 1)
                en <- paste(paste(utils::head(unlist(strsplit(						en, "")), -1), collapse = ""), ".", 						sep = "")
              graphics::text(mean(c(from, to)), mean(yy), en, 						cex = cex, col = col2[i])
            }
        }
    }
    xgrid
  } # end of plotGeo

  #
  # ----------------------------
  plotGeoToLeft <- function(gts, unit, yy){

    id <- which(names(gts) %in% c(unit, "fromMA", "toMA"))
    gts <- gts[id]
    names(gts) <- c("unit", "from", "to")

    stages <- unique(gts$unit)

    # manage colors
    # --------------
    if (length(col) == 1) col <- rep(col, 2)
    col1 <- rep(col, length(stages))
    col1 <- utils::head(col1, length(col1) / 2)
    if (length(texcol) == 1) texcol <- rep(texcol, 2)
    col2 <- rep(texcol, length(stages))
    col2 <- utils::head(col2, length(col2) / 2)

    xgrid <- NULL
    for (i in seq(along = stages)){
      from <- maxage - max(gts[gts$unit == stages[i], 2])
      to <- maxage - min(gts[gts$unit == stages[i], 3])
      graphics::rect(yy[2], from, yy[1], to, col = col1[i], 				border = "black", lwd = 0.5)
      xgrid <- c(xgrid, from, to)

      # graphics::text labels
      en <- as.character(stages[i])
      yxr <- (max(lastPP$y.lim) - min(lastPP$y.lim)) / 				(max(lastPP$x.lim) - min(lastPP$x.lim)) * 1.5
      if ((to - from) > graphics::strwidth(en, cex = cex * yxr))
        graphics::text(mean(yy), mean(c(from, to)), en, 					cex = cex, col = col2[4], srt = 90)			else {
          asp <- (to - from) / yxr
          thiscex <- adjustCex(asp, en, cex) * 0.95
          if (1.5 * thiscex >= cex)
            graphics::text(mean(yy), mean(c(from, to)), en, 						cex = thiscex, col = col2[i], srt = 90)				else{
              while (nchar(en) > 0 & 								graphics::strwidth(en, cex = cex * yxr) 						>= (to - from))
                en <- paste(utils::head(unlist(strsplit(en, "")), 					-1), collapse = "")
              if (nchar(en) > 1)
                en <- paste(paste(utils::head(unlist(strsplit(						en, "")), -1), collapse = ""), ".", 						sep = "")
              graphics::text(mean(yy), mean(c(from, to)),  en, 						cex = cex, col = col2[i], srt = 90)
            }
        }
    }
    xgrid
  } # end of plotGeoLeft

  bh <- - graphics::strheight("Ap", cex = cex) * 1.5
  if ( lastPP$direction == "rightwards" ){
    if ( ages ) yy <- c(bh, 2 * bh) else yy <- c(0, bh)
    for ( j in seq_along(unit) ){
      message(paste0("\nPlot unit:", unit[j]))
      if ( j == 1 ) xgrid <- plotGeo(gts, unit[j], yy)
      else plotGeo(gts, unit[j], yy)
      yy <- yy + bh
    }
  }

  if (lastPP$direction == "upwards"){
    if (ages) yy <- c(bh, 2*bh) else yy <- c(0, bh)
    for (j in seq(along = unit)){
      if (j == 1)
        xgrid <- plotGeoToLeft(gts, unit[j], yy)
      else plotGeoToLeft(gts, unit[j], yy)
      yy <- yy + bh
    }
  }

  # draw grid and ages:
  # -------------------
  xgrid <- unique(sort(xgrid, decreasing = TRUE))
  label <- maxage -xgrid

  id <- TRUE
  for (k in seq(along = xgrid)){
    if (lastPP$direction == "rightwards")
      graphics::lines(rep(xgrid[k], 2), c(0, ntips + 1),
            lty = gridty, col = gridcol)
    if (lastPP$direction == "upwards")
      graphics::lines(c(0, ntips + 1), rep(xgrid[k], 2),
            lty = gridty, col = gridcol)

    if (k < length(xgrid)) {
      spneeded <- graphics::strwidth(label[k], cex = cex * 0.8)/2
      spavailable <- xgrid[k] - xgrid[k + 1]
      if (spavailable  < spneeded * 1.5) id <- c(id, FALSE) 			else id <- c(id, TRUE)
    }
  }

  id <- c(id, TRUE)
  if (ages) {
    xgrid <- xgrid[id]
    label <- label[id]
    if (lastPP$direction == "rightwards")
      graphics::text(xgrid, -0.2, round(label, digits = 1), cex = cex * 0.8)
    if (lastPP$direction == "upwards")
      graphics::text(-0.2, xgrid, round(label, digits = 1), cex = cex * 0.8)

  }

  graphics::par(xpd = FALSE)
}
