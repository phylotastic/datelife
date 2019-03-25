# hosting phyloch functions so cran won't complain
#' Function to add interval bars to nodes on a phylo plot.
#' @param phy A phylo object with two extra elements specifying the MAX and MIN limits of intervals to plot.
#' @param label A character vector corresponding to the name used to label intervals stored in phy
#' @param tab Something I'm not sure of, should ask Christoph Heibl.
#' @param nodes Numeric. Use this if you wanna plot bars on certain nodes only.
#' @param col Color for the bars, default to "grey75"
#' @param lwd Line width for the bars, default to 5.
#' @param broken Something I'm not sure of, should ask Christoph Heibl.
#' @inheritDotParams graphics::segments
#' @return Adds bars to nodes on a phylo plot.
#' @export
#' @details Make sure ape is loaded otherwise it won't find .PlotPhyloEnv
HPDbars <-
function(phy, label = "height_95%_HPD", tab = NULL, nodes, col, lwd, broken = FALSE, ...){

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
