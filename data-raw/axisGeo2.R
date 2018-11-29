axisGeo2 <- 
function (GTS, tip.time = 0, unit = c("epoch", "period"), ages = TRUE, 
    cex = 1, col = "white", texcol = "black", gridty = 0, gridcol = "black") 
{
    adjustCex <- function(space, string, cex) {
        while (strwidth(string, cex = cex) >= space & cex > 0.001) cex <- cex - 
            0.001
        cex
    }
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    ntips <- lastPP$Ntip
    root <- ntips + 1
    if (lastPP$direction == "rightwards") 
        maxage <- max(lastPP$xx) + tip.time
    if (lastPP$direction == "upwards") 
        maxage <- max(lastPP$yy) + tip.time
    gts <- GTS
    maid <- grep("MA", names(gts))
    gts <- cbind(gts[, 1:maid], c(0, head(gts[, maid], -1)), 
        gts[, ((maid + 1):dim(gts)[2])])
    names(gts)[maid:(maid + 1)] <- c("fromMA", "toMA")
    if (sum(gts[1, maid:(maid + 1)]) == 0) 
        gts <- gts[-1, ]
    ind <- which(gts$fromMA <= maxage & gts$toMA >= tip.time)
    gts <- gts[c(ind, max(ind) + 1), ]
    gts$toMA[1] <- tip.time
    gts$fromMA[dim(gts)[1]] <- maxage
    par(xpd = NA)
    plotGeo <- function(gts, unit, yy) {
        id <- which(names(gts) %in% c(unit, "fromMA", "toMA"))
        gts <- gts[id]
        names(gts) <- c("unit", "from", "to")
        stages <- unique(gts$unit)
        if (length(col) == 1) 
            col <- rep(col, 2)
        col1 <- rep(col, length(stages))
        col1 <- head(col1, length(col1)/2)
        if (length(texcol) == 1) 
            texcol <- rep(texcol, 2)
        col2 <- rep(texcol, length(stages))
        col2 <- head(col2, length(col2)/2)
        xgrid <- NULL
        for (i in seq(along = stages)) {
            cat("\\nStage ", i, ": ", stages[i], sep = "")
            from <- maxage - max(gts[gts$unit == stages[i], 2])
            to <- maxage - min(gts[gts$unit == stages[i], 3])
            rect(from, yy[1], to, yy[2], col = col1[i], border = "black", 
                lwd = 0.5)
            xgrid <- c(xgrid, from, to)
            en <- as.character(stages[i])
            if ((to - from) > strwidth(en, cex = cex)) 
                text(mean(c(from, to)), mean(yy), en, cex = cex, 
                  col = col2[4])
            else {
                thiscex <- adjustCex(to - from, en, cex) * 0.95
                if (2 * thiscex >= cex) 
                  text(mean(c(from, to)), mean(yy), en, cex = thiscex, 
                    col = col2[i])
                else {
                  while (nchar(en) > 0 & strwidth(en, cex = cex) >= 
                    (to - from)) en <- paste(head(unlist(strsplit(en, 
                    "")), -1), collapse = "")
                  if (nchar(en) > 1) 
                    en <- paste(paste(head(unlist(strsplit(en, 
                      "")), -1), collapse = ""), ".", sep = "")
                  text(mean(c(from, to)), mean(yy), en, cex = cex, 
                    col = col2[i])
                }
            }
        }
        xgrid
    }
    plotGeoToLeft <- function(gts, unit, yy) {
        id <- which(names(gts) %in% c(unit, "fromMA", "toMA"))
        gts <- gts[id]
        names(gts) <- c("unit", "from", "to")
        stages <- unique(gts$unit)
        if (length(col) == 1) 
            col <- rep(col, 2)
        col1 <- rep(col, length(stages))
        col1 <- head(col1, length(col1)/2)
        if (length(texcol) == 1) 
            texcol <- rep(texcol, 2)
        col2 <- rep(texcol, length(stages))
        col2 <- head(col2, length(col2)/2)
        xgrid <- NULL
        for (i in seq(along = stages)) {
            from <- maxage - max(gts[gts$unit == stages[i], 2])
            to <- maxage - min(gts[gts$unit == stages[i], 3])
            rect(yy[2], from, yy[1], to, col = col1[i], border = "black", 
                lwd = 0.5)
            xgrid <- c(xgrid, from, to)
            en <- as.character(stages[i])
            yxr <- (max(lastPP$y.lim) - min(lastPP$y.lim))/(max(lastPP$x.lim) - 
                min(lastPP$x.lim)) * 1.5
            if ((to - from) > strwidth(en, cex = cex * yxr)) 
                text(mean(yy), mean(c(from, to)), en, cex = cex, 
                  col = col2[4], srt = 90)
            else {
                asp <- (to - from)/yxr
                thiscex <- adjustCex(asp, en, cex) * 0.95
                if (1.5 * thiscex >= cex) 
                  text(mean(yy), mean(c(from, to)), en, cex = thiscex, 
                    col = col2[i], srt = 90)
                else {
                  while (nchar(en) > 0 & strwidth(en, cex = cex * 
                    yxr) >= (to - from)) en <- paste(head(unlist(strsplit(en, 
                    "")), -1), collapse = "")
                  if (nchar(en) > 1) 
                    en <- paste(paste(head(unlist(strsplit(en, 
                      "")), -1), collapse = ""), ".", sep = "")
                  text(mean(yy), mean(c(from, to)), en, cex = cex, 
                    col = col2[i], srt = 90)
                }
            }
        }
        xgrid
    }
    
    if(ntips == 2){
    	bh <- strheight("Mp", cex = 1) * 1.5
    } else {
	    bh <- -strheight("Mp", cex = cex) * 1.5	
    }
    if (lastPP$direction == "rightwards") {
        if (ages) 
            yy <- c(bh, 2 * bh)
        else yy <- c(0, bh)
        for (j in seq_along(unit)) {
            cat("\\nPlot unit:", unit[j])
            if (j == 1) 
                xgrid <- plotGeo(gts, unit[j], yy)
            else plotGeo(gts, unit[j], yy)
            yy <- yy + bh
        }
    }
    if (lastPP$direction == "upwards") {
        if (ages) 
            yy <- c(bh, 2 * bh)
        else yy <- c(0, bh)
        for (j in seq(along = unit)) {
            if (j == 1) 
                xgrid <- plotGeoToLeft(gts, unit[j], yy)
            else plotGeoToLeft(gts, unit[j], yy)
            yy <- yy + bh
        }
    }
    xgrid <- unique(sort(xgrid, decreasing = TRUE))
    label <- maxage - xgrid
    id <- TRUE
    for (k in seq(along = xgrid)) {
        if (lastPP$direction == "rightwards") 
            lines(rep(xgrid[k], 2), c(0, ntips + 1), lty = gridty, 
                col = gridcol)
        if (lastPP$direction == "upwards") 
            lines(c(0, ntips + 1), rep(xgrid[k], 2), lty = gridty, 
                col = gridcol)
        if (k < length(xgrid)) {
            spneeded <- strwidth(label[k], cex = cex * 0.8)/2
            spavailable <- xgrid[k] - xgrid[k + 1]
            if (spavailable < spneeded * 1.5) 
                id <- c(id, FALSE)
            else id <- c(id, TRUE)
        }
    }
    id <- c(id, TRUE)
    if (ages) {
        xgrid <- xgrid[id]
        label <- label[id]
        if (lastPP$direction == "rightwards") 
            text(xgrid, -0.2, round(label, digits = 1), cex = cex * 
                0.8)
        if (lastPP$direction == "upwards") 
            text(-0.2, xgrid, round(label, digits = 1), cex = cex * 
                0.8)
    }
    par(xpd = FALSE)
}
# <bytecode: 0x7feebcc0bcd8>
# <environment: namespace:phyloch>
