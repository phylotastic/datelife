# Martes 2 de Enero 2018
# # Runtime tests: get a good plot
library(ggplot2)
setwd("~/GoogleDrive/datelife/runtime_tests")
ninput <- c(10,100,200,300,400,500,700,1000,1500,2000,3000,5000,7000,8000,9000,10000)
res <- c()
for(i in ninput){
	xname <- paste0("aves",i,".1.gfr.runtime_2017.12.29")
	x <- paste0(xname, ".RData")
	load(x)
	res <- rbind(res, get(xname))
}
length(res)
str(res)
load('~/Google Drive/datelife/runtime_tests/aves.all.gfr.runtime_2017.12.29.RData')
aves.all.gfr.runtime_2017.12.29
res <- rbind(res, aves.all.gfr.runtime_2017.12.29)
autoplot(res)


levels(aves.all.gfr.runtime_2017.12.29$expr) <- "all Aves spp names (12750)"
res <- c()
for(i in ninput){
	xname <- paste0("aves",i,".1.gfr.runtime_2017.12.29")
	res <- rbind(res, get(xname))
}
res <- rbind(res, aves.all.gfr.runtime_2017.12.29)
autoplot(res)
# change time units
# change spp to x axis
str(res)
object <- res
levels(res$expr) <- as.character(c(ninput, "all Aves (12750)"))
n.f <- levels(object$expr)[length(levels(object$expr)):1]
res$expr <- factor(res$expr, levels=n.f)
length(levels(res$time))
microbenchmark:::autoplot.microbenchmark
?convert_to_unit
# next piece of code works, maybe there's a cleaner way, but trying later:
autoplot.gfr <- function (object, ..., log = TRUE, y_max = 1.05 * max(object$time)) {
    y_min <- 0
    object$Time <- microbenchmark:::convert_to_unit(object$time, "t") #changing the name of the element itself is the easiest way to make it appear as axis label
    # object$'Query Length' <- object$expr #changing for a name with spaces won't work...
    plt <- ggplot2::ggplot(object, ggplot2::aes_string(x = "expr", y = "Time"))
    plt <- plt + coord_cartesian(ylim = c(y_min, y_max))
    plt <- plt + stat_ydensity()
    # plt <- plt + xlim(levels(object$expr)[length(levels(object$expr)):1])
    plt <- plt + scale_x_discrete(name = "")
    plt <- plt + theme(axis.text.x = element_text(angle=270))
    plt <- plt + theme(axis.text.y = element_text(angle=315))
    plt <- if (log) {
        # plt + scale_y_log10(name = sprintf("", attr(object$ntime, "unit"))) # this does not work...
        # plt + scale_y_log10(name = sprintf("Time", attr(object$ntime, "unit"))) # this does not work...
        # plt + scale_y_log10(name = "Seconds") # this does not work either...
		plt + scale_y_log10(breaks=c(1e+03, 1e+035, 1e+04, 1e+045, 1e+05), labels=c("1e+03"="1s", "1e+035"="", "1e+04"="10s", "1e+045"="", "1e+05"="100s"), position="top")
    }
    else {
        plt + scale_y_continuous(name = sprintf("Time [%s]", 
            attr(object$ntime, "unit")))
    }
    plt <- plt + ggplot2::coord_flip() # if I inactivate this, I get the following Warning message:
# Transformation introduced infinite values in continuous y-axis 
# Need to figure out how to transform time so I won't get this warning
    plt
}
autoplot.gfr(res) # saved a version of this in ~/Google Drive/datelife/runtime_tests/get_datelife_result_microbenchmark_2018.01.02.pdf
res.plt <- autoplot.gfr(res)
res.plt + scale_y_log10(name = "Seconds") # this does not work
res.plt + scale_y_log10(labels=c("1e+03"="1s", "1e+035"="", "1e+04"="10s", "1e+045"="", "1e+05"="100s")) # this works to change the labels, but I need to specify breaks to be sure it works appropriately
labels(res.plt)
res.plt$labels
str(res.plt)
autoplot.mb(res, log=FALSE)



