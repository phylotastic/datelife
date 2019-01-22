autoplot.gfr <- function (object, ..., log = TRUE, y_max = 1.05 * max(object$time)) {
    y_min <- 0
    object$Time <- microbenchmark:::convert_to_unit(object$time, "t") #changing the name of the element itself is the easiest way to make it appear as axis label
    # object$'Query Length' <- object$expr #changing for a name with spaces won't work...
    plt <- ggplot2::ggplot(object, ggplot2::aes_string(x = "expr", y = "Time"))
    plt <- plt + ggplot2::coord_cartesian(ylim = c(y_min, y_max))
    plt <- plt + ggplot2::stat_ydensity()
    # plt <- plt + xlim(levels(object$expr)[length(levels(object$expr)):1])
    plt <- plt + ggplot2::scale_x_discrete(name = "")
    plt <- plt + ggplot2::theme(axis.text.x = ggplot2::element_text(angle=270))
    plt <- plt + ggplot2::theme(axis.text.y = ggplot2::element_text(angle=315))
    plt <- if (log) {
        # plt + scale_y_log10(name = sprintf("", attr(object$ntime, "unit"))) # this does not work...
        # plt + scale_y_log10(name = sprintf("Time", attr(object$ntime, "unit"))) # this does not work...
        # plt + scale_y_log10(name = "Seconds") # this does not work either...
		plt + ggplot2::scale_y_log10(breaks=c(1e+03, 1e+035, 1e+04, 1e+045, 1e+05), labels=c("1e+03"="1s", "1e+035"="", "1e+04"="10s", "1e+045"="", "1e+05"="100s"), position="top")
    }
    else {
        plt + ggplot2::scale_y_continuous(name = sprintf("Time [%s]", 
            attr(object$ntime, "unit")))
    }
    plt <- plt + ggplot2::coord_flip() # if I inactivate this, I get the following Warning message:
# Transformation introduced infinite values in continuous y-axis 
# Need to figure out how to transform time so I won't get this warning
    plt
}
autoplot.gfr2 <- function (object, ..., log = TRUE, y_max = 1.05 * max(object$time)) {
    y_min <- 0
    object$Time <- microbenchmark:::convert_to_unit(object$time, "t") #changing the name of the element itself is the easiest way to make it appear as axis label
    # object$'Query Length' <- object$expr #changing for a name with spaces won't work...
    plt <- ggplot2::ggplot(object, ggplot2::aes_string(x = "expr", y = "Time"))
    plt <- plt + ggplot2::coord_cartesian(ylim = c(y_min, y_max))
    plt <- plt + ggplot2::stat_ydensity()
    # plt <- plt + xlim(levels(object$expr)[length(levels(object$expr)):1])
    plt <- plt + ggplot2::scale_x_discrete(name = "Query Length")
    plt <- plt + ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45))
    plt <- plt + ggplot2::theme(axis.text.y = ggplot2::element_text(angle=0))

	plt <- 	plt + ggplot2::scale_y_log10(breaks=c(1e+03, 1e+035, 1e+04, 1e+045, 1e+05), labels=c("1e+03"="1s", "1e+035"="", "1e+04"="10s", "1e+045"="", "1e+05"="100s"), position="left")

    plt
}
