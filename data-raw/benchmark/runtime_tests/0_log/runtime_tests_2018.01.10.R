# Miercoles 10 de Enero 2018
# # Runtime tests for get_datelife_result and plot:
devtools::install_github("phylotastic/datelife", force=TRUE)
library(datelife)
library(ape)
library(geiger)
library(testthat)
library(microbenchmark)

# # # Randomize spp names for each run (100 times)
setwd("~/Google Drive/datelife/runtime_tests")
# load(file="seeds.RData")
# length(seeds)  # 100 seeds; don't use them anymore
load(file="aves.spp.RData")
aves.spp$cleaned.names
ninput <- c(10,100,200,300,400,500,700,1000,1500,2000,3000,4000, 5000, 6000,7000,8000,9000,10000)
setwd("~/Google Drive/datelife/runtime_tests/name_samples")
for(i in ninput){
	x <- vector(mode="list")
	for(j in 1:100){
		x <- c(x, list(sample(aves.spp$cleaned.names, i)))
	}
	xname <- paste0("random_sample_",i, "_aves_spp")
	assign(xname, x)
	save(list=xname, file=paste0(xname,".RData"))
}

# # # Launch a first run that is not recorded, so cache is already loaded when evaluating the first set.
asd <- "((Zea mays,Oryza sativa),((Arabidopsis thaliana,(Glycine max,Medicago sativa)),Solanum lycopersicum)Pentapetalae);"
get_datelife_result(input=asd, process_input=TRUE) # it works. It's defnitely slower when run for the first time. Evaluate cache uploading?
# now start the test runs:
# screen -r 2367
ninput <- c(10,100,200,300,400,500,700,1000,1500,2000,3000,4000, 5000, 6000,7000,8000,9000,10000)
for(i in ninput){
	xname <- paste0("random_sample_",i, "_aves_spp")
	setwd("~/Google Drive/datelife/runtime_tests/1_name_samples")
	load(file=paste0(xname,".RData"))
	y <- microbenchmark(get_datelife_result(input=get(xname)[[1]], process_input=TRUE), times=1L)  # input must be processed :)
	levels(y$expr)[1] <- paste0(i, " names")
	for(j in 2:100){
		yy <- microbenchmark(get_datelife_result(input=get(xname)[[j]], process_input=TRUE), times=1L)
		levels(yy$expr)[1] <- paste0(i, " names")
		y <- rbind(y, yy)
	}	
	rm(list=xname)
	xnameobj <- paste0("gfr_runtime_2018.01.10_", i,"_aves_spp")
	assign(xnameobj, y)
	setwd("~/Google Drive/datelife/runtime_tests/2_tests/2_random_spp_names/1_gfr")
	save(list=xnameobj, file=paste0(xnameobj,".RData"))
	rm(list=xnameobj)
}
setwd("~/Google Drive/datelife/runtime_tests")
load(file="aves.spp.RData")
setwd("~/Google Drive/datelife/runtime_tests/2_tests/0_all_names")
yy <- microbenchmark(get_datelife_result(input=aves.spp$cleaned.names), times=100L)
levels(yy$expr)[1] <- "12750"
aves.all.gfr.runtime_2018.01.12 <- yy
save(aves.all.gfr.runtime_2018.01.12, file="aves.all.gfr.runtime_2018.01.12.RData")

# # # make a plot:
setwd("~/Google Drive/datelife/runtime_tests/2_tests/2_random_spp_names/1_gfr")
ninput <- c(10,100,200,300,400,500,700,1000,1500,2000,3000,4000, 5000, 6000,7000,8000,9000,10000)
length(ninput)
res <- c()
for(i in ninput){
	xname <- paste0("gfr_runtime_2018.01.10_",i,"_aves_spp")
	x <- paste0(xname, ".RData")
	load(x)
	res <- rbind(res, get(xname))
}
length(res)
str(res)
load('~/Google Drive/datelife/runtime_tests/2_tests/0_all_names/aves.all.gfr.runtime_2017.12.29.RData')
aves.all.gfr.runtime_2017.12.29
levels(aves.all.gfr.runtime_2017.12.29$expr) <- "all Aves (12750)"
res <- rbind(res, aves.all.gfr.runtime_2017.12.29)
ggplot2::autoplot(res)
autoplot.gfr(res) # saved a version of this in ~/Google Drive/datelife/runtime_tests/gfr_runtime/2_random_spp_names/get_datelife_result_microbenchmark_2018.01.10_ytime.pdf
levels(res$expr) <- as.character(c(ninput, "all Aves (12750)"))
# n.f <- levels(res$expr)[length(levels(res$expr)):1]
# res$expr <- factor(res$expr, levels=n.f)
autoplot.gfr2(res) 
object <- res
res <- object
# I'm doing the plot witout a function, cause for some reason it will always give the following error message if we do not flip ccords:
# Warning message:
# Transformation introduced infinite values in continuous y-axis 
# I got it!! the problem is that the y limit minimum is 0 and when transforming to log it goes to Inf value...
y_min <- 50
y_max <- 1.05 * max(res$Time)
    res$Time <- microbenchmark:::convert_to_unit(res$time, "t") #changing the name of the element itself is the easiest way to make it appear as axis label
    # object$'Query Length' <- object$expr #changing for a name with spaces won't work...
    plt <- ggplot2::ggplot(res, ggplot2::aes_string(x = "expr", y = "Time"))
    plt <- plt + ggplot2::coord_cartesian(ylim = c(y_min, y_max))
    plt <- plt + ggplot2::stat_ydensity()
    plt <- plt + ggplot2::scale_x_discrete(name = "Query Length")
    plt <- plt + ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1))
    plt <- plt + ggplot2::theme(axis.text.y = ggplot2::element_text(angle=0))
	plt <- 	plt + ggplot2::scale_y_log10(name="Time (seconds)", breaks=c(1e+03, 1e+035, 1e+04, 1e+045, 1e+05), labels=c("1e+03"="1 s", "1e+035"="", "1e+04"="10 s", "1e+045"="", "1e+05"="100 s"), position="left")
    plt
# saved as get_datelife_result_microbenchmark_2018.01.10_xtime.pdf and get_datelife_result_microbenchmark_2018.01.10_xtime_good_labels.pdf

# now change the y axis to log. DOESN'T WORK, CAUSE REMOVES CATEGORIES AND CAN'T CONSTRUCT A HISTOGRAM
# # also, Brian and the team think it's not worthy
# res$names <- res$expr
# levels(res$names) <- c(ninput, 12750)
# as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}  # git this from https://stackoverflow.com/questions/3418128/how-to-convert-a-factor-to-an-integer-numeric-without-a-loss-of-information
# res$names <- as.numeric.factor(res$expr)
    # plt <- ggplot2::ggplot(res, ggplot2::aes_string(x = "names", y = "Time"))
    # plt <- plt + ggplot2::scale_x_log10(name = "Query Length")
# plt

# the plot with the latest runtime tests:
load(file="~/Google Drive/datelife/runtime_tests/2_tests/0_all_names/aves.all.gfr.runtime_2018.01.12.RData")
aves.all.gfr.runtime_2018.01.12
str(aves.all.gfr.runtime_2018.01.12)
res <- rbind(res, aves.all.gfr.runtime_2018.01.12)
setwd("~/Google Drive/datelife/runtime_tests/2_tests/2_random_spp_names/1_gfr")
pdf(file="get_datelife_result_microbenchmark_2018.01.12_xtime_good_labels.pdf", width = 5, height = 3)
y_min <- 200
y_max <- 1e+5
res$Time <- microbenchmark:::convert_to_unit(res$time, "t") #changing the name of the element itself is the easiest way to make it appear as axis label
# object$'Query Length' <- object$expr #changing for a name with spaces won't work...
plt <- ggplot2::ggplot(res, ggplot2::aes_string(x = "expr", y = "Time"))
plt <- plt + ggplot2::coord_cartesian(ylim = c(y_min, y_max))
plt <- plt + ggplot2::stat_ydensity()
# plt <- plt + ggplot2::scale_x_discrete(name = "Query Length (log10)", 
#                                        labels=c("10 names" = "1", 
#                                                 "100 names" = expression(10^2), 
#                                                 "200 names" = expression(2*"x"*10^2), 
#                                                 "300 names" = expression(3*"x"*10^2), 
#                                                 "400 names" = expression(4*"x"*10^2), 
#                                                 "500 names" = expression(5*"x"*10^2), 
#                                                 "700 names" = expression(7*"x"*10^2), 
#                                                 "1000 names" = expression(10^3),
#                                                 "1500 names" = expression(1.5*"x"*10^3), 
#                                                 "2000 names" = expression(2*"x"*10^3), 
#                                                 "3000 names" = expression(3*"x"*10^3), 
#                                                 "4000 names" = expression(4*"x"*10^3), 
#                                                 "5000 names" = expression(5*"x"*10^3), 
#                                                 "6000 names" = expression(6*"x"*10^3), 
#                                                 "7000 names" = expression(7*"x"*10^3), 
#                                                 "8000 names" = expression(8*"x"*10^3), 
#                                                 "9000 names" = expression(9*"x"*10^3), 
#                                                 "10000 names" = expression(10^4), 
#                                                 "12750" = expression(1.275*"x"*10^4) 
#                                        )
# )
plt <- plt + ggplot2::scale_x_discrete(name = "Query Length", 
                                       labels=c("10 names" = "10", 
                                                "100 names" = "100", 
                                                "200 names" = "200", 
                                                "300 names" = "300", 
                                                "400 names" = "400", 
                                                "500 names" = "500", 
                                                "700 names" = "700", 
                                                "1000 names" = expression(1~0*0*0),
                                                "1500 names" = expression(1~500), 
                                                "2000 names" = expression(2~0*0*0), 
                                                "3000 names" = expression(3~0*0*0), 
                                                "4000 names" = expression(4~0*0*0), 
                                                "5000 names" = expression(5~0*0*0), 
                                                "6000 names" = expression(6~0*0*0), 
                                                "7000 names" = expression(7~0*0*0), 
                                                "8000 names" = expression(8~0*0*0), 
                                                "9000 names" = expression(9~0*0*0), 
                                                "10000 names" = expression(10~0*0*0), 
                                                "12750" = expression(12~750) 
                                       )
)
plt <- plt + ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1))
plt <- plt + ggplot2::theme(axis.text.y = ggplot2::element_text(angle=0))
plt <- 	plt + ggplot2::scale_y_log10(name="Time (seconds)", breaks=c(1e+03, 1e+04, 3e+04, 1e+05), 
                                     labels=c("1e+03"="1 s", "1e+04"="10 s", "3e+04"="30 s", "1e+05"="100 s"), 
                                     position="left",
                                     sec.axis = ggplot2::sec_axis(~ . *1, name="Time (minutes)", 
                                                                  breaks=c(6e+03, 3e+04, 6e+04, 9e+04), 
                                                                  labels=c("6e+03"="0.1 min", "3e+04"="0.5 min", "6e+04"="1 min", "9e+04"="1.5 min")
                                     )
)
plt
dev.off()
