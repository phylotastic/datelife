# Jueves 11 de Enero 2018
# Runtime tests for:
make_bold_otol_tree
summarizeResults
make_datelife_query √
evaluate cache loading, with different amount of trees?

devtools::install_github("phylotastic/datelife")
library(datelife)
library(ape)
library(geiger)
library(testthat)
library(microbenchmark)

# # # Use Randomized spp names for each run (100 times)
setwd("~/Google Drive/datelife/runtime_tests/name_samples")
for(i in ninput){
	xname <- paste0("random_sample_",i, "_aves_spp")
	load(file=paste0(xname,".RData"))
}
load('~/Google Drive/datelife/runtime_tests/aves.spp.RData')
# # # Launch a first run that is not recorded, so cache is already loaded when evaluating the first set.
asd <- "((Zea mays,Oryza sativa),((Arabidopsis thaliana,(Glycine max,Medicago sativa)),Solanum lycopersicum)Pentapetalae);"
asd.pi <- make_datelife_query(asd)
asd.names <- asd.pi$cleaned.names
make_bold_otol_tree(input=asd.names) 
get_datelife_result(asd, process_input=TRUE)
# now start the test runs:
# screen -r 6417
ninput <- c(10,100,200,300,400,500,700,1000,1500,2000,3000,4000, 5000, 6000,7000,8000,9000,10000)
for(i in ninput){
	xname <- paste0("random_sample_",i, "_aves_spp")
	setwd("~/Google Drive/datelife/runtime_tests/1_name_samples")
	load(file=paste0(xname,".RData"))
	y <- microbenchmark(make_datelife_query(input=get(xname)[[1]]), times=1L)
	levels(y$expr)[1] <- as.character(i)
	for(j in 2:100){
		yy <- microbenchmark(make_datelife_query(input=get(xname)[[j]]), times=1L)
		levels(yy$expr)[1] <- as.character(i)
		y <- rbind(y, yy)
	}	
	rm(list=xname)
	xnameobj <- paste0("make_datelife_query_runtime_2018.01.11_", i,"_aves_spp")
	assign(xnameobj, y)
	setwd("~/Google Drive/datelife/runtime_tests/2_tests/2_random_spp_names/2_make_datelife_query")
	save(list=xnameobj, file=paste0(xnameobj,".RData"))
	rm(list=xnameobj)
}

load('~/Google Drive/datelife/runtime_tests/aves.spp.RData')
length(aves.spp$cleaned.names)
yy <- microbenchmark(make_datelife_query(input=aves.spp$cleaned.names), times=100L)
levels(yy$expr)[1] <- "12750"
assign("make_datelife_query_runtime_2018.01.11_12750_aves_spp", yy)
setwd("~/Google Drive/datelife/runtime_tests/2_tests/0_all_names")
save(make_datelife_query_runtime_2018.01.11_12750_aves_spp, file="make_datelife_query_runtime_2018.01.11_12750_aves_spp.Rdata")

# # # make a plot:
setwd("~/Google Drive/datelife/runtime_tests/2_tests/2_random_spp_names/2_make_datelife_query")
ninput <- c(10,100,200,300,400,500,700,1000,1500,2000,3000,4000, 5000, 6000,7000,8000,9000,10000)
length(ninput)
res <- c()
for(i in ninput){
	x <- paste0("processinput_runtime_2018.01.11_",i,"_aves_spp")
	xname <- paste0("make_datelife_query_runtime_2018.01.11_",i,"_aves_spp.RData")
	load(xname)
	res <- rbind(res, get(x))
}
length(res)
str(res)
load('~/Google Drive/datelife/runtime_tests/2_tests/0_all_names/processinput_runtime_2018.01.11_12750_aves_spp.RData')
processinput_runtime_2018.01.11_12750_aves_spp
# levels(make_datelife_query_runtime_2018.01.11_12750_aves_spp$expr) <- "all Aves (12750)"
res <- rbind(res, processinput_runtime_2018.01.11_12750_aves_spp)
ggplot2::autoplot(res)
setwd("~/Google Drive/datelife/runtime_tests/2_tests/2_random_spp_names/2_make_datelife_query")
pdf(file="make_datelife_query_microbenchmark_2018.01.11_xtime_good_labels.pdf", width = 5, height = 3)
y_min <- 200
y_max <- 1e+5
    res$Time <- microbenchmark:::convert_to_unit(res$time, "t") #changing the name of the element itself is the easiest way to make it appear as axis label
    # object$'Query Length' <- object$expr #changing for a name with spaces won't work...
    plt <- ggplot2::ggplot(res, ggplot2::aes_string(x = "expr", y = "Time"))
    plt <- plt + ggplot2::coord_cartesian(ylim = c(y_min, y_max))
    plt <- plt + ggplot2::stat_ydensity()
    # plt <- plt + ggplot2::scale_x_discrete(name = "Query Length (log10)", 
    #                                        labels=c("10" = "1", 
    #                                                 "100" = expression(10^2), 
    #                                                 "200" = expression(2*"x"*10^2), 
    #                                                 "300" = expression(3*"x"*10^2), 
    #                                                 "400" = expression(4*"x"*10^2), 
    #                                                 "500" = expression(5*"x"*10^2), 
    #                                                 "700" = expression(7*"x"*10^2), 
    #                                                 "1000" = expression(10^3),
    #                                                 "1500" = expression(1.5*"x"*10^3), 
    #                                                 "2000" = expression(2*"x"*10^3), 
    #                                                 "3000" = expression(3*"x"*10^3), 
    #                                                 "4000" = expression(4*"x"*10^3), 
    #                                                 "5000" = expression(5*"x"*10^3), 
    #                                                 "6000" = expression(6*"x"*10^3), 
    #                                                 "7000" = expression(7*"x"*10^3), 
    #                                                 "8000" = expression(8*"x"*10^3), 
    #                                                 "9000" = expression(9*"x"*10^3), 
    #                                                 "10000" = expression(10^4), 
    #                                                 "12750" = expression(1.275*"x"*10^4) 
    #                                        )
    # )
    plt <- plt + ggplot2::scale_x_discrete(name = "Query Length", 
                                           labels=c("10" = "10", 
                                                    "100" = "100", 
                                                    "200" = "200", 
                                                    "300" = "300", 
                                                    "400" = "400", 
                                                    "500" = "500", 
                                                    "700" = "700", 
                                                    "1000" = expression(1~0*0*0),
                                                    "1500" = expression(1~500), 
                                                    "2000" = expression(2~0*0*0), 
                                                    "3000" = expression(3~0*0*0), 
                                                    "4000" = expression(4~0*0*0), 
                                                    "5000" = expression(5~0*0*0), 
                                                    "6000" = expression(6~0*0*0), 
                                                    "7000" = expression(7~0*0*0), 
                                                    "8000" = expression(8~0*0*0), 
                                                    "9000" = expression(9~0*0*0), 
                                                    "10000" = expression(10~0*0*0), 
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
