# Monday May 7 2018
# Runtime tests plots:
devtools::install_github("phylotastic/datelife", force=TRUE)
library(datelife)
library(ape)
library(geiger)
library(testthat)
library(microbenchmark)
ninput <- c(10,100,200,300,400,500,600,700,800,900, 1000,2000,3000,4000, 5000, 6000,7000,8000,9000,10000)
length(ninput)

# 1. make_datelife_query_tnrs:
setwd("~/Google Drive/datelife/runtime_tests/2_tests/2_random_spp_names/2_make_datelife_query/")
res <- c()
for(i in ninput){
  x <- paste0("make_datelife_query_runtime_tnrs_2018.03.29_",
    i, "_aves_spp")
  xname <- paste0(x, ".RData")
  load(xname)
  res <- rbind(res, get(x))
}
ls()
length(res)
str(res)
str(res$expr)

ggplot2::autoplot(res)
y_min <- 100
y_max <- 15e+05  # 25 minutes
res$Time <- microbenchmark:::convert_to_unit(res$time, "t")  #changing the name of the element itself is the easiest way to make it appear as axis label
plt <- ggplot2::ggplot(res, ggplot2::aes_string(x = "expr", y = "Time"))
plt <- plt + ggplot2::coord_cartesian(ylim = c(y_min, y_max))
plt <- plt + ggplot2::stat_ydensity()
plt <- plt + ggplot2::theme(axis.text.y = ggplot2::element_text(angle=0, size=7))
plt <- plt + ggplot2::scale_y_log10(name="Time (seconds)", breaks=c(1e+03, 1e+04, 3e+04, 6e+04), 
                                     labels=c("1e+03"="1 s", "1e+04"="10 s", "3e+04"="30 s", "6e+04"="60 s"), 
                                     position="left",
                                     sec.axis = ggplot2::sec_axis(~ . *1, name = "Time (minutes)",
                                                                  breaks=c(6e+04, 18e+04, 30e+04, 6e+05, 9e+05, 15e+05), 
                                                                  labels=c("6e+04"="1 min", "18e+04"="3 min", "30e+04"="5 min", "6e+05"="10 min", "9e+05"="15 min", "15e+05"="25 min")
                                     )
)
pdf(file="make_datelife_query_runtime_tnrs_2018.03.29_xtime_all_labels.pdf", width = 6, height = 3)
plt <- plt + ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1))
plt <- plt + ggplot2::scale_x_discrete(name = "Query Length", 
                                       labels=c("10" = "10", 
                                                "100" = "100", 
                                                "200" = "200", 
                                                "300" = "300", 
                                                "400" = "400", 
                                                "500" = "500", 
                                                "600" = "600", 
                                                "700" = "700", 
                                                "800" = "800", 
                                                "900" = "900", 
                                                "1000" = expression(1~0*0*0),
                                                "2000" = expression(2~0*0*0), 
                                                "3000" = expression(3~0*0*0), 
                                                "4000" = expression(4~0*0*0), 
                                                "5000" = expression(5~0*0*0), 
                                                "6000" = expression(6~0*0*0), 
                                                "7000" = expression(7~0*0*0), 
                                                "8000" = expression(8~0*0*0), 
                                                "9000" = expression(9~0*0*0), 
                                                "10000" = expression(10~0*0*0) 
                                       )
)
plt
dev.off()

pdf(file="make_datelife_query_runtime_tnrs_2018.03.29_xtime_some_labels1.pdf", width = 6, height = 3)
plt <- plt + ggplot2::theme(axis.text.x = ggplot2::element_text(angle=0, hjust=0.5, size = 7))
plt <- plt + ggplot2::scale_x_discrete(name = "Query Length", 
                                       labels=c("10" = "10", 
                                                "100" = "100", 
                                                "200" = "", 
                                                "300" = "", 
                                                "400" = "", 
                                                "500" = "500", 
                                                "600" = "", 
                                                "700" = "", 
                                                "800" = "", 
                                                "900" = "", 
                                                "1000" = expression(1~0*0*0),
                                                "2000" = "", 
                                                "3000" = "", 
                                                "4000" = "", 
                                                "5000" = expression(5~0*0*0), 
                                                "6000" = "", 
                                                "7000" = "", 
                                                "8000" = "", 
                                                "9000" = "", 
                                                "10000" = expression(10~0*0*0) 
                                       )
)
plt
dev.off()

# # # # # # # # # # # # # # # # # # # # # #
# 2. make_datelife_query without tnrs:
setwd("~/Google Drive/datelife/runtime_tests/2_tests/2_random_spp_names/2_make_datelife_query/")
res2 <- c()
for(i in ninput){
  x <- paste0("make_datelife_query_runtime_2018.04.07_",
              i, "_aves_spp")
  xname <- paste0(x, ".RData")
  load(xname)
  res2 <- rbind(res2, get(x))
}
ls()
length(res2)
str(res2)
str(res2$expr)

ggplot2::autoplot(res2)
y_min <- 100
y_max <- 15e+05  # 25 minutes
res2$Time <- microbenchmark:::convert_to_unit(res2$time, "t")  #changing the name of the element itself is the easiest way to make it appear as axis label
plt <- ggplot2::ggplot(res2, ggplot2::aes_string(x = "expr", y = "Time"))
plt <- plt + ggplot2::coord_cartesian(ylim = c(y_min, y_max))
plt <- plt + ggplot2::stat_ydensity()
plt <- plt + ggplot2::theme(axis.text.y = ggplot2::element_text(angle=0, size=7))
plt <- plt + ggplot2::scale_y_log10(name="Time (seconds)", breaks=c(1e+03, 1e+04, 3e+04, 6e+04), 
                                    labels=c("1e+03"="1 s", "1e+04"="10 s", "3e+04"="30 s", "6e+04"="60 s"), 
                                    position="left",
                                    sec.axis = ggplot2::sec_axis(~ . *1, name = "Time (minutes)",
                                                                 breaks=c(6e+04, 18e+04, 30e+04, 6e+05, 9e+05, 15e+05), 
                                                                 labels=c("6e+04"="1 min", "18e+04"="3 min", "30e+04"="5 min", "6e+05"="10 min", "9e+05"="15 min", "15e+05"="25 min")
                                    )
)
pdf(file="make_datelife_query_runtime_2018.04.07_xtime_all_labels.pdf", width = 6, height = 3)
plt <- plt + ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1))
plt <- plt + ggplot2::scale_x_discrete(name = "Query Length", 
                                       labels=c("10" = "10", 
                                                "100" = "100", 
                                                "200" = "200", 
                                                "300" = "300", 
                                                "400" = "400", 
                                                "500" = "500", 
                                                "600" = "600", 
                                                "700" = "700", 
                                                "800" = "800", 
                                                "900" = "900", 
                                                "1000" = expression(1~0*0*0),
                                                "2000" = expression(2~0*0*0), 
                                                "3000" = expression(3~0*0*0), 
                                                "4000" = expression(4~0*0*0), 
                                                "5000" = expression(5~0*0*0), 
                                                "6000" = expression(6~0*0*0), 
                                                "7000" = expression(7~0*0*0), 
                                                "8000" = expression(8~0*0*0), 
                                                "9000" = expression(9~0*0*0), 
                                                "10000" = expression(10~0*0*0) 
                                       )
)
plt
dev.off()

pdf(file="make_datelife_query_runtime_2018.04.07_xtime_some_labels1.pdf", width = 6, height = 3)
plt <- plt + ggplot2::theme(axis.text.x = ggplot2::element_text(angle=0, hjust=0.5, size = 7))
plt <- plt + ggplot2::scale_x_discrete(name = "Query Length", 
                                       labels=c("10" = "10", 
                                                "100" = "100", 
                                                "200" = "", 
                                                "300" = "", 
                                                "400" = "", 
                                                "500" = "500", 
                                                "600" = "", 
                                                "700" = "", 
                                                "800" = "", 
                                                "900" = "", 
                                                "1000" = expression(1~0*0*0),
                                                "2000" = "", 
                                                "3000" = "", 
                                                "4000" = "", 
                                                "5000" = expression(5~0*0*0), 
                                                "6000" = "", 
                                                "7000" = "", 
                                                "8000" = "", 
                                                "9000" = "", 
                                                "10000" = expression(10~0*0*0) 
                                       )
)
plt
dev.off()


# get_datelife_result:
# make_bold_otol_tree:
