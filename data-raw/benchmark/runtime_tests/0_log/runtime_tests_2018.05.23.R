# Wednesday May 23, Thursday 24, Tuesday 29 2018
# Runtime tests plots:
devtools::install_github("phylotastic/datelife", force=TRUE)
library(datelife)
library(ape)
library(geiger)
library(testthat)
library(microbenchmark)
ninput <- c(10,100,200,300,400,500,600,700,800,900, 1000,2000,3000,4000, 5000, 6000,7000,8000,9000,10000)
length(ninput)

setwd("~/Google Drive/datelife/runtime_tests/2_tests/2_random_spp_names/2_make_datelife_query/")
# 1. make_datelife_query_tnrs:
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

# 2. make_datelife_query:
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
length(res2$expr)
length(res2$time)
head(res2$time)


ggplot2::autoplot(res)
ggplot2::autoplot(res2)
res$Time <- microbenchmark:::convert_to_unit(res$time, "t")  #changing the name of the element itself is the easiest way to make it appear as axis label
res2$Time <- microbenchmark:::convert_to_unit(res2$time, "t")  #changing the name of the element itself is the easiest way to make it appear as axis label

res$tnrs <- rep("Using TNRS", 2000)
res2$tnrs <- rep("Not using TNRS", 2000)
res_all <- dplyr::bind_rows(res, res2)
names(res_all)
lapply(res_all, length)
lapply(res_all, is.factor)
res_all$tnrs <- as.factor(res_all$tnrs)
y_min <- 100
y_max <- 15e+05  # 25 minutes

rm(plt)
plt
pdf(file="make_datelife_query_runtime_BOTH_Legend_xtime_all_labels_2018.05.23_scale_COUNT_identity_gray_bkgd.pdf", width = 6, height = 3)
plt <- ggplot2::ggplot(res_all, ggplot2::aes_string(x = "expr", y = "Time"))
plt <- plt + ggplot2::coord_cartesian(ylim = c(y_min, y_max))
# plt <- plt + ggplot2::stat_ydensity()
# # this plots both at the same time:
plt <- plt + ggplot2::stat_ydensity(ggplot2::aes(color = tnrs, fill = tnrs), scale = "count", position = "identity")  # this and follow line work the same
# color controls the outer line color; fill control the inside color.
# plt <- plt + ggplot2::geom_violin(ggplot2::aes(fill = tnrs))  # fill = factor(tnrs) if tnrs is not a factor
plt <- plt + ggplot2::scale_color_manual(values = c("black", "grey")) +
          ggplot2::scale_fill_manual(values = c("black", "grey"))
plt <- plt + ggplot2::theme(axis.text.y = ggplot2::element_text(angle=0, size=7), legend.position = "none")
plt <- plt + ggplot2::scale_y_log10(name="Time (seconds)", breaks=c(1e+03, 1e+04, 3e+04, 6e+04), 
                                     labels=c("1e+03"="1 s", "1e+04"="10 s", "3e+04"="30 s", "6e+04"="60 s"), 
                                     position="left",
                                     sec.axis = ggplot2::sec_axis(~ . *1, name = "Time (minutes)",
                                                                  breaks=c(6e+04, 18e+04, 30e+04, 6e+05, 9e+05, 15e+05), 
                                                                  labels=c("6e+04"="1 min", "18e+04"="3 min", "30e+04"="5 min", "6e+05"="10 min", "9e+05"="15 min", "15e+05"="25 min")
                                     )
)
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
plt <- plt + ggplot2::theme(legend.position = c(0.07,0.87), 
                            legend.title = ggplot2::element_text(size = 0, face = "bold"), # size = 0 eliminate legend title
                            legend.text = ggplot2::element_text(size = 5),
                            legend.background= ggplot2::element_rect(fill= "white", size = 5),  # legend.key does not remove all legend background, only key background
                            legend.key.width= ggplot2::unit(0.4,"line"),
                            legend.key.height= ggplot2::unit(0.3,"line") 
                            # legend.key.size = ggplot2::unit(1,"line")
                            )
# plt <- plt + ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=1)))
# plt <- plt + ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white"),
#                             panel.grid.major.y = ggplot2::element_line(color = "grey", linetype = "solid"),
#                             panel.grid.minor.y = ggplot2::element_line(color = "grey", linetype = "solid"),
#                             panel.grid.major.x = ggplot2::element_blank(),
#                             panel.grid.minor.x = ggplot2::element_blank())
plt
dev.off()