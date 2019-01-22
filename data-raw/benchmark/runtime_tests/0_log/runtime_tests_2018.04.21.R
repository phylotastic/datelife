# Sat, April 21st 2018

# 1. relaunch the 6000 one

# screen -r 18:50:58 got an otol error
# again Apr 27 2018, screen -r 18:50:58

setwd("~/Desktop/datelife")
devtools::load_all()
library(ape)
library(geiger)
library(testthat)
library(microbenchmark)
asd <- "((Zea mays,Oryza sativa),((Arabidopsis thaliana,(Glycine max,Medicago sativa)),Solanum lycopersicum)Pentapetalae);"
asd.pi <- make_datelife_query(asd)
asd.names <- asd.pi$cleaned_names
make_bold_otol_tree(input=asd.names)
get_datelife_result(asd)

ninput <- c(6000)

for(i in ninput){
  xname <- paste0("random_sample_",i, "_aves_spp")
  setwd("~/Google Drive/datelife/runtime_tests/1_name_samples")
  load(file=paste0(xname,".RData"))
  y <- microbenchmark(make_bold_otol_tree(input=get(xname)[[1]]), times=1L)
  levels(y$expr)[1] <- as.character(i)
  for(j in 60:100){
    yy <- microbenchmark(make_bold_otol_tree(input=get(xname)[[j]]), times=1L)
    levels(yy$expr)[1] <- as.character(i)
    y <- rbind(y, yy)
    print(j)
  }
  rm(list=xname)
  xnameobj <- paste0("gbot_runtime_2018.03.26_", i,"_60to100_aves_spp")
  assign(xnameobj, y)
  setwd("~/Google Drive/datelife/runtime_tests/2_tests/2_random_spp_names/3_gbot")
  save(list=xnameobj, file=paste0(xnameobj,".RData"))
  rm(list=xnameobj)
  print(i)
}

# finished ok, Apr 28 00:05
# but I need to rename the object, because it starts at 59 and not 60, duu

load('~/Google Drive/datelife/runtime_tests/2_tests/2_random_spp_names/3_gbot/gbot_runtime_2018.03.26_6000_60to100_aves_spp.RData')
ls()
str(gbot_runtime_2018.03.26_6000_60to100_aves_spp)
gbot_runtime_2018.03.26_6000_60to100_aves_spp$expr
gbot_runtime_2018.03.26_6000_60to100_aves_spp$time
gbot_runtime_2018.03.26_6000_59to100_aves_spp <- gbot_runtime_2018.03.26_6000_60to100_aves_spp
save(gbot_runtime_2018.03.26_6000_59to100_aves_spp, file = "~/Google Drive/datelife/runtime_tests/2_tests/2_random_spp_names/3_gbot/gbot_runtime_2018.03.26_6000_59to100_aves_spp.RData")

load('~/Google Drive/datelife/runtime_tests/2_tests/2_random_spp_names/3_gbot/gbot_runtime_2018.03.26_6000_58_aves_spp.RData')
gbot_runtime_2018.03.26_6000_58_aves_spp$time
length(gbot_runtime_2018.03.26_6000_58_aves_spp$time)  # 58, ok

# TO DO
# A. retake the one I stopped at i = 8000
# j = 10


# B. show a figure with different amount of data in chronograms cache
# there are at least 2 bird names in 68 out of 202 chronograms
# that's 68/202 = 33.66% of all chronograms
# chr_coverage <-  c(10, 20, 30, 40, 50)
# also, amount of names in chronogram might be an important parameter
# BRIAN CONSIDERS THIS UNUSEFUL FOR NOW. He suggests WE CAN THINK OF BETTER SEARCHING ALGORITHMS WHEN OTOL DATABASE GROWS
