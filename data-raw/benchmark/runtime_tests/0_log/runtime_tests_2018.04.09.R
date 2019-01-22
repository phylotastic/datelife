# Mon, April 9th 2018

# 1. retake runtime tests for gbot
# screen -r 17:02:54
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

ninput <- c(7000, 8000, 9000, 10000)

for(i in ninput){
  xname <- paste0("random_sample_",i, "_aves_spp")
  setwd("~/Google Drive/datelife/runtime_tests/1_name_samples")
  load(file=paste0(xname,".RData"))
  y <- microbenchmark(make_bold_otol_tree(input=get(xname)[[1]]), times=1L)
  levels(y$expr)[1] <- as.character(i)
  for(j in 2:100){
    yy <- microbenchmark(make_bold_otol_tree(input=get(xname)[[j]]), times=1L)
    levels(yy$expr)[1] <- as.character(i)
    y <- rbind(y, yy)
  }
  rm(list=xname)
  xnameobj <- paste0("gbot_runtime_2018.03.26_", i,"_aves_spp")
  assign(xnameobj, y)
  setwd("~/Google Drive/datelife/runtime_tests/2_tests/2_random_spp_names/3_gbot")
  save(list=xnameobj, file=paste0(xnameobj,".RData"))
  rm(list=xnameobj)
  print(i)
}

# finished 7000
# stopped 8000 at 65, it was taking too long
