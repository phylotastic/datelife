# Jueves 29 de Marzo 2018

# 1. remake runtime tests for make datelife query with tnrs
# use same name number scheme I decided last time.
# screen -r 
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

ninput <- c(10,100,200,300,400,500,600,700,800, 900, 1000, 2000, 3000,4000, 5000, 6000, 7000, 8000, 9000, 10000)
for(i in ninput){
  xname <- paste0("random_sample_",i, "_aves_spp")
  setwd("~/Google Drive/datelife/runtime_tests/1_name_samples")
  load(file=paste0(xname,".RData"))
  y <- microbenchmark(make_datelife_query(input=get(xname)[[1]], use_tnrs = TRUE), times=1L)
  levels(y$expr)[1] <- as.character(i)
  for(j in 2:100){
    yy <- microbenchmark(make_datelife_query(input=get(xname)[[j]], use_tnrs = TRUE), times=1L)
    levels(yy$expr)[1] <- as.character(i)
    y <- rbind(y, yy)
  }	
  rm(list=xname)
  xnameobj <- paste0("make_datelife_query_runtime_tnrs_2018.03.29_", i,"_aves_spp")
  assign(xnameobj, y)
  setwd("~/Google Drive/datelife/runtime_tests/2_tests/2_random_spp_names/2_make_datelife_query")
  save(list=xnameobj, file=paste0(xnameobj,".RData"))
  rm(list=xnameobj)
  print(i)
}

# TO DO
# A. get_datelife_result() with already processed names and current number of chronograms in otol
# There's 202 chronograms in otol
# screen -r 
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

ninput <- c(10,100,200,300,400,500,600,700,800, 900, 1000, 2000, 3000,4000, 5000, 6000, 7000, 8000, 9000, 10000)
for(i in ninput){
  xname <- paste0("random_sample_",i, "_aves_spp")
  setwd("~/Google Drive/datelife/runtime_tests/1_name_samples")
  load(file=paste0(xname,".RData"))
  mdq <- make_datelife_query(input=get(xname)[[1]])
  y <- microbenchmark(get_datelife_result(input = mdq), times=1L)
  levels(y$expr)[1] <- as.character(i)
  for(j in 2:100){
    mdq <- make_datelife_query(input=get(xname)[[i]])
    yy <- microbenchmark(get_datelife_result(input = mdq), times=1L)
    levels(yy$expr)[1] <- as.character(i)
    y <- rbind(y, yy)
  }	
  rm(list=xname)
  xnameobj <- paste0("gfr_runtime_NOmdq_2018.03._", i,"_aves_spp")
  assign(xnameobj, y)
  setwd("~/Google Drive/datelife/runtime_tests/2_tests/2_random_spp_names/1_gfr")
  save(list=xnameobj, file=paste0(xnameobj,".RData"))
  rm(list=xnameobj)
  print(i)
}

# show a figure with different amount of data in chronograms cache
# there are at least 2 bird names in 68 out of 202 chronograms
# that's 68/202 = 33.66% of all chronograms
# chr_coverage <-  c(10, 20, 30, 40, 50)
# also, amount of names in chronogram might be an important parameter

