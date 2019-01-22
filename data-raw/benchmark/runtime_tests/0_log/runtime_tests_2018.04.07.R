# Sat, April 7th 2018

# 1. remake runtime tests for make datelife query without tnrs and
# get_datelife_result() with already processed names
# use current number of chronograms in otol: 202 chronograms 
# screen -r 17:18
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
  y <- microbenchmark(mdq <- make_datelife_query(input=get(xname)[[1]], use_tnrs = FALSE), times=1L)
  levels(y$expr)[1] <- as.character(i)
  g <- microbenchmark(get_datelife_result(input = mdq), times=1L)
  levels(g$expr)[1] <- as.character(i)
  for(j in 2:100){
    yy <- microbenchmark(mdq <- make_datelife_query(input=get(xname)[[j]], use_tnrs = FALSE), times=1L)
    levels(yy$expr)[1] <- as.character(i)
    y <- rbind(y, yy)
    gg <- microbenchmark(get_datelife_result(input = mdq), times=1L)
    levels(gg$expr)[1] <- as.character(i)
    g <- rbind(g, gg)
  }	
  rm(list=xname)
  xnameobjY <- paste0("make_datelife_query_runtime_2018.04.07_", i,"_aves_spp")
  assign(xnameobjY, y)
  setwd("~/Google Drive/datelife/runtime_tests/2_tests/2_random_spp_names/2_make_datelife_query")
  save(list=xnameobjY, file=paste0(xnameobjY,".RData"))
  rm(list=xnameobjY)
  xnameobjG <- paste0("gfr_runtime_2018.04.07_", i,"_aves_spp")
  assign(xnameobjG, g)
  setwd("~/Google Drive/datelife/runtime_tests/2_tests/2_random_spp_names/1_gfr")
  save(list=xnameobjG, file=paste0(xnameobjG,".RData"))
  rm(list=xnameobjG)
  print(i)
}

# finished around 18:40

# 2. retake gbot runtime tests
# screen -r 22:44
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

ninput <- c(6000, 7000, 8000, 9000, 10000)

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
#kohana_error h1,
#kohana_error h2 { margin: 0; padding: 1em; font-size: 1em; font-weight: normal; background: #911; color: #fff; }
#kohana_error h1 a,
#kohana_error h2 a { color: #fff; }
#kohana_error h2 { background: #222; }
#kohana_error h3 { margin: 0; padding: 0.4em 0 0; font-size: 1em; font-weight: normal; }
#kohana_error p { margin: 0; padding: 0.2em 0; }
#kohana_error a { color: #1b323b; }
#kohana_error pre { overflow: auto; white-space: pre-wrap; }
#kohana_error table { width: 100%; display: block; margin: 0 0 0.4em; padding: 0; border-collapse: collapse; background: #fff; }
#kohana_error table td { border: solid 1px #ddd; text-align: left; vertical-align: top; padding: 0.4em; }
#kohana_error div.content { padding: 0.4em 1em 1em; overflow: hidden; }
#kohana_error pre.source { margin: 0 0 1em; padding: 0.4em; background: #fff; border: d
# In addition: There were 50 or more warnings (use warnings() to see the first 50)
# i is 6000
# j is 59

xnameobj <- paste0("gbot_runtime_2018.03.26_", i,"_58_aves_spp")
assign(xnameobj, y)
setwd("~/Google Drive/datelife/runtime_tests/2_tests/2_random_spp_names/3_gbot")
save(list=xnameobj, file=paste0(xnameobj,".RData"))

