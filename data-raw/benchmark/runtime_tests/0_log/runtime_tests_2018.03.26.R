# Lunes 26 de Marzo 2018

# 1. Continue Runtime tests for make_bold_otol_tree
# I was having a lot of problems while running this tests: internet connections and then resposiveness of otol
# so I'm repeating every test since the beggining
# I'm also doing a different name sampling strategy. 10, 100:1000 by 100; 1000:10 000 by 1000; 10 000: 1 million by power 10
setwd("~/Google Drive/datelife/runtime_tests")
load(file="aves.spp.RData")
aves.spp$cleaned.names
ninput <- c(600, 800, 900)
setwd("~/Google Drive/datelife/runtime_tests/1_name_samples")
for(i in ninput){
  x <- vector(mode="list")
  for(j in 1:100){
    x <- c(x, list(sample(aves.spp$cleaned.names, i)))
  }
  xname <- paste0("random_sample_",i, "_aves_spp")
  assign(xname, x)
  save(list=xname, file=paste0(xname,".RData"))
}

# screen -r Martes 27, 6:15 am
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
# stopped it at i =5000
# j = 10
xnameobj <- paste0("gbot_runtime_2018.03.26_", i,"9_aves_spp")
assign(xnameobj, y)
setwd("~/Google Drive/datelife/runtime_tests/2_tests/2_random_spp_names/3_gbot")
save(list=xnameobj, file=paste0(xnameobj,".RData"))
rm(list=xnameobj)
print(i)

