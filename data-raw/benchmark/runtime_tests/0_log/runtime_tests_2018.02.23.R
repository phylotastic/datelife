# Viernes 23 de Febrero 2018

# 1. Continue Runtime tests for:
# make_bold_otol_tree(input = c("Rhea americana",  "Struthio camelus", "Gallus gallus"), usetnrs = FALSE, approximatematch = TRUE, marker = "COI", otol_version = "v2", chronogram = TRUE, doML = FALSE, sppfromtaxon=FALSE, verbose=FALSE)
# first with doML=FALSE
# finished 2000, but internet was disconnected
Error in curl::curl_fetch_memory(x$url$url, handle = x$url$handle) :
  transfer closed with outstanding read data remaining
In addition: There were 50 or more warnings (use warnings() to see the first 50)
# rerun from 3000 names; screen -r 14:45 hrs
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

ninput <- c(3000,4000, 5000, 6000,7000,8000,9000,10000)
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
	xnameobj <- paste0("gbot_runtime_2018.02.22_", i,"_aves_spp")
	assign(xnameobj, y)
	setwd("~/Google Drive/datelife/runtime_tests/2_tests/2_random_spp_names/3_gbot")
	save(list=xnameobj, file=paste0(xnameobj,".RData"))
	rm(list=xnameobj)
	print(i)
}
# 2. continue with this: checking for missing documentation entries ... WARNING
# Undocumented code objects:
#   ‘author.pretty’ ‘author.results’ ‘curator.pretty’ ‘curator.results’
# ‘missed_doi’ ‘plant_bold_otol_tree’ ‘tb.author.pretty’
# ‘tb.author.results’
# Undocumented data sets:
#   ‘author.pretty’ ‘author.results’ ‘curator.pretty’ ‘curator.results’
# ‘missed_doi’ ‘plant_bold_otol_tree’ ‘tb.author.pretty’
# ‘tb.author.results’
# All user-level objects in a package should have documentation entries.
# See chapter ‘Writing R documentation files’ in the ‘Writing R
# Extensions’ manual.
plant_bold_otol_tree

# 5. solve codecoverage travis problem
