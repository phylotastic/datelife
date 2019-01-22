# Viernes 16 de Febrero 2018
# Continue Runtime tests for:
# make_bold_otol_tree(input = c("Rhea americana",  "Struthio camelus", "Gallus gallus"), usetnrs = FALSE, approximatematch = TRUE, marker = "COI", otol_version = "v2", chronogram = TRUE, doML = FALSE, sppfromtaxon=FALSE, verbose=FALSE)
# # first with doML=FALSE

devtools::install_github("phylotastic/datelife")
library(datelife)
devtools::load_all()
library(ape)
library(geiger)
library(testthat)
library(microbenchmark)

# Error in curl::curl_fetch_memory(url, handle = handle) :
# Operation was aborted by an application callback

# # # Launch a first run that is not recorded, so cache is already loaded when evaluating the first set.
asd <- "((Zea mays,Oryza sativa),((Arabidopsis thaliana,(Glycine max,Medicago sativa)),Solanum lycopersicum)Pentapetalae);"
asd.pi <- make_datelife_query(asd)
asd.names <- asd.pi$cleaned.names
make_bold_otol_tree(input=asd.names)
get_datelife_result(asd)
i <- 300
xname <- paste0("random_sample_",i, "_aves_spp")
setwd("~/Google Drive/datelife/runtime_tests/1_name_samples")
load(file=paste0(xname,".RData"))
make_bold_otol_tree(input=random_sample_300_aves_spp[[1]])
# Warning messages:
# 1: In check_tnrs(res) : Chrysococcyx rufomerus are not matched
# 2: In rotl::tnrs_match_names(names = input[xx[i]:yy[i]]) :
#   NAs introduced by coercion
# 3: In make_bold_otol_tree(input = random_sample_300_aves_spp[[1]]) :
#   	Negative branch lengths in BOLD chronogram.
rotl::tnrs_match_names(names = "Chrysococcyx rufomerus")
# restart the test from 200: screen -r 16:47 hrs
ninput <- c(200,300,400,500,700,1000,1500,2000,3000,4000, 5000, 6000,7000,8000,9000,10000)
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
	xnameobj <- paste0("gbot_runtime_2018.02.09_", i,"_aves_spp")
	assign(xnameobj, y)
	setwd("~/Google Drive/datelife/runtime_tests/2_tests/2_random_spp_names/3_gbot")
	save(list=xnameobj, file=paste0(xnameobj,".RData"))
	rm(list=xnameobj)
	print(i)
}
