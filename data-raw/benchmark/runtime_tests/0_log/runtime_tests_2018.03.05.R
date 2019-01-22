# Lunes 5 de Marzo 2018

# 1. Continue Runtime tests for:
# make_bold_otol_tree(input = c("Rhea americana",  "Struthio camelus", "Gallus gallus"), usetnrs = FALSE, approximatematch = TRUE, marker = "COI", otol_version = "v2", chronogram = TRUE, doML = FALSE, sppfromtaxon=FALSE, verbose=FALSE)
# first with doML=FALSE
# [1] 5000
# [1] 6000
# Error: HTTP failure: 502
# <!DOCTYPE HTML PUBLIC "-//IETF//DTD HTML 2.0//EN">
#   <html><head>
#   <title>502 Proxy Error</title>
#   </head><body>
#   <h1>Proxy Error</h1>
#   <p>The proxy server received an invalid
# response from an upstream server.<br />
#   The proxy server could not handle the request <em><a href="/v3/tnrs/match_names">POST&nbsp;/v3/tnrs/match_names</a></em>.<p>
#   Reason: <strong>Error reading from remote server</strong></p></p>
#   <hr>
#   <address>Apache/2.4.10 (Debian) Server at api.opentreeoflife.org Port 443</address>
#   </body></html>
#   In addition: There were 50 or more warnings (use warnings() to see the first 50)
# rerun from 7000 names; screen -r 14:24:47 hrs
# this time, Brian suggested I connected through ethernet, so I stop getting curl errors
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
# se quedó en j = 55 de i = 7000

	for(j in 2:100){
		yy <- microbenchmark(make_bold_otol_tree(input=get(xname)[[j]]), times=1L)
		levels(yy$expr)[1] <- as.character(i)
		y <- rbind(y, yy)
	}
	xnameobj <- paste0("gbot_runtime_2018.02.22_", i,"54_aves_spp")
	assign(xnameobj, y)
	setwd("~/Google Drive/datelife/runtime_tests/2_tests/2_random_spp_names/3_gbot")
	save(list=xnameobj, file=paste0(xnameobj,".RData"))
	rm(list=xnameobj)
 # saved the results and quit R
	
	# open new screen:
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
	
	load("~/Google Drive/datelife/runtime_tests/2_tests/2_random_spp_names/3_gbot/gbot_runtime_2018.02.22_7000_54_aves_spp.RData")
y <- gbot_runtime_2018.02.22_7000_54_aves_spp
i <- 7000
xname <- paste0("random_sample_",i, "_aves_spp")
setwd("~/Google Drive/datelife/runtime_tests/1_name_samples")
load(file=paste0(xname,".RData"))
for(j in 55:100){
	  yy <- microbenchmark(make_bold_otol_tree(input=get(xname)[[j]]), times=1L)
	  levels(yy$expr)[1] <- as.character(i)
	  y <- rbind(y, yy)
	}
	xnameobj <- paste0("gbot_runtime_2018.02.22_", i,"aves_spp")
	assign(xnameobj, y)
	setwd("~/Google Drive/datelife/runtime_tests/2_tests/2_random_spp_names/3_gbot")
	save(list=xnameobj, file=paste0(xnameobj,".RData"))
	
	
	
	ninput <- c(5000, 6000,7000,8000,9000,10000)
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


