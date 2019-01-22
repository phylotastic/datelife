# Martes 27 de Febrero 2018

# 1. Continue Runtime tests for:
# make_bold_otol_tree(input = c("Rhea americana",  "Struthio camelus", "Gallus gallus"), usetnrs = FALSE, approximatematch = TRUE, marker = "COI", otol_version = "v2", chronogram = TRUE, doML = FALSE, sppfromtaxon=FALSE, verbose=FALSE)
# first with doML=FALSE
[1] 3000
[1] 4000
Error in curl::curl_fetch_memory(x$url$url, handle = x$url$handle) :
  transfer closed with outstanding read data remaining
In addition: There were 50 or more warnings (use warnings() to see the first 50)
# rerun from 5000 names; screen -r 14:24:47 hrs
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
# 2. continue checking for missing documentation entries ... WARNING √ FINISHED √
# Undocumented code objects:
author.pretty
author.results 
curator.pretty
curator.results
missed_doi 
plant_bold_otol_tree 
tb.author.pretty
tb.author.results

str(curator.results)
str(missed_doi)
str(plant_bold_otol_tree)
str(tb.author.pretty)
str(tb.author.results)
levels(names(tb.author.results$person))
names(tb.author.results$person)
factor(names(tb.author.results$person))

# 3. Warning in strptime(xx, f <- "%Y-%m-%d %H:%M:%OS", tz = tz) :
  # unknown timezone 'zone/tz/2018c.1.0/zoneinfo/America/New_York'
?Sys.timezone
Sys.timezone()

OlsonNames() # vector of 594 timezones
Sys.getenv("TZ") # it's empty, we need to set it:
Sys.setenv(TZ = "America/New_York") #typing this in the Rsession 
# where I'm running devtools::check() worked and removed this warning; 
# I'll have to see if the warning appears again in a new session.

# 3b. checking package dependencies ... ERROR
# Package suggested but not available: ‘phylocomr’
# The suggested packages are required for a complete check.
# Checking can be attempted without them by setting the environment
# variable _R_CHECK_FORCE_SUGGESTS_ to a false value.
?check()

# 4. * checking for hidden files and directories ... NOTE
# Found the following hidden files and directories:
#   .tmp_PATHd8.infile
# .tmp_PATHd8.pathd8.orig.out
# .tmp_PATHd8.pathd8.out
# .travis.yml
# These were most likely included in error. See section ‘Package
# structure’ in the ‘Writing R Extensions’ manual.
 
# 5. * checking top-level files ... NOTE
# Non-standard files/directories found at top level:
#   ‘LICENSE.txt’ ‘codecov.yml’ ‘execute’ ‘headers’
# ‘opentree_chronograms.RData’

Checking top-level files. Only specified files and directories are allowed at the top level of the package (e.g. DESCRIPTION, R/, src/). To include other files, you have two choices:
  
  If they don’t need to be installed (i.e. they’re only used in the source package): add them to .Rbuildignore with devtools::use_build_ignore().

If they need to be installed: move them into inst/. They’ll be moved back to the top-level package directory when installed.

# So, I removed execute and opentree_chronograms.RData files, which shouldn't be in the repo and just there respectively
# also, opentree_chronograms.RData was just an old version of it, so I just deleted it

# 6. solve codecoverage travis problem
