# Jueves 22 de Febrero 2018

# 1. Continue Runtime tests for:
# make_bold_otol_tree(input = c("Rhea americana",  "Struthio camelus", "Gallus gallus"), usetnrs = FALSE, approximatematch = TRUE, marker = "COI", otol_version = "v2", chronogram = TRUE, doML = FALSE, sppfromtaxon=FALSE, verbose=FALSE)
# first with doML=FALSE
# rerun from 2000 names; screen -r 16:32 hrs
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

ninput <- c(2000,3000,4000, 5000, 6000,7000,8000,9000,10000)
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
# 2. fix changed warnings to errors travis CI build problem √
# checked with Brian. It's not a problem. It's supposed to fail cause now warnings are errors. The CI buils will pass when there are no more warnings

# 3. fix other rcmd check problem: datelife_result_sdm.Rd: non-ASCII input and no declared encoding √
  # problem found in ‘datelife_result_sdm.Rd’
# we need to declare encoding in Desciprion file. Just added the line: 'Encoding: latin1' in the DESCRIPTION file and it worked :)
??devtools::encoding
library(devtools)
?devtools
??encoding

# 4. fix this: checking for missing documentation entries ... WARNING
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
# # so, we need to make .Rd files for all objects
# # maybe there's a function to do them in roxygen2
??roxygen2
# found it! you have to make a description block as for a function and the roxygenise() fn will write the .Rd files for you
# from https://stackoverflow.com/questions/2310409/how-can-i-document-data-sets-with-roxygen/22598293
#' This is data to be included in my package
#'
#' @name data-name
#' @docType data
#' @author My Name \email{blahblah@@roxygen.org}
#' @references \url{data_blah.com}
#' @keywords data
NULL

# from https://cran.r-project.org/web/packages/roxygen2/vignettes/rd.html
#' Prices of 50,000 round cut diamonds.
#'
#' A dataset containing the prices and other attributes of almost 54,000
#' diamonds.
#'
#' @format A data frame with 53940 rows and 10 variables:
#' \describe{
#'   \item{price}{price, in US dollars}
#'   \item{carat}{weight of the diamond, in carats}
#'   ...
#' }
#' @source \url{http://www.diamondse.info/}
"diamonds"

# 5. solve codecoverage travis problem
