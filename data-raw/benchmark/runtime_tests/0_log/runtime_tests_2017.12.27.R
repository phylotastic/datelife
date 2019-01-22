# Miercoles 27 a Viernes 28 de Diciembre 2017
# testing fns run time
install.packages("pbapply") 
library(pbapply)
?pbapply #adding progress bar to apply functions

# # Rprof() function

setwd("~/Google Drive/datelife/runtime_tests")
Rprof(filename = "make_bold_otol_tree_Rprof_1", append = FALSE, interval = 0.02,
       memory.profiling = TRUE, gc.profiling = FALSE, 
       line.profiling = TRUE, numfiles = 100L, bufsize = 10000L)
asd <- "((Zea mays,Oryza sativa),((Arabidopsis thaliana,(Glycine max,Medicago sativa)),Solanum lycopersicum)Pentapetalae);"
asd.bt <- make_bold_otol_tree(asd)
Rprof(NULL)

# # running profiles graphically
install.packages("GUIProfiler")
# # # install GUIProfiler depedencies not in cran:
source("https://bioconductor.org/biocLite.R")
biocLite("Rgraphviz")
# # # # retry installation of packages that had non-sero exit status 
install.packages("RCurl")
install.packages("rmetasim")
biocLite("graph")
library(GUIProfiler)
?GUIProfiler
x <- summaryRprof(filename = "make_bold_otol_tree_Rprof_1")
ls(x)
# $by.self
                 # self.time self.pct total.time total.pct
# "rawToChar"           0.50    43.10       0.50     43.10
# ".Call"               0.20    17.24       0.20     17.24
# "system"              0.12    10.34       0.12     10.34
# "paste0"              0.06     5.17       0.56     48.28
# "newick"              0.04     3.45       0.04      3.45
# "file.exists"         0.02     1.72       0.02      1.72
# "gsub"                0.02     1.72       0.02      1.72
# "integer"             0.02     1.72       0.02      1.72
# "match.fun"           0.02     1.72       0.02      1.72
# "readRDS"             0.02     1.72       0.02      1.72
# "regexec"             0.02     1.72       0.02      1.72
# "rep"                 0.02     1.72       0.02      1.72
# "scan"                0.02     1.72       0.02      1.72
# "strptime"            0.02     1.72       0.02      1.72
# "strsplit"            0.02     1.72       0.02      1.72
# "textConnection"      0.02     1.72       0.02      1.72
# "unlist"              0.02     1.72       0.02      1.72

# $by.total
                             # total.time total.pct self.time self.pct
# "make_bold_otol_tree"                  1.16    100.00      0.00     0.00
# "bold::bold_seqspec"               0.76     65.52      0.00     0.00
# "paste0"                           0.56     48.28      0.06     5.17
# "rawToChar"                        0.50     43.10      0.50    43.10
# ".Call"                            0.20     17.24      0.20    17.24
# "curl::curl_fetch_memory"          0.20     17.24      0.00     0.00
# "ape::multi2di"                    0.16     13.79      0.00     0.00
# "phangorn::as.phyDat"              0.16     13.79      0.00     0.00
# "rotl::tol_induced_subtree"        0.16     13.79      0.00     0.00
# ".tnrs_match_names"                0.14     12.07      0.00     0.00
# ".tol_induced_subtree"             0.14     12.07      0.00     0.00
# "b_GET"                            0.14     12.07      0.00     0.00
# "cli$get"                          0.14     12.07      0.00     0.00
# "ips::mafft"                       0.14     12.07      0.00     0.00
# "otl_POST"                         0.14     12.07      0.00     0.00
# "rotl::tnrs_match_names"           0.14     12.07      0.00     0.00
# "system"                           0.12     10.34      0.12    10.34
# "httr::POST"                       0.12     10.34      0.00     0.00
# "private$make_request"             0.12     10.34      0.00     0.00
# "request_perform"                  0.12     10.34      0.00     0.00
# "crul_fetch"                       0.10      8.62      0.00     0.00
# "request_fetch.write_memory"       0.10      8.62      0.00     0.00
# "request_fetch"                    0.10      8.62      0.00     0.00
# "read.table"                       0.06      5.17      0.00     0.00
# "utils::read.delim"                0.06      5.17      0.00     0.00
# "newick"                           0.04      3.45      0.04     3.45
# "ape::collapse.singles"            0.04      3.45      0.00     0.00
# "lapply"                           0.04      3.45      0.00     0.00
# "phytools::read.newick"            0.04      3.45      0.00     0.00
# "make_datelife_query"                     0.04      3.45      0.00     0.00
# "ProcessPhy"                       0.04      3.45      0.00     0.00
# "file.exists"                      0.02      1.72      0.02     1.72
# "gsub"                             0.02      1.72      0.02     1.72
# "integer"                          0.02      1.72      0.02     1.72
# "match.fun"                        0.02      1.72      0.02     1.72
# "readRDS"                          0.02      1.72      0.02     1.72
# "regexec"                          0.02      1.72      0.02     1.72
# "rep"                              0.02      1.72      0.02     1.72
# "scan"                             0.02      1.72      0.02     1.72
# "strptime"                         0.02      1.72      0.02     1.72
# "strsplit"                         0.02      1.72      0.02     1.72
# "textConnection"                   0.02      1.72      0.02     1.72
# "unlist"                           0.02      1.72      0.02     1.72
# "%in%"                             0.02      1.72      0.00     0.00
# "<Anonymous>"                      0.02      1.72      0.00     0.00
# "as.phyDat.DNAbin"                 0.02      1.72      0.00     0.00
# "as.POSIXct"                       0.02      1.72      0.00     0.00
# "c_time"                           0.02      1.72      0.00     0.00
# "data.frame"                       0.02      1.72      0.00     0.00
# "do.call"                          0.02      1.72      0.00     0.00
# "force"                            0.02      1.72      0.00     0.00
# "FUN"                              0.02      1.72      0.00     0.00
# "headers_parse"                    0.02      1.72      0.00     0.00
# "HttpResponse$new"                 0.02      1.72      0.00     0.00
# "ifelse"                           0.02      1.72      0.00     0.00
# "internal_make_phylo"              0.02      1.72      0.00     0.00
# "make_ua"                          0.02      1.72      0.00     0.00
# "make.names"                       0.02      1.72      0.00     0.00
# "matrix"                           0.02      1.72      0.00     0.00
# "order"                            0.02      1.72      0.00     0.00
# "packageDescription"               0.02      1.72      0.00     0.00
# "parse_http_date"                  0.02      1.72      0.00     0.00
# "phyDat.DNA"                       0.02      1.72      0.00     0.00
# "phylo_from_otl"                   0.02      1.72      0.00     0.00
# "read.fas"                         0.02      1.72      0.00     0.00
# "rncl::read_newick_phylo"          0.02      1.72      0.00     0.00
# "rncl"                             0.02      1.72      0.00     0.00
# "sapply"                           0.02      1.72      0.00     0.00
# "suppressWarnings"                 0.02      1.72      0.00     0.00
# "tolower"                          0.02      1.72      0.00     0.00
# "utils::packageVersion"            0.02      1.72      0.00     0.00
# "vapply"                           0.02      1.72      0.00     0.00
# "withCallingHandlers"              0.02      1.72      0.00     0.00

# $sample.interval
# [1] 0.02

# $sampling.time
# [1] 1.16

x <- make_datelife_query(input="Thraupidae", sppfromtaxon=TRUE) #532 spp

x <- make_datelife_query(input="Aves", sppfromtaxon=TRUE) # 12750 spp names
ls(x)
length(x[[2]])
x[[1]] #NA 

# # microbenchmark
install.packages("microbenchmark")
library(microbenchmark)
?microbenchmark
f <- function() NULL
res1 <- microbenchmark(NULL, f(), times=1000L)
res1
autoplot(res1)
my_check <- function(values) {
  all(sapply(values[-1], function(x) identical(values[[1]], x)))
}

f2 <- function(a, b)
  2 + 2

a <- 2
res2 <- microbenchmark(2 + 2, 2 + a, f2(2, a), f2(2, 2), check=my_check)
res2
if (require("ggplot2")) {
  autoplot(res2)
}

res <- rbind(res1,res2)
res
autoplot(res)
resx <- microbenchmark(EstimateDates(input="Thraupidae", sppfromtaxon=TRUE), times=1L)
resx
# # # screen -r 27846
# # # # 2017.12.28
devtools::install_github("phylotastic/datelife")
library(datelife)
library(ape)
library(geiger)
library(testthat)
library(microbenchmark)
# time to get all available phylogenies, as phylo objects (output.format="phylo.all"):
res <- microbenchmark(EstimateDates(input="Thraupidae", sppfromtaxon=TRUE), times=100L)
setwd("~/Google Drive/datelife/runtime_tests")
thraupidae.ed.runtime_2017.12.28 <- res
save(thraupidae.ed.runtime_2017.12.28, file="thraupidae.ed.runtime_2017.12.28.RData")

# # # screen -r 28095
# # # # 2017.12.28
devtools::install_github("phylotastic/datelife")
library(datelife)
library(ape)
library(geiger)
library(testthat)
library(microbenchmark)
# time to get all available phylogenies, as phylo objects (output.format="phylo.all"):
aves.ed.runtime_2017.12.28 <- microbenchmark(EstimateDates(input="Aves", sppfromtaxon=TRUE), times=1L)
setwd("~/Google Drive/datelife/runtime_tests")
save(aves.ed.runtime_2017.12.28, file="aves.ed.runtime_2017.12.28.RData")
load('~/Google Drive/datelife/runtime_tests/aves.ed.runtime_2017.12.28.RData')
aves.ed.runtime_2017.12.28
# Unit: seconds
                                               # expr      min       lq     mean
 # EstimateDates(input = "Aves", sppfromtaxon = TRUE) 2813.751 2813.751 2813.751
   # median       uq      max neval
 # 2813.751 2813.751 2813.751     1


load(file="aves.spp.RData")
length(aves.spp$cleaned.names) #12750
x <- get_datelife_result(input=aves.spp$cleaned.names)
aves.gfr_2017.12.28 <- x
save(aves.gfr_2017.12.28, file="aves.gfr_2017.12.28.RData")
# # # # finished ok

# # # play with all aves gfr results
load("/Users/luna/Google Drive/datelife/runtime_tests/aves.gfr_2017.12.28.RData")
length(aves.gfr_2017.12.28) #68 phylogenies
aves_sample <- sapply(aves.gfr_2017.12.28, function(x) dim(x)[1])
plot(aves_sample)
aves_sample[which.max(aves_sample)] # Hedges et al. 2015 has 9562 aves spp
aves_sample[which(aves_sample>5000)]
names(aves_sample)
grep("Jetz", names(aves_sample))
aves_sample[46:47] # Jetz et al. 2012 study has only 6658 aves spp



# # # Check when the following happens:
# # # # spp <- make_datelife_query(input="Thraupidae", sppfromtaxon=TRUE)
# # # # Error in open.connection(con, "rb") : HTTP error 500.
# # # # # it happens if sppfromtaxon=TRUE, at rphylotastic::taxon_get_species

# # # tests to loop:
x <- make_datelife_query(input="Aves", sppfromtaxon=TRUE) # 12750 spp names
aves.spp <- x
aves.spp[[2]]
save("aves.spp", file="aves.spp.RData")
load(file="aves.spp.RData")

set.seed(10)
seeds <- runif(100, 1, 1e9) # set.seed only accepts number up to 9 integers-ish:
# # # # set.seed(2140000000)
# # # # works with numbers <=2.14e+09
seeds
save(seeds, file="seeds.RData")
set.seed(seeds[1])
spp400.1 <- sample(aves.spp$cleaned.names, 400)
aves400.1.gfr.runtime_2017.12.28 <- microbenchmark(get_datelife_result(input=spp400.1), times=100L)
save(aves400.1.gfr.runtime_2017.12.28, file="aves400.1.gfr.runtime_2017.12.28.RData")

load('~/Google Drive/datelife/runtime_tests/thraupidae.ed.runtime_2017.12.28.RData')
load('~/Google Drive/datelife/runtime_tests/aves400.1.gfr.runtime_2017.12.28.RData')
str(aves400.1.gfr.runtime_2017.12.28)
res <- rbind(aves400.1.gfr.runtime_2017.12.28, thraupidae.ed.runtime_2017.12.28)
str(res)
# Classes ‘microbenchmark’ and 'data.frame':	200 obs. of  2 variables:
 # $ expr: Factor w/ 2 levels "get_datelife_result(input = spp400.1)",..: 1 1 1 1 1 1 1 1 1 1 ...
 # $ time: num  1.99e+09 9.45e+08 7.66e+08 6.95e+08 6.38e+08 ...
res$expr
res$time
if (require("ggplot2")) {
  autoplot(res)
}
ls(res)
ninput <- c(10,100,200,300,400,500,700,1000,1500,2000,3000,5000,7000)
for (i in ninput){
	set.seed(seeds[1])
	x <- sample(aves.spp$cleaned.names, i)
	xname <- paste0("spp",i)
	assign(xname, x)
	save(list=xname, file=paste0(xname,".RData"))
}
ninput <- c(8000,9000,10000)
for (i in ninput){
	set.seed(seeds[1])
	x <- sample(aves.spp$cleaned.names, i)
	xname <- paste0("spp",i)
	assign(xname, x)
	save(list=xname, file=paste0(xname,".RData"))
}


# # # screen -r 27846
setwd("~/Google Drive/datelife/runtime_tests")
ninput <- c(10,100,200,300,400,500,700,1000)
for(i in ninput){
	x <- paste0("spp",i,".RData")
	load(x)
}
devtools::install_github("phylotastic/datelife")
library(datelife)
library(ape)
library(geiger)
library(testthat)
library(microbenchmark)
aves10.1.gfr.runtime_2017.12.28 <- microbenchmark(get_datelife_result(input=spp10), times=100L)
aves100.1.gfr.runtime_2017.12.28 <- microbenchmark(get_datelife_result(input=spp100), times=100L)
aves200.1.gfr.runtime_2017.12.28 <- microbenchmark(get_datelife_result(input=spp200), times=100L)
aves300.1.gfr.runtime_2017.12.28 <- microbenchmark(get_datelife_result(input=spp300), times=100L)
aves400.1.gfr.runtime_2017.12.28 <- microbenchmark(get_datelife_result(input=spp400), times=100L)
aves500.1.gfr.runtime_2017.12.28 <- microbenchmark(get_datelife_result(input=spp500), times=100L)
aves700.1.gfr.runtime_2017.12.28 <- microbenchmark(get_datelife_result(input=spp700), times=100L)
aves1000.1.gfr.runtime_2017.12.28 <- microbenchmark(get_datelife_result(input=spp1000), times=100L)
for(i in ninput){
	xname <- paste0("aves",i,".1.gfr.runtime_2017.12.28")
	save(list=xname, file=paste0(xname, ".RData"))
}
# # # # finished ok

# # # screen -r 28611
library(datelife)
library(ape)
library(geiger)
library(testthat)
library(microbenchmark)
setwd("~/Google Drive/datelife/runtime_tests")
load("spp1500.RData")
aves1500.1.gfr.runtime_2017.12.28 <- microbenchmark(get_datelife_result(input=spp1500), times=100L)
save(aves1500.1.gfr.runtime_2017.12.28, file="aves1500.1.gfr.runtime_2017.12.28.RData")
# # # # finished ok

# # # screen -r 28624
library(datelife)
library(ape)
library(geiger)
library(testthat)
library(microbenchmark)
setwd("~/Google Drive/datelife/runtime_tests")
load("spp2000.RData")
aves2000.1.gfr.runtime_2017.12.28 <- microbenchmark(get_datelife_result(input=spp2000), times=100L)
save(aves2000.1.gfr.runtime_2017.12.28, file="aves2000.1.gfr.runtime_2017.12.28.RData")
# # # # finished ok

# # # screen -r 28636
library(datelife)
library(ape)
library(geiger)
library(testthat)
library(microbenchmark)
setwd("~/Google Drive/datelife/runtime_tests")
load("spp3000.RData")
aves3000.1.gfr.runtime_2017.12.28 <- microbenchmark(get_datelife_result(input=spp3000), times=100L)
save(aves3000.1.gfr.runtime_2017.12.28, file="aves3000.1.gfr.runtime_2017.12.28.RData")
# # # # finished ok

# # # screen -r 28644
library(datelife)
library(ape)
library(geiger)
library(testthat)
library(microbenchmark)
setwd("~/Google Drive/datelife/runtime_tests")
load("spp5000.RData")
aves5000.1.gfr.runtime_2017.12.28 <- microbenchmark(get_datelife_result(input=spp5000), times=100L)
save(aves5000.1.gfr.runtime_2017.12.28, file="aves5000.1.gfr.runtime_2017.12.28.RData")
# # # # finished ok

# # # screen -r 28650
library(datelife)
library(ape)
library(geiger)
library(testthat)
library(microbenchmark)
setwd("~/Google Drive/datelife/runtime_tests")
load("spp7000.RData")
aves7000.1.gfr.runtime_2017.12.28 <- microbenchmark(get_datelife_result(input=spp7000), times=100L)
save(aves7000.1.gfr.runtime_2017.12.28, file="aves7000.1.gfr.runtime_2017.12.28.RData")
# # # # finished ok

# # # screen -r 28660
library(datelife)
library(ape)
library(geiger)
library(testthat)
library(microbenchmark)
setwd("~/Google Drive/datelife/runtime_tests")
load(file="aves.spp.RData")
aves.all.gfr.runtime_2017.12.29 <- microbenchmark(get_datelife_result(input=aves.spp$cleaned.names), times=100L)
save(aves.all.gfr.runtime_2017.12.29, file="aves.all.gfr.runtime_2017.12.29.RData")
# # # # finished ok

# # # screen -r 31197
library(datelife)
library(ape)
library(geiger)
library(testthat)
library(microbenchmark)
setwd("~/Google Drive/datelife/runtime_tests")
load("spp8000.RData")
aves8000.1.gfr.runtime_2017.12.28 <- microbenchmark(get_datelife_result(input=spp8000), times=100L)
save(aves8000.1.gfr.runtime_2017.12.28, file="aves8000.1.gfr.runtime_2017.12.28.RData")

# # # screen -r 29268
library(datelife)
library(ape)
library(geiger)
library(testthat)
library(microbenchmark)
setwd("~/Google Drive/datelife/runtime_tests")
load("spp9000.RData")
aves9000.1.gfr.runtime_2017.12.28 <- microbenchmark(get_datelife_result(input=spp9000), times=100L)
save(aves9000.1.gfr.runtime_2017.12.28, file="aves9000.1.gfr.runtime_2017.12.28.RData")
# # # # finished ok

# # # screen -r 31204
library(datelife)
library(ape)
library(geiger)
library(testthat)
library(microbenchmark)
setwd("~/Google Drive/datelife/runtime_tests")
load("spp10000.RData")
aves10000.1.gfr.runtime_2017.12.28 <- microbenchmark(get_datelife_result(input=spp10000), times=100L)
save(aves10000.1.gfr.runtime_2017.12.28, file="aves10000.1.gfr.runtime_2017.12.28.RData")


# # # analyse results
ninput <- c(10,100,200,300,400,500,700,1000,1500,2000,3000,5000,7000, 8000, 9000, 10000)
res <- c()
for(i in ninput){
	xname <- paste0("aves",i,".1.gfr.runtime_2017.12.28")
	x <- paste0(xname, ".RData")
	load(x)
	res <- rbind(res, get(xname))
}
length(res)
str(res)

load('~/Google Drive/datelife/runtime_tests/aves.all.gfr.runtime_2017.12.29.RData')
aves.all.gfr.runtime_2017.12.29
res <- rbind(res, aves.all.gfr.runtime_2017.12.29)
autoplot(res)
plot(res)
?autoplot
str(res)
res$time <- res$time*0.001 # converts units to microseconds, so not useful
p <- autoplot(res)
str(p)
aves7000.gfr <- get_datelife_result(input=spp7000)
class(aves7000.gfr)
length(aves7000.gfr)
str(aves7000.gfr[1])
dim(aves7000.gfr)[1]

# # # Rerun everything, cause weirdly flat:
# # # # tests for changing level names
str(aves10.1.gfr.runtime_2017.12.29)
length(aves10.1.gfr.runtime_2017.12.29)
names(aves10.1.gfr.runtime_2017.12.29)
levels(aves10.1.gfr.runtime_2017.12.29$expr)
# # # # screen -r 
library(datelife)
library(ape)
library(geiger)
library(testthat)
library(microbenchmark)
ninput <- c(10,100,200,300,400,500,700,1000,1500,2000,3000,5000,7000,8000,9000,10000)
setwd("~/Google Drive/datelife/runtime_tests")
for(i in ninput){
	xname <- paste0("spp",i)
	load(paste0(xname,".RData"))
	x <- microbenchmark(get_datelife_result(input=get(xname), process_input=TRUE), times=100L) # input must be processed :)
	# y <- levels(x$expr)
	# levels(x$expr)[levels(x$expr==y)] <- paste0(i, " names")
	levels(x$expr)[1] <- paste0(i, " names")
	xnameobj <- paste0("aves",i,".1.gfr.runtime_2017.12.29")
	assign(xnameobj, x)
	save(list=xnameobj, file=paste0(xnameobj,".RData"))
	rm(list=xnameobj)
}

# # # We could do 9 more runs with different spp samples
# # # # Ask the team:
# # # # They said 100 repetitions is ok, but that names should be randomized every run. Done later...
# setwd("~/Google Drive/datelife/runtime_tests")
# load(file="seeds.RData")
# seeds
# load('~/Google Drive/datelife/runtime_tests/aves.spp.RData')
# ninput <- c(10,100,200,300,400,500,700,1000,1500,2000,3000,5000,7000,8000,9000,10000)
# for (j in 2:10){	
	# for (i in ninput){
		# set.seed(seeds[j])
		# x <- sample(aves.spp$cleaned.names, i)
		# xname <- paste0("spp",i, ".",j)
		# assign(xname, x)
		# save(list=xname, file=paste0(xname,".RData"))
	# }
# }


