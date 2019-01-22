# Miercoles 28 de Marzo 2018

# 1. make up more chronograms for a fake database to test get_datelife_result runtime
# with different amounts of chronograms in cache
length(opentree_chronograms)
ls(opentree_chronograms)
length(opentree_chronograms$trees)

plot(opentree_chronograms$trees[[1]])
otol_chr_length <- sapply(opentree_chronograms$trees, function(x) length(x$tip.label))
names(otol_chr_length)
class(otol_chr_length)
max(otol_chr_length)  #48017
opentree_chronograms$trees[which.max(otol_chr_length)]  #tree of life
# $`Hedges, S. Blair, Julie Marin, Michael Suleski, Madeline Paymer, Sudhir Kumar. 2015. Tree of life reveals clock-like speciation and diversification. Molecular Biology and Evolution 32 (4): 835-845`
# 
# Phylogenetic tree with 48017 tips and 48000 internal nodes.
# 
# Tip labels:
#   Thermotoga maritima, Deinococcus radiodurans, Leifsonia xyli, Bifidobacterium longum, Frankia sp. CcI3, Corynebacterium glutamicum, ...
# 
# Rooted; includes branch lengths.

min(otol_chr_length)  #6
opentree_chronograms$trees[which.min(otol_chr_length)]  # flamingos
# $`Torres, Chris R, Lisa M Ogawa, Mark AF Gillingham, Brittney Ferrari, Marcel van Tuinen. 2014. A multi-locus inference of the evolutionary diversification of extant flamingos (Phoenicopteridae). BMC Evolutionary Biology 14 (1): 36.`
# 
# Phylogenetic tree with 6 tips and 5 internal nodes.
# 
# Tip labels:
#   [1] "Phoenicopterus ruber"     "Phoenicopterus roseus"   
# [3] "Phoenicopterus chilensis" "Phoenicoparrus minor"    
# [5] "Phoenicoparrus andinus"   "Phoenicoparrus jamesi"   
# 
# Rooted; includes branch lengths.

hist(otol_chr_length)
ocl_table <- table(otol_chr_length)
plot(ocl_table)
plot(ocl_table[-length(ocl_table)])
ocl_table*10/2
names(which.max(ocl_table))
sum(ocl_table*10/2)
ocl_table["6670"]
x <- 2:500
i <- "6670"
prob <- 1:3
size = 100
sampleChr <- function(size, x, prob = NULL, otol_chr=opentree_chronograms){
  if(is.null(prob)) {
    a <- rep(0.01, length(x))
    names(a) <- x
    otol_chr_length <- sapply(otol_chr$trees, function(x) length(x$tip.label))
    ocl_table <- table(otol_chr_length)
    m <- as.numeric(names(which.max(ocl_table)))
    b <- ceiling(abs(rnorm(n = size, mean = m, sd = m/2)))
    b_table <- table(b)
    a[names(ocl_table)] <- a[names(ocl_table)] + ocl_table
    a[names(b_table)] <- a[names(b_table)] + b_table
    a <- a[!is.na(a)]
    prob <- a/sum(a)
  }
  print(length(x))
  print(length(prob))
  res <- sample(x = x, size = size, replace = T, prob = prob)
  return(res)
}
length(b_table)
length(ocl_table)
cache_size <- c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000)
x <- 5:50000
a <- rep(1, length(x))
names(a) <- x
otol_chr_length <- sapply(otol_chr$trees, function(x) length(x$tip.label))
ocl_table <- table(otol_chr_length)
a[names(ocl_table)] <- a[names(ocl_table)] + ocl_table
a <- a[!is.na(a)]
length(a)
prob <- a/sum(a)
sampleChr(size = 100, x = x)


sampleChr2 <- function(size){
  s90 <- floor(0.9*size)
  s5 <- floor(0.5*(size-s90))
  # s4 <- floor(0.4*(size-s90))
  # s3 <- floor(0.5*(size-s90-s4))
  s3 <- floor(0.6*(size-s90-s5))
  # s2 <- floor(0.67*(size-s90-s4-s3))
  # s1 <- size-s90-s4-s3-s2
  s1a <- floor(0.5*(size-s90-s5-s3))
  s1b <- size-s90-s5-s3-s1a
  
  s90 <- ceiling(abs(rnorm(n = s90, mean = 100, sd = 300))) + 1
  ceiling(abs(rnorm(n = s4, mean = 1500, sd = 500)))
  # s4 <- sample(x = 1001:3000, size = s4)
  s5 <- sample(x = 1001:3000, size = s5)
  s3 <- sample(x = 3001:10000, size = s3)
  s1a <- sample(x = 10001:30000, size = s1a)
  s1b <- sample(x = 30001:50000, size = s1b)
  # s1 <- sample(x = 30001:50000, size = s1)  # usaré el que ya esta hecho
  # ss <- c(s90, s4, s3, s2, s1)
  ss <- c(s90, s5, s3, s1a, s1b)
  return(ss)
}
sampleChr2(100)
library(TreeSim)
?TreeSim
t1 <- sim.bd.taxa(n = 50000, numbsim = 1, lambda = 1, mu = 0)
# it was quite fast, so let's simulate the big trees too
pdf(file = "~/Google Drive/datelife/runtime_tests/chr_50k_tips.pdf", height = 200, width = 20)
plot(t1[[1]], cex = 0.5, edge.width = 0.3)
dev.off()
max(branching.times(t1[[1]]))  #[1] 9.646046
# plot it to see if lambda 1 looks good. It looks good.

sim_chronograms_size <- c()
for (i in 1:10){
  sim_chronograms_size <- c(sim_chronograms_size, sampleChr2(size = 100))
}

sapply(sim_chronograms_size[1:10], function (x) sim.bd.taxa(n = x, numbsim = 1, lambda = 1, mu = 0))
sim_chronograms <- sapply(sim_chronograms_size, function (x) sim.bd.taxa(n = x, numbsim = 1, lambda = 1, mu = 0))
sim_chronograms_2018.03.28 <- sim_chronograms
save(sim_chronograms_2018.03.28, file = "~/Google Drive/datelife/runtime_tests/sim_chronograms_2018.03.28.RData")

# To DO:
# A. remake runtime tests for make datelife query with tnrs
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
  xnameobj <- paste0("make_datelife_query_runtime_tnrs_2018.03._", i,"_aves_spp")
  assign(xnameobj, y)
  setwd("~/Google Drive/datelife/runtime_tests/2_tests/2_random_spp_names/2_make_datelife_query")
  save(list=xnameobj, file=paste0(xnameobj,".RData"))
  rm(list=xnameobj)
}

# B. get_datelife_result() with already processed names and current number of chronograms in otol
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
}

# show a figure with different amount of data in chronograms cache
# there are at least 2 bird names in 68 out of 202 chronograms
# that's 68/202 = 33.66% of all chronograms
# chr_coverage <-  c(10, 20, 30, 40, 50)
# also, amount of names in chronogram might be an important parameter

