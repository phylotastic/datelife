library(plyr)
library(stringr)
#this is a temporary script to just put files with the final names. Replace when cleaning actually happens
data.dir<-"/Library/WebServer/Sites/datelife.org/phyloorchard_new/pkg/data/"
cleaned.dir<-"/Library/WebServer/Sites/datelife.org/datelife/data/"
treefilenames <- dir(data.dir)

# Load all trees into the workspace
# CHANGE PATH!!!!
#l_ply(treefilenames, function(x) load(paste(data.dir,x,sep=""), .GlobalEnv))
setwd(data.dir)

# Get tree names in the workspace
trees <- sapply(treefilenames, function(x) str_replace(x, ".rda", ""), USE.NAMES=F)
for (i in sequence(length(treefilenames))) {
  load(treefilenames[i]) 
  trees_cleaned<-eval(parse( text=trees[i] ))
  if(class(trees_cleaned)=="phylo") {
    trees_cleaned<-c(trees_cleaned) 
  }
  save(trees_cleaned, file=paste(cleaned.dir, trees[i], ".rda", sep=""), compress=FALSE)
  rm(trees_cleaned)
}