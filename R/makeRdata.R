citations<-c()
embargoed<-c()
studies<-list()
all.taxa<-c()
ntrees<-0
ntrees.per.study<-c()
ntax.per.study<-c()
max.trees.per.study<-100 #some studies will have a huge number of trees. They will be consistently sampled up to this number
max.array.size<-14000 #Any dataset for which ntax * ntrees is bigger than this will be left in smaller, but slower, multiphylo format

GetFieldFromRdFile <- function (file, field="source", sep="") {
  all.text<-read.delim(file=file, stringsAsFactors=FALSE)[,1]
  focal.start<-which(all.text==paste("\\",field,"{",sep=""))+1
  if(length(focal.start)==1) { #it is defined
     final.text<-""
     while(all.text[focal.start]!="}") {
        final.text<-paste(final.text, all.text[focal.start], sep=sep)
        focal.start<-focal.start+1
     }
     return(final.text)
  }
  return(NULL)
}

cleaned.data.dir <- "/Library/WebServer/Sites/datelife.org/datelife/data/"
man.data.dir <- "/Library/WebServer/Sites/datelife.org/phyloorchard/pkg/man/"
#cleaned.data.dir <- "/Users/bomeara/Desktop/phyloorchard/pkg/cleaned/"
#man.data.dir <- "/Users/bomeara/Desktop/phyloorchard/pkg/man/"

#treefilenames <- dir(cleaned.data.dir)
treefilenames <- system(paste("ls -1S ", cleaned.data.dir), intern=TRUE)

trees <- sapply(treefilenames, function(x) str_replace(x, ".rda", ""), USE.NAMES=F)
for (i in sequence(length(treefilenames))) {
  load(paste(cleaned.data.dir,treefilenames[i], sep="")) 
  if(class(trees_cleaned)=="phylo") {
    trees_cleaned<-c(trees_cleaned) 
  }
  all.taxa<-append(all.taxa, trees_cleaned[[1]]$tip.label)
  ntrees<-ntrees+min(max.trees.per.study, length(trees_cleaned))
  new.pos <- length(studies)+1
  trees_cleaned<-trees_cleaned[unique(round(seq(from=1, to=length(trees_cleaned), length.out=max.trees.per.study)))]
  ntrees.per.study[new.pos] <- length(trees_cleaned)
  ntax.per.study[new.pos] <- Ntip(trees_cleaned[[1]])
  
  print(paste("now processing ",treefilenames[i]," with ",length(trees_cleaned)," trees each with ",Ntip(trees_cleaned[[1]])," taxa", sep=""))
  
  studies[[new.pos]]<-NA
  if (Ntip(trees_cleaned[[1]]) * length(trees_cleaned) <= max.array.size) {
  	try(studies[[new.pos]] <- BindMatrices(mclapply(trees_cleaned,ComputePatristicDistance, mc.cores=min(1, detectCores()-2))), silent=TRUE)
  }
  if (is.na(studies[[new.pos]])) {
    studies[[new.pos]] <- trees_cleaned #couldn't make the array, perhaps too big for memory or for R. Stick to slower trees
  } 
  print("loaded")
  new.citation<-NULL
  
  try(new.citation <- GetFieldFromRdFile(file=paste(man.data.dir, trees[i], ".Rd", sep="", collapse="")), silent=TRUE)
  if (is.null(new.citation)) {
    matching.files<-system(paste("grep -l ", trees[i], " ", man.data.dir, "*", sep=""), intern=TRUE)
    try(new.citation <- GetFieldFromRdFile(file=matching.files[1]), silent=TRUE)
  }
  citations[new.pos] <- new.citation
  embargoed[new.pos] <- FALSE
  if (grepl("Unpub",trees[i])) {
  	embargoed[new.pos] <- TRUE
  }
  rm(trees_cleaned)
}

unique.taxa <- unique(all.taxa)

print(citations)

print(paste("number of studies = ",length(studies)))

print(paste("number of taxa = ", length(unique.taxa)))

print(paste("number of trees = ", ntrees))

save(list=ls(), file=paste(cleaned.data.dir,"StartFiles.Data"), compress=TRUE) #Note that this is saved with a different extension