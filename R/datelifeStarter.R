citations<-c()
embargoed<-c()
studies<-list()
original.wd<-getwd()
setwd("/Library/WebServer/Sites/datelife.org/datelife/data/")
#setwd("/Users/bomeara/Documents/MyDocuments/Active/datelife/data/")
print(system("ls",intern=TRUE))

print("beginning to load data")


citations<-c(citations,"Heath et al. 2012")
embargoed<-c(embargoed,FALSE)
data(HeathEtAl2012)
print(ls())
studies[[length(studies)+1]]<-HeathEtAl2012_tree
rm(HeathEtAl2012_tree)
print(citations)

citations<-c(citations,"10K trees Perissodactyla, V1")
embargoed<-c(embargoed,FALSE)
harvard10k.perissodactyla.trees<-read.nexus("TreeBlock_10kTrees_Perissodactyla_Version1.nex")
print(ls())
studies[[length(studies)+1]]<-harvard10k.perissodactyla.trees
rm(harvard10k.perissodactyla.trees)
print(citations)

citations<-c(citations,"Bininda-Emonds et al. 2007")
embargoed<-c(embargoed,FALSE)
data(BinindaEmondsEtAl2007)
print(ls())
studies[[length(studies)+1]]<-BinindaEmondsEtAl2007
rm(BinindaEmondsEtAl2007)
print(citations)

citations<-c(citations,"10K trees Carnivores, V1")
embargoed<-c(embargoed,FALSE)
harvard10k.carnivores.trees<-read.nexus("TreeBlock_10kTrees_Carnivora_Version1.nex")
print(ls())
studies[[length(studies)+1]]<-harvard10k.carnivores.trees
rm(harvard10k.carnivores.trees)
print(citations)

citations<-c(citations,"10K trees Primates, V3")
embargoed<-c(embargoed,FALSE)
harvard10k.primates.trees<-read.nexus("TreeBlock_10kTrees_Primates_Version3.nex")
print(ls())
studies[[length(studies)+1]]<-harvard10k.primates.trees
rm(harvard10k.primates.trees)
print(citations)

citations<-c(citations,"Alfaro et al. 2009")
embargoed<-c(embargoed,FALSE)
data(AlfaroEtAl2009)
print(ls())
studies[[length(studies)+1]]<-AlfaroEtAl2009_tree
rm(AlfaroEtAl2009_tree)
print(citations)

citations<-c(citations,"Bergmann and Irschick 2012")
embargoed<-c(embargoed,FALSE)
data(BergmannEtAl2012)
print(ls())
studies[[length(studies)+1]]<-BergmannEtAl2012
rm(BergmannEtAl2012)
print(citations)

# Eastman et al. unpublished tree of salamanders
citations<-c(citations,"embargoed")
embargoed<-c(embargoed,TRUE)
data(Unpub1)
print(ls())
studies[[length(studies)+1]]<-Unpub1
rm(Unpub1)
print(citations)

citations<-c(citations,"Oaks 2011")
embargoed<-c(embargoed,FALSE)
data(Oaks2012)
print(ls())
studies[[length(studies)+1]]<-Oaks2011
rm(Oaks2011)
print(citations)

citations<-c(citations,"Jaffe et al. 2011")
embargoed<-c(embargoed,FALSE)
data(JaffeEtAl2011)
print(ls())
studies[[length(studies)+1]]<-JaffeEtAl2011
rm(JaffeEtAl2011)
print(citations)

citations<-c(citations,"Zhang and Wake 2009")
embargoed<-c(embargoed,FALSE)
data(ZhangandWake2009)
print(ls())
studies[[length(studies)+1]]<-ZhangandWake2009
print(citations)

print(paste("number of studies = ",length(studies)))

print(paste("number of species = ",length(unique(unlist(sapply(studies,rownames))))))

function(ifelse(class(x)=="multiPhylo", length(x), 1L)
print(paste("number of trees total = ",sum(sapply(studies,dim)[3,])))

setwd(original.wd)

