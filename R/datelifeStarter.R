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
load("heath.rda")
print(ls())
studies[[length(studies)+1]]<-heath.all
rm(heath.all)
print(citations)

citations<-c(citations,"10K trees Perissodactyla, V1")
embargoed<-c(embargoed,FALSE)
load("harvard10k.perissodactyla.rda")
print(ls())
studies[[length(studies)+1]]<-harvard10k.perissodactyla.all
rm(harvard10k.perissodactyla.all)
print(citations)

if (1>2) {
citations<-c(citations,"Bininda-Emonds et al. 2007")
embargoed<-c(embargoed,FALSE)
load("bininda.emonds.rda")
print(ls())
studies[[length(studies)+1]]<-studies,bininda.emonds.all
rm(bininda.emonds.all)
print(citations)

citations<-c(citations,"10K trees Carnivores, V1")
embargoed<-c(embargoed,FALSE)
load("harvard10k.carnivores.rda")
print(ls())
studies[[length(studies)+1]]<-harvard10k.carnivores.all
rm(harvard10k.carnivores.all)
print(citations)

citations<-c(citations,"10K trees Primates, V3")
embargoed<-c(embargoed,FALSE)
load("harvard10k.primates.rda")
print(ls())
studies[[length(studies)+1]]<-harvard10k.primates.all
rm(harvard10k.primates.all)
print(citations)
}

setwd(original.wd)