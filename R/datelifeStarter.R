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

citations<-c(citations,"Bininda-Emonds et al. 2007")
embargoed<-c(embargoed,FALSE)
load("bininda.emonds.rda")
print(ls())
studies[[length(studies)+1]]<-bininda.emonds.all
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

citations<-c(citations,"Alfaro et al. 2009")
embargoed<-c(embargoed,FALSE)
load("alfaro2009.rda")
print(ls())
studies[[length(studies)+1]]<-alfaro2009.all
rm(alfaro2009.all)
print(citations)

citations<-c(citations,"Bergmann and Irschick 2012")
embargoed<-c(embargoed,FALSE)
load("bergmann.rda")
print(ls())
studies[[length(studies)+1]]<-bergmann.all
rm(bergmann.all)
print(citations)

citations<-c(citations,"embargoed")
embargoed<-c(embargoed,TRUE)
load("embargoed2.rda")
print(ls())
studies[[length(studies)+1]]<-embargoed2.all
rm(embargoed2.all)
print(citations)

citations<-c(citations,"Oaks 2011")
embargoed<-c(embargoed,FALSE)
load("oaks2011.rda")
print(ls())
studies[[length(studies)+1]]<-oaks2011.all
rm(oaks2011.all)
print(citations)

citations<-c(citations,"Jaffe et al. 2011")
embargoed<-c(embargoed,FALSE)
load("jaffe.rda")
print(ls())
studies[[length(studies)+1]]<-jaffe.all
rm(jaffe.all)
print(citations)

citations<-c(citations,"Zhang and Wake 2009")
embargoed<-c(embargoed,FALSE)
load("zhangwake2009.rda")
print(ls())
studies[[length(studies)+1]]<-zhangwake2009.all
rm(zhangwake2009.all)
print(citations)

print(paste("number of studies = ",length(studies)))

print(paste("number of species = ",length(unique(unlist(sapply(studies,rownames))))))

print(paste("number of trees total = ",sum(sapply(studies,dim)[3,])))

setwd(original.wd)