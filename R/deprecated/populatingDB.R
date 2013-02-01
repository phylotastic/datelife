setwd("/Library/WebServer/Sites/datelife.org/datelife/data")
library(PhyloOrchard) #install.packages("PhyloOrchard", repos="http://R-Forge.R-project.org")
library(phylobase)
source("../R/datelife.R")


#harvard10k.carnivores.trees<-read.nexus("TreeBlock_10kTrees_Carnivora_Version1.nex")
#harvard10k.perissodactyla.trees<-read.nexus("TreeBlock_10kTrees_Perissodactyla_Version1.nex")
#harvard10k.primates.trees<-read.nexus("TreeBlock_10kTrees_Primates_Version3.nex")
#harvard10k.carnivores.all<-BindMatrices(lapply(harvard10k.carnivores.trees,ComputePatristicDistance))
#harvard10k.perissodactyla.all<-BindMatrices(lapply(harvard10k.perissodactyla.trees,ComputePatristicDistance))
#harvard10k.primates.all<-BindMatrices(lapply(harvard10k.primates.trees,ComputePatristicDistance))
#save(harvard10k.carnivores.all,file="harvard10k.carnivores.rda",compress=FALSE)
#save(harvard10k.perissodactyla.all,file="harvard10k.perissodactyla.rda",compress=FALSE)
#save(harvard10k.primates.all,file="harvard10k.primates.rda",compress=FALSE)
#data(BinindaEmondsEtAl2007)
#bininda.emonds.all<-BindMatrices(list(
#  ComputePatristicDistance(reorder(as(BinindaEmondsEtAl2007$trees$best,"phylo")) ),
#  ComputePatristicDistance(reorder(as(BinindaEmondsEtAl2007$trees$upper,"phylo")) ),
#  ComputePatristicDistance(reorder(as(BinindaEmondsEtAl2007$trees$lower,"phylo")) )
#  ))
#save(bininda.emonds.all,file="bininda.emonds.rda",compress=FALSE)
#data(Heath)
#heath.all<-BindMatrices(lapply(sample(Heath$trees,1000),ComputePatristicDistance))
#save(heath.all,file="heath.rda",compress=FALSE)


#load('/Library/WebServer/Sites/datelife.org/phyloorchard/pkg/branches/data/RaboskyEtAlUnpublished.rda')
#embargoed1.all<-BindMatrices(list(ComputePatristicDistance(RaboskyEtAlUnpublished_tree)))
#save(embargoed1.all,file="embargoed1.rda",compress=FALSE)

load('/Library/WebServer/Sites/datelife.org/phyloorchard/pkg/branches/data/ZhangandWake2009.rda')
zhangwake2009.all<-BindMatrices(list(ComputePatristicDistance(as(ZhangandWake2009_tree,"phylo"))))
save(zhangwake2009.all,file="zhangwake2009.rda",compress=FALSE)

load('/Library/WebServer/Sites/datelife.org/phyloorchard/pkg/branches/data/EastmanEtAlUnpublished.rda')
embargoed2.all<-BindMatrices(list(ComputePatristicDistance(EastmanEtAlUnpublished_tree)))
save(embargoed2.all,file="embargoed2.rda",compress=FALSE)

load('/Library/WebServer/Sites/datelife.org/phyloorchard/pkg/branches/data/Oaks2011.rda')
oaks2011.all<-BindMatrices(list(ComputePatristicDistance(as(Oaks2011_tree,"phylo"))))
save(oaks2011.all,file="oaks2011.rda",compress=FALSE)

load('/Library/WebServer/Sites/datelife.org/phyloorchard/pkg/branches/data/BergmannEtAl2012.rda')
bergmann.all<-BindMatrices(list(ComputePatristicDistance(as(BergmannEtAl2012_tree,"phylo"))))
save(bergmann.all,file="bergmann.rda",compress=FALSE)

load('/Library/WebServer/Sites/datelife.org/phyloorchard/pkg/branches/data/AlfaroEtAl2009.rda')
alfaro2009.all<-BindMatrices(list(ComputePatristicDistance(as(AlfaroEtAl2009_tree,"phylo"))))
save(alfaro2009.all,file="alfaro2009.rda",compress=FALSE)

load('/Library/WebServer/Sites/datelife.org/phyloorchard/pkg/branches/data/JaffeEtAl2011.rda')
jaffe.all<-BindMatrices(list(ComputePatristicDistance(as(JaffeEtAl2011_tree,"phylo"))))
save(jaffe.all,file="jaffe.rda",compress=FALSE)

load('/Library/WebServer/Sites/datelife.org/phyloorchard/pkg/branches/data/megatree.rda')
embargoed3.all<-BindMatrices(list(ComputePatristicDistance(as(megaTree_tree,"phylo"))))
save(embargoed3.all,file="embargoed3.rda",compress=FALSE)

