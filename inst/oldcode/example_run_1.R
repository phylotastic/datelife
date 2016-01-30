library(PhyloOrchard) #install.packages("PhyloOrchard", repos="http://R-Forge.R-project.org")
library(phylobase)
data(BinindaEmondsEtAl2007)
best.phylo4<-BinindaEmondsEtAl2007$trees$best
best.phylo3<-reorder(as(best.phylo4,"phylo"))

patristic.distance <- ComputePatristicDistance(best.phylo3)

# this should work
mrca.age <- GetAge(GetSubsetMatrix( patristic.distance, c("Homo_sapiens","Pan_paniscus"))$patristic.matrix)
print(mrca.age)

# so will this, with three taxa
mrca.age <- GetAge(GetSubsetMatrix( patristic.distance, c("Homo_sapiens", "Pan_paniscus", "Mus_musculus"))$patristic.matrix)
print(mrca.age)

# if we want all the pairwise dates, just print the output matrix, halved (patristic distance = twice the depth)
print(GetSubsetMatrix( patristic.distance, c("Homo_sapiens", "Pan_paniscus", "Mus_musculus"))$patristic.matrix/2)

# now try to find a group with a missing species:
result<-GetSubsetMatrix( patristic.distance, c("Homo_sapiens", "Pan_paniscus","Unicorn"))
print(result$problem)
print(GetAge(result$patristic.matrix))

# now have insufficient overlap even for an estimate
result<-GetSubsetMatrix( patristic.distance, c("Homo_sapiens", "Griffin","Unicorn"))
print(result$problem)
print(GetAge(result$patristic.matrix))

# now do all three Bininda-Emonds trees:
bininda.emonds.all<-BindMatrices(list(
  ComputePatristicDistance(reorder(as(BinindaEmondsEtAl2007$trees$best,"phylo")) ),
  ComputePatristicDistance(reorder(as(BinindaEmondsEtAl2007$trees$upper,"phylo")) ),
  ComputePatristicDistance(reorder(as(BinindaEmondsEtAl2007$trees$lower,"phylo")) )
  ))
result<-GetSubsetArray (bininda.emonds.all, c("Homo_sapiens", "Pan_paniscus", "Mus_musculus"))
print(GetAges(result$patristic.matrix.array))

# load 20K primate trees from Tracy Heath
data(Heath)
heath.all<-BindMatrices(lapply(sample(Heath$trees,500),ComputePatristicDistance))
#heath.all<-BindMatrices(lapply(Heath$trees,ComputePatristicDistance)) #will be slow
result<-GetSubsetArray (heath.all, c("Lemur_catta", "Galago_crassicaudatus"))
print(GetAges(result$patristic.matrix.array))
print(GetQuantiles(GetAges(result$patristic.matrix.array)))
