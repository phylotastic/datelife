data(BinindaEmondsEtAl2007)
best.phylo4<-BinindaEmondsEtAl2007$trees$best
best.phylo3<-reorder(as(best.phylo4,"phylo"))

patristic.distance <- ComputePatristicDistance(best.phylo3)

