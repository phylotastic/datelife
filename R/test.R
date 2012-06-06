run<-function(taxa="cat", ...) {
	cleaned.names<-strsplit( gsub("\\s","",taxa), ",")[[1]]
	oprint(cleaned.names)
	results<-GetSubsetMatrix( patristic.distance, cleaned.names)
	oprint(paste("Problems: ",results$problem,sep=""))
	oprint(paste("Age (MY): ",GetAge(results$patristic.matrix),sep=""))
}
