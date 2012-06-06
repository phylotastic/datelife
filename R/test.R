run<-function(taxa="cat", ...) {
	cleaned.names<-strsplit( gsub("\\s","",taxa), ",")[[1]]
	oprint(cleaned.names)
	results<-GetSubsetMatrix( patristic.distance, cleaned.names)
	oprint(paste("Problems: ",results$problem,sep=""))
	oprint(paste("Age (MY): ",GetAge(results$patristic.matrix),sep=""))
	if (dim(results$patristic.matrix)[1]>2) {
		oprint(write.tree(PatristicMatrixToTree( results$patristic.matrix )))
	}
	done()
}
