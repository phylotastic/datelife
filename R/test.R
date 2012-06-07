run<-function(taxa=c("Homo_sapiens","Mus_musculus"), format="html", partial="liberal") {
	cleaned.names<-strsplit( gsub("\\s","",taxa), ",")[[1]]
	results<-GetSubsetMatrix( patristic.distance, cleaned.names)
  if (format=="html") {
    out("<!doctype html><html lang='en'><head><title>DateLife</title></head>")
    out("<body>")
  	oprint(cleaned.names)
  	oprint(paste("Problems: ",results$problem,sep=""))
  	oprint(paste("Age (MY): ",GetAge(results$patristic.matrix),sep=""))
  	if (dim(results$patristic.matrix)[1]>2) {
		  oprint(write.tree(PatristicMatrixToTree( results$patristic.matrix )))
	  }
    out("</body></html>")
    return(done())
  }
  if (format=="newick") {
    if (dim(results$patristic.matrix)[1]>2) {
      out(write.tree(PatristicMatrixToTree( results$patristic.matrix )))
      return(done())
    }
    else {
       out(paste("Error: only ",dim(results$patristic.matrix)[1]," taxa returned",sep=""))
       return(done())
    }
  }
  if (format=="bestguess") {
    out(GetAge(results$patristic.matrix))
    return(done())
  }
}
