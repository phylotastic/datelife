run<-function(taxa=c("Homo_sapiens","Mus_musculus"), format="html", partial="liberal",useembargoed="yes") {
  #remember we have from datelifeStarter.R the vectors citations and embargoed and the list of patristic.matrix.arrays
  #  studies
	cleaned.names<-strsplit( gsub("\\s","",taxa), ",")[[1]]
	results.list<-lapply(studies,GetSubsetArray, taxa=cleaned.names)
  if (format=="html") {
    out("<!doctype html><html lang='en'><head><title>DateLife</title></head>")
    out("<body><p>")
    out("<table border='1'><tr><td>Median</td><td>Min</td><td>2.5% quantile</td><td>97.5% quantile</td><td>Max</td><td<NTrees</td><td>Problems</td><td>Citation</td></tr>")
    for (i in sequence(length(studies))) {
      result<-results.list[[i]]
      num.matching<-0
      if (!is.na(result$patristic.matrix.array)) { #happens if 0 return
        num.matching<-dim(SplitArray(result$patristic.matrix.array)[[1]])[1] #take just the first array
      }
      display.result<-FALSE
      if(num.matching == length(cleaned.names)) {
        display.result<-TRUE 
      }
      else {
        if ( ( num.matching > 1)  && (partial=="liberal") ) {
          display.result<-TRUE 
        }
      }
      if ((useembargoed=="no") && (embargoed[i]==TRUE)) {
        display.result<-FALSE 
      }
      if (display.result) {
        ages<-GetAges(result$patristic.matrix.array)
        out("\n")
        probs<-c(0.5,0,0.025,0.975,1)
        out(paste("<tr>",VectorToTableRow(GetQuantiles(ages,probs)),"<td>",length(ages),"</td><td>",result$problem,"</td><td>",citations[i],"</td></tr>",sep="",collapse=""))
      }
    }
    out("</table></p></body></html>")
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
