run<-function(taxa=c("Rhinoceros_unicornis","Equus_caballus"), format="html", partial="liberal",useembargoed="yes", uncertainty=100) {
  #remember we have from datelifeStarter.R the vectors citations and embargoed and the list of patristic.matrix.arrays
  #  studies
	cleaned.names<-strsplit( gsub("\\s","",taxa), ",")[[1]]
	results.list<-lapply(studies,GetSubsetArray, taxa=cleaned.names)
  ages.matrix<-c() #will hold median, and 95% CI
  uncertainty<-as.numeric(uncertainty)/100 #make percentage
  if (format=="html") {
    out("<!doctype html><html>")
    pagestart<-scan('/Library/WebServer/Sites/datelife.org/datelife/php/pagestart.html',"raw",sep="\n",quiet=TRUE)
    for(i in sequence(length(pagestart))) {
      out(pagestart[i])
      out("\n")
    }
    out("<p>")
    out("<table border='1'><tr><th>Median</th><th>Min</th><th>2.5% quantile</th><th>97.5% quantile</th><th>Max</th><th>NTrees</th><th>Problems</th><th>Citation</th></tr>")
  }
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
      if (length(ages)==1) {
        ages.matrix<-rbind(ages.matrix,matrix(c(ages[1],ages[1]-uncertainty*ages[1],ages[1]+uncertainty*ages[1]) ,nrow=1))
      }
      else {
        ages.matrix<-rbind(ages.matrix,matrix(quantile(ages,c(0.5,0.025,0.975)),nrow=1))
      }
      if (format=="html") {
        probs<-c(0.5,0,0.025,0.975,1)
        out(paste("\n<tr>",VectorToTableRow(GetQuantiles(ages,probs)),"<td>",length(ages),"</td><td>",result$problem,"</td><td>",citations[i],"</td></tr>",sep="",collapse=""))
      }
    }
  }
  if (format=="html") {
    out("</table></p>")
    out(paste("<p>The best guess (median) for the estimate is ",median(ages.matrix[,1])," MY, ",sep=""))
    out(paste("but the median uncertainty for age goes from ",median(ages.matrix[,2])," to ",median(ages.matrix[,3]),sep=""))
    out(paste("MY and the maximum uncertainty goes from ",min(ages.matrix[,2])," to ",max(ages.matrix[,3])," MY.",sep=""))
    pageend<-scan('/Library/WebServer/Sites/datelife.org/datelife/php/pageend.html',"raw",sep="\n",quiet=TRUE)
    for(i in sequence(length(pageend))) {
      out(pageend[i])
      out("\n")
    }
  }
  if (format=="bestguess") {
    out(median(ages.matrix[,1])) 
  }
	if (format=="bestguessuncert") {
	  out(paste(median(ages.matrix[,1]),median(ages.matrix[,2]),median(ages.matrix[,3]),sep=",")) 
	}
	
  return(done())

}
