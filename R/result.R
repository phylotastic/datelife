run<-function(input=c("Rhinoceros_unicornis","Equus_caballus"), format="html", partial="liberal",useembargoed="yes", uncertainty=100, randomtreesperstudy=0, plot.width=600, plot.height=600, usetnrs="no", tnrssource="NCBI", version="stable") {
  if(version=="stable") {
    #remember we have from datelifeStarter.R the vectors citations and embargoed and the list of patristic.matrix.arrays
    #  studies
    phy<-NULL
    input<-gsub("\\+"," ",input)
    input<-str_trim(input, side = "both")
 
    if(grepl('\\(', input) & grepl('\\)', input) & (substr(input,nchar(input),nchar(input))==";")) { #our test for newick
      phy<-read.tree(text=input)
    }
    cleaned.names<-""
    if(!is.null(phy)) {
      if(usetnrs=="yes") {
        phy <- suppressMessages(checknames(phylo=phy, source_=tnrssource, byfilename=FALSE))
      }
      cleaned.names<-phy$tip.label 
    } else {
      cleaned.names<-strsplit( gsub("\\s","",input), ",")[[1]]
      #cleaned.names<-lapply(cleaned.names, str_trim, side="both")
      if (usetnrs=="yes") {
        cleaned.names <- checknames(charvector=cleaned.names, source_=tnrssource)
      }
    }
    
    if(format=="newick1000") {
      randomtreesperstudy<-1000
    }
    tree.list<-list()
    results.list<-lapply(studies,GetSubsetArrayDispatch, taxa=cleaned.names, phy=phy) 
    median.patristic.matrices<-list()
    ages.matrix<-c() #will hold median, and 95% CI
    uncertainty<-as.numeric(uncertainty)/100 #make percentage
    used.studies<-c()
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
      used.studies<-c(used.studies,display.result)
      if (display.result) { 
        ages<-GetAges(result$patristic.matrix.array)
        if(result$problem=="none") { #rather than dealing with missing taxa
          median.patristic.matrices[[length(median.patristic.matrices)+1]] <- SummaryPatristicMatrix(result$patristic.matrix.array)
        }
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
        if (format=="newick1000") {
          SamplePatristicMatrix <- function(patristic.matrix.array, uncertainty) {
            # if (dim(patristic.matrix.array)[3] == 1) {
            # 	patristic.matrix<-patristic.matrix.array[,,1]
            # 	#need order of node depths, from just the upper triangular and diagonal part of the matrix
            # 	element.order<-order(patristic.matrix[upper.tri(patristic.matrix,diag=FALSE)],decreasing=TRUE)
            # 	new.patristic.matrix<-patristic.matrix*0
            # 	cur.val<-patristic.matrix[upper.tri(patristic.matrix,diag=FALSE)][element.order[1]]
            #   new.patristic.matrix[upper.tri(new.patristic.matrix,diag=FALSE)][element.order[1]] <- cur.val + runif(1, -cur.val*uncertainty/100, cur.val*uncertainty/100)
            #	element.order<-element.order[-1]
            #  	for (i in sequence(length(element.order))) {
            #  		cur.val<-patristic.matrix[upper.tri(patristic.matrix,diag=FALSE)][element.order[i]]
            #  		new.patristic.matrix[upper.tri(new.patristic.matrix,diag=FALSE)][element.order[i]] <- cur.val + runif(1, -cur.val*uncertainty/100, min(cur.val*uncertainty/100, min( ))
            #  	}
            #  }
            #  else {
            return(patristic.matrix<-patristic.matrix.array[,,sample.int(1, size=dim(patristic.matrix.array)[3], replace=TRUE )] )
            # }
          }
          for (rep in sequence(5) ) {
            #out("try")
            out(write.tree(PatristicMatrixToTree(SamplePatristicMatrix(result$patristic.matrix.array, uncertainty=uncertainty))))
          }
        }
      }
    }
    median.patristic.matrix<-matrix(ncol=1,nrow=1)
    if(length(median.patristic.matrices)>0) {
      median.patristic.matrix<-SummaryPatristicMatrix(BindMatrices(median.patristic.matrices))
    }
    if (format=="html") {
      out("</table></p>")
      out(paste("<p>The best guess (median) for the estimate is ",median(ages.matrix[,1])," MY, ",sep=""))
      out(paste("but the median uncertainty for age goes from ",median(ages.matrix[,2])," to ",median(ages.matrix[,3]),sep=""))
      out(paste("MY and the maximum uncertainty goes from ",min(ages.matrix[,2])," to ",max(ages.matrix[,3])," MY.",sep=""))
      if (dim(median.patristic.matrix)[1]>2) {
        out("<p>Newick tree string: based on median tree from each study (only those with no problems), then median of those:</p>")
        out(write.tree(PatristicMatrixToTree( median.patristic.matrix )))
        out("<p>")
        p<-WebPlot(as.numeric(plot.width), as.numeric(plot.height))
        plot(PatristicMatrixToTree( median.patristic.matrix ))
        axisPhylo()
        out(p)
        out("</p>")
      }
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
    if (format=="newickmed") {
      if (dim(median.patristic.matrix)[1]>2) {
        out(write.tree(PatristicMatrixToTree( median.patristic.matrix )))
      }
    }
    if (format=="png") {
      p<-WebPlot(as.numeric(plot.width), as.numeric(plot.height))
      plot(PatristicMatrixToTree( median.patristic.matrix ))
      axisPhylo()
      out(p)
    }
    return(done())
  }
  if (version=="bleedingedge") {
    source("/Library/WebServer/Sites/datelife.org/datelife/R/bleedingedge.R", local=TRUE) #has commands to use
  }
}
