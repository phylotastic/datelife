#' Core function to generate results
#' @param input A newick string or vector of taxa
#' @param format The output format
#' @param partial How to deal with trees that have a subset of taxa in the query
#' @param uncertainty How much to multiply the range by to represent uncertainty
#' @param randomtreesperstudy IDK
#' @param plot.width Width in pixels for output plot
#' @param plot.height Height in pixels for output plot
#' @param usetnrs Whether to use OpenTree's TNRS for the input
#' @param approximatematch IDK
#' @param pruneonmatch IDK
#' @param datelife.cache The list of lists containing the input trees and other info
#' @return Depends on options
#' @export
run<-function(input=c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"), format="html", partial="liberal", uncertainty=100, randomtreesperstudy=0, plot.width=600, plot.height=600, usetnrs="no", approximatematch="yes", prunenonmatch="yes", datelife.cache=datelife.cache, do.out = TRUE) {
  out.vector <- ""
 
    input.processed <- ProcessInput(input, usetnrs, approximatematch)
    phy <- input.processed$phy
    cleaned.names <- input.processed$cleaned.names


    tree.list<-list()
    results.list<-lapply(datelife.cache$trees,GetSubsetArrayDispatch, taxa=cleaned.names, phy=phy)
    filtered.results <- ProcessResultsList(results.list, taxa, partial)
    
    
    #NOW START PROCESSING THE FILTERED.RESULTS: GET MAX AGE, GET TREES, ETC.
    
    
    
    
    
    
    
    
    median.patristic.matrices<-list()
    ages.matrix<-c() #will hold median, and 95% CI
    uncertainty<-as.numeric(uncertainty)/100 #make percentage
    used.studies<-c()
    if (format=="html") {
    	out.vector <- append(out.vector, "<!doctype html><html>")
    	pagestart <- "START"
      #pagestart<-scan('/Library/WebServer/Sites/datelife.org/datelife/php/pagestart.html',"raw",sep="\n",quiet=TRUE)
      for(i in sequence(length(pagestart))) {
        out.vector <- append(out.vector, pagestart[i])
        out.vector <- append(out.vector, "\n")
      }
      out.vector <- append(out.vector, "<p>")
      out.vector <- append(out.vector, "<table border='1'><tr><th>Median</th><th>Min</th><th>2.5% quantile</th><th>97.5% quantile</th><th>Max</th><th>NTrees</th><th>Problems</th><th>Citation</th></tr>")

    }
    for (i in sequence(length(datelife.cache$trees))) {
      result<-results.list[[i]]
      num.matching<-0
      if (!is.na(result$patristic.matrix.array)) { #happens if 0 return
 #       num.matching<-dim(SplitArray(result$patristic.matrix.array)[[1]])[1] #take just the first array
 		num.matching <- dim(result$patristic.matrix.array)[1]
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
          out.vector <- append(out.vector, paste("\n<tr>",VectorToTableRow(GetQuantiles(ages,probs)),"<td>",length(ages),"</td><td>",result$problem,"</td><td>",names(datelife.cache$trees)[i],"</td></tr>",sep="",collapse=""))
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
            #out.vector <- append(out.vector, "try")
            out.vector <- append(out.vector, write.tree(PatristicMatrixToTree(SamplePatristicMatrix(result$patristic.matrix.array, uncertainty=uncertainty))))
          }
        }
      }
    }
    median.patristic.matrix<-matrix(ncol=1,nrow=1)
    if(length(median.patristic.matrices)>0) {
      median.patristic.matrix<-SummaryPatristicMatrix(BindMatrices(median.patristic.matrices))
    }
    if (format=="html") {
      out.vector <- append(out.vector, "</table></p>")
      out.vector <- append(out.vector, paste("<p>The best guess (median) for the estimate is ",median(ages.matrix[,1])," MY, ",sep=""))
      out.vector <- append(out.vector, paste("but the median uncertainty for age goes from ",median(ages.matrix[,2])," to ",median(ages.matrix[,3]),sep=""))
      out.vector <- append(out.vector, paste("MY and the maximum uncertainty goes from ",min(ages.matrix[,2])," to ",max(ages.matrix[,3])," MY.",sep=""))
      if (dim(median.patristic.matrix)[1]>2) {
        out.vector <- append(out.vector, "<p>Newick tree string: based on median tree from each study (only those with no problems), then median of those:</p>")
        out.vector <- append(out.vector, write.tree(PatristicMatrixToTree( median.patristic.matrix )))
        out.vector <- append(out.vector, "<p>")
        p<-WebPlot(as.numeric(plot.width), as.numeric(plot.height))
        plot(PatristicMatrixToTree( median.patristic.matrix ))
        axisPhylo()
        out.vector <- append(out.vector, p)
        out.vector <- append(out.vector, "</p>")
      }
      pageend <- "END"
      #pageend<-scan('/Library/WebServer/Sites/datelife.org/datelife/php/pageend.html',"raw",sep="\n",quiet=TRUE)
      for(i in sequence(length(pageend))) {
        out.vector <- append(out.vector, pageend[i])
        out.vector <- append(out.vector, "\n")
      }
    }
    if (format=="bestguess") {
      out.vector <- append(out.vector, median(ages.matrix[,1]))
    }
    if (format=="bestguessuncert") {
      out.vector <- append(out.vector, paste(median(ages.matrix[,1]),median(ages.matrix[,2]),median(ages.matrix[,3]),sep=","))
    }
    if (format=="newickmed") {
      if (dim(median.patristic.matrix)[1]>2) {
        out.vector <- append(out.vector, write.tree(PatristicMatrixToTree( median.patristic.matrix )))
      }
    }
    if (format=="png") {
      p<-WebPlot(as.numeric(plot.width), as.numeric(plot.height))
      plot(PatristicMatrixToTree( median.patristic.matrix ))
      axisPhylo()
      out.vector <- append(out.vector, p)
    }
    if(do.out) {
    	out(out.vector)
    	return(done())
    } else {
    	return(out.vector)
    }
  }