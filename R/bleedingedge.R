    #remember we have from datelifeStarter.R the vectors citations and embargoed and the list of patristic.matrix.arrays
    #  studies

source("/Library/WebServer/Sites/datelife.org/datelife/R/datelife.R")
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
    return(lapply(studies, GetSubsetArrayDispatch, taxa=cleaned.names))
