# Including TNRS in Datelife

#' Find only rows in a data.frame where the values don't match in a row
#' 
#' @import plyr
#' @param x A data.frame
#' @examples
#' dff <- data.frame(FabreEtAl2009$tip.label, out$tip.label)
#' nonmatching_by_row(dff)
nonmatching_by_row <- function(df){
	does_match <- function(y){ ifelse(y[,1] %in% y[,2], 0, 1) }
	df2 <- data.frame(old=sort(df[,1]), new=sort(df[,2]))
	temp <- adply(df2, 1, does_match)
	outt <- temp[temp$V1 == 1,][,-3]
	if(nrow(outt)==0){
		message("all taxa matched, no changes")
	} else
	{ outt }
}

#' Function to break up a vector by certain N sized chunks.
#' 
#' @import plyr
#' @param input Input character vector
#' @param by Number by which to split the input character vector
slice <- function(input, by = 2) {
	starts <- seq(1, length(input), by)
	tt <- lapply(starts, function(y) input[y:(y + (by - 1))])
	llply(tt, function(x) x[!is.na(x)])
}

#' Function to detect messy names and remove them, and remove sp.'s
#' 
#' @import stringr
#' @param taxon Taxonomic name, or a vector of taxonomic names.
#' @examples \dontrun{
#' detect_unid(taxon = Pyron2011[[1]]$tip.label)
#' }
detect_unid <- function(taxon){
	# replace underscores
	taxon2 <- str_replace_all(taxon, "_", " ")
	
	# split species names to find duplicates in next step
	if(all(sapply(taxon2, function(x) length(str_split(x, " ")[[1]]), USE.NAMES=F) == 1)){
		dups <- taxon2
	} else
	{
		split_collapse <- function(x){
			y <- str_split(x, " ")[[1]]
			if(length(y)>1){
				paste(y[1:2], collapse=" ",sep="")	
			} else { y }
		}
		dups <- sapply(taxon2, split_collapse, USE.NAMES=F)
	}
	
	# remove duplicates
	numdups <- split(dups, as.factor(dups))
	makenull <- function(x){
		if(length(x)>1) {
			x[1:length(x)-1] <- rep("NULL NULL", length(x)-1)
			x
		} else{ x }
	}
	dupsout <- unlist(lapply(numdups, makenull), use.names=F)
	# 	dupsout <- unique(dups)
	# remove epithets with just "cf"
	# 	nocfs <- dupsout[grep("^cf$", sapply(dupsout, function(x) str_split(x, " ")[[1]][[2]], USE.NAMES=F), invert=T)]
	
	# remove sp.'s, etc. 
	# 	remove_sps <- gsub(" sp[^aeiouh]+| sp\\.+| sp$", "", dupsout)
	
	# remove epithets with just one letter
	# 	remove_onechar_epithet <- gsub(" [A-Za-z0-9]{1}$", "", remove_sps)
	
	# remove genera (without epithets) if matches in the vector
	# 	remove_sps_u <- unique(remove_onechar_epithet)
	# 	justfirst <- sapply(remove_sps_u, function(x) str_split(x, " ")[[1]][[1]], USE.NAMES=F)
	# 	matches <- remove_sps_u[is.na(match(remove_sps_u, unique(justfirst)))]
	
	# 	unlist(lapply(split(remove_onechar_epithet, as.factor(remove_onechar_epithet)), makenull), use.names=F)
	
	return(dupsout)
}

#' Function to check names on a species list or phylogeny.
#' 
#' @import taxize
#' @param phylo Phylogeny object, class phylo
#' @param source_ Source to match names to. One of ncbi, iPlant
#' @param splitby Length of chunks by which to split species list
#' @examples \dontrun{
#' # A phylogeny as input
#' library(doMC)
#' out <- checknames(phylo="BergmannEtAl2012", source_="NCBI", splitby=200)
#' out <- checknames(phylo=A"lfaroEtAl2009_tree", source_="NCBI", splitby=100)
#' out <- checknames(phylo="EastmanEtAlUnpublished_tree", source_="NCBI", splitby=100)
#' out <- checknames(phylo="HeathEtAl2012_tree", source_="NCBI")
#' out <- checknames(phylo="JaffeEtAl2011", source_="NCBI", splitby=100)
#' out <- checknames(phylo="Chaetodontidae2011", source_="NCBI")
#' out <- checknames(phylo="Drosophila2012_large", source_="NCBI", splitby=50, writefile=TRUE)
#' out <- checknames(phylo="Drosophila2012_small", source_="NCBI", splitby=50)
#' 
#' out <- checknames(phylo="Apogonidae2011", source_="NCBI")
#' out <- checknames(phylo="BinindaEmondsEtAl2007", source_="NCBI", splitby=500)
#' out <- checknames(phylo="HardyCook2012", source_="NCBI", writefile=TRUE)
#' out <- checknames(phylo="FabreEtAl2009", source_="NCBI", splitby=50)
#' out <- checknames(phylo="Pyron2011", source_="NCBI")
#' out <- checknames(phylo="SantiniEtAl2009", source_="NCBI", splitby=50)
#' 
#' # A character vector as input (of species names)
#' out <- checknames(charvector=AlfaroEtAl2009_tree$tip.label, source_="NCBI", splitby=50)
#' out <- checknames(charvector=AlfaroEtAl2009_tree$tip.label[1:20], source_="NCBI")
#' }
checknames <- function(phylo=NULL, charvector=NULL, source_ = "NCBI", 
	splitby = NULL, writefile = FALSE, byfilename = TRUE)
{	
	if(!is.null(phylo)){
		if(byfilename){			
			obj <- eval(parse( text=phylo ))
		} else { obj <- phylo }
		
		# if multiphylo compress the multiphylo so just one set of tip labels
		if(class(obj)=="multiPhylo") {
			obj <- .compressTipLabel(obj)
			orig_tips <- attr(obj,"TipLabel")
			tree_tips <- str_replace_all(orig_tips, "_", " ") # replace underscores
			message(paste("Cleaning names on ", phylo, " - ", length(tree_tips), " taxa", sep=""))
		} else
		{
			orig_tips <- obj$tip.label	
			tree_tips <- str_replace_all(obj$tip.label, "_", " ") # replace underscores
			message(paste("Cleaning names on ", phylo, " - ", length(tree_tips), " taxa", sep=""))
		}
	} else
		if(!is.null(charvector)) {
			tree_tips <- str_replace_all(charvector, "_", " ") # replace underscores
			orig_tips <- tree_tips
			message(paste("Cleaning names on your vector ", " - ", length(tree_tips), " taxa", sep=""))
		} else stop("please provide an object in either phylo or charvector")
	
	# Run tnrs function from taxize to clean names
	if(is.null(splitby)){
		out <- tnrs(tree_tips, getpost="POST", source_=source_) # get TNRS data
		if(nrow(out)==0){ # if no rows in data.frame return X
			out <- data.frame(submittedName=tree_tips, 
												acceptedName=tree_tips, 
												sourceId=rep("nomatch",length(tree_tips)), 
												score=rep(0,length(tree_tips)))
		} else
		{
			out <- out[,1:4]
		}
	} else
	{
		registerDoMC(cores=4)
		out <- ldply(slice(tree_tips, splitby), function(x) tnrs(x, getpost="POST", source_=source_)[,1:4], .parallel=TRUE) # get TNRS data
		if(nrow(out)==0){ # if no rows in data.frame return X
			out <- data.frame(submittedName=tree_tips, 
												acceptedName=tree_tips, 
												sourceId=rep("nomatch",length(tree_tips)), 
												score=rep(0,length(tree_tips)))
		} else
		{
			out <- out[,1:4]
		}
	}
	
# 	temp <- out[out$sourceId %in% source_, ] # get data for user specified source only
	temp <- out[!duplicated(out),] # remove duplicates, if any
	
	names_match <- function(x){  # function to grab submitted name if matched or accepted if not
		if(x$submittedName %in% x$acceptedName){
			return(as.character(x$submittedName))
		} else
		{
			return(as.character(x$acceptedName))
		}
	} 
	registerDoMC(4)
	temp2 <- ddply(temp, .(submittedName), names_match, .parallel=TRUE)
	
	nosourcematch <- unique(out$submittedName[!out$submittedName %in% temp$submittedName]) # no match to source
	temp22 <- rbind(temp2, data.frame(submittedName=nosourcematch, V1=nosourcematch)) # add in species for which there was no match
	
	# replace spaces with underscores in both columns
	temp22$submittedName <- str_replace_all(temp22$submittedName, " ", "_")
	temp22$V1 <- str_replace_all(temp22$V1, " ", "_")
	
	notnrsmatch <- orig_tips[!orig_tips %in% as.character(temp22$submittedName)] # no match at all
	temp33 <- rbind(temp22, data.frame(submittedName=notnrsmatch, V1=notnrsmatch)) # add notnrsmatches to data.frame
	order_ <- sapply(orig_tips, function(x) match(x, temp33$submittedName), USE.NAMES=F)
	temp3 <- temp33[order_,] # reorder data.frame to order of tip.labels
	
	# Keep names that match genera, but not those that don't match genera
	names_match_bin <- function(x){ 
		if(as.character(x$old) %in% as.character(x$new)) { 0 } else { 1 }
	}
	df_ <- data.frame(old=orig_tips, new=temp3$V1)
	df_2 <- ddply(df_, .(old), names_match_bin)
	df_2_merge <- merge(df_, df_2, by="old")
	df_3 <- droplevels(df_2_merge[df_2_merge$V1 == 1, ])
	
	genus_match <- function(x){
		one <- str_split(x$old, "_")[[1]]
		two <- str_split(x$new, "_")[[1]]
		ifelse(one[[1]] %in% two[[1]], as.character(x$new), as.character(x$old))
	}
	asdf <- adply(df_3, 1, genus_match)
	
	if(!nrow(asdf)==0){
		order_2 <- sapply(temp3[temp3$submittedName %in% asdf$old, "submittedName"], function(x) match(x, as.character(asdf$old)), USE.NAMES=F)
		asdf2 <- asdf[order_2,]
		temp3[temp3$submittedName %in% asdf$old, "V1"] <- asdf2$V1
	} else
	{
		temp3 <- temp3
	}	
	
	# clean names, remove unid. species, etc., 
	# when taxa match after cleaning, keep only the one taxon with longest branch
	tips_ <- detect_unid(temp3$V1)
	
	# replace spaces with underscores
	tips_2 <- str_replace_all(tips_, " ", "_")
	
	if(class(obj)=="phylo"){
		obj$tip.label <- tips_2 # assign new names to tip.labels on phylo object
		temppp <- drop.tip(obj, tip="NULL NULL") # drop tips with "NULL NULL"
		
		if(writefile){
			trees_cleaned <- temppp
			save(trees_cleaned, file=paste(phylo,"_cleaned.rda",sep=""))
		} else
		{ return(temppp) }
	} else
		if(class(obj)=="multiPhylo"){
			attr(obj,"TipLabel") <- tips_2
			obj <- .uncompressTipLabel(obj)
			
			temppp <- lapply(obj, function(x) drop.tip(x, tip="NULL NULL") )
			class(temppp) <- "multiPhylo"
			
			if(writefile){
				trees_cleaned <- temppp
				save(trees_cleaned, file=paste("/Library/WebServer/Sites/datelife.org/datelife/data/",phylo,".rda",sep=""))
# 				save(trees_cleaned, file=paste(phylo,"_cleaned.rda",sep=""))
			} else
			{ return(temppp) }
		} else
		{
			charvector <- tips_2
			return(charvector)
		}
}