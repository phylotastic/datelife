# Including TNRS in Datelife

#' Find only rows in a data.frame where the values don't match in a row
#' 
#' @pararm x A data.frame
#' @examples
#' dff <- data.frame(FabreEtAl2009[[1]][[1]]$tip.label, out$tip.label)
#' nonmatching_by_row(treename=dff)
#' 
nonmatching_by_row <- function(treename, cleanedtree){
	dff <- data.frame(eval(pasrse(text=treename))$tip.label, cleanedtree$tip.label)
	foo <- function(y){ ifelse(y[,1] %in% y[,2], 0, 1) }
	temp <- adply(dff, 1, foo)
	temp[temp$V1 == 1,]
}

#' Function to break up a vector by certain N sized chunks.
#' 
#' @import plyr
#' @param input Input character vector
#' @pararm by Number by which to split the input character vector
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
#' detect_unid(taxon = Drosophila2012_large[[1]]$tip.label)
#' }
detect_unid <- function(taxon){
	# replace underscores
	taxon2 <- str_replace_all(taxon, "_", " ") # replace underscores
	# detect duplicates
	dups <- sapply(taxon2, function(x) paste(str_split(x, " ")[[1]][1:2], collapse=" ",sep=""), USE.NAMES=F)
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
#' out <- checknames(phylo=Apogonidae2011, source_="NCBI")
#' out <- checknames(phylo="BinindaEmondsEtAl2007", source_="NCBI", splitby=500)
#' out <- checknames(phylo=FabreEtAl2009, source_="NCBI", splitby=50)
#' out <- checknames(phylo=HardyCook2012, source_="NCBI", writefile=TRUE)
#' out <- checknames(phylo="FabreEtAl2009", source_="NCBI", splitby=50, writefile=TRUE)
#' 
#' # A character vector as input (of species names)
#' out <- checknames(charvector=AlfaroEtAl2009_tree$tip.label, source_="NCBI", splitby=50)
#' out <- checknames(charvector=AlfaroEtAl2009_tree$tip.label[1:20], source_="NCBI")
#' }
checknames <- function(phylo=NULL, charvector=NULL, source_ = "NCBI", splitby = NULL, writefile = FALSE)
{	
	if(!is.null(phylo)){
		obj <- eval(parse( text=phylo ))
		
		# if multiphylo compress the multiphylo so just one set of tip labels
		if(class(obj)=="multiPhylo") {
			obj <- .compressTipLabel(obj)
			orig_tips <- attr(obj,"TipLabel")
			tree_tips <- str_replace_all(orig_tips, "_", " ") # replace underscores
			message(paste("Cleaning names on ", phylo, "- ", length(tree_tips), " taxa", sep=""))
		} else
		{
			orig_tips <- obj$tip.label	
			tree_tips <- str_replace_all(obj$tip.label, "_", " ") # replace underscores
			message(paste("Cleaning names on ", phylo, "- ", length(tree_tips), " taxa", sep=""))
		}
	} else
	if(!is.null(charvector)) {
		tree_tips <- str_replace_all(charvector, "_", " ") # replace underscores
		orig_tips <- tree_tips
		message(paste("Cleaning names on your vector ", "- ", length(tree_tips), " taxa", sep=""))
	} else stop("please provide an object in either phylo or charvector")
	
	
# 	tree_tips <- detect_unid(tree_tips) # clean names, remove unid. species, etc.
	
	if(is.null(splitby)){
		out <- tnrs(tree_tips, getpost="POST")[,1:4] # get TNRS data
	} else
	{
		registerDoMC(cores=4)
		out <- ldply(slice(tree_tips, splitby), function(x) tnrs(x, getpost="POST")[,1:4], .parallel=TRUE) # get TNRS data
	}
	
	temp <- out[out$sourceId %in% source_, ] # get data for user specified source only
	temp <- temp[!duplicated(temp),] # remove duplicates, if any
	
	foo2 <- function(x){  # function to grab submitted name if matched or accepted if not
		if(x$submittedName %in% x$acceptedName){
			return(as.character(x$submittedName))
		} else
		{
			return(as.character(x$acceptedName))
		}
	} 
	registerDoMC(4)
	temp2 <- ddply(temp, .(submittedName), foo2, .parallel=TRUE)
	
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
	foooo <- function(x){ 
		if(as.character(x$old) %in% as.character(x$new)) { 0 } else { 1 }
	}
	df_ <- data.frame(old=orig_tips, new=temp3$V1)
	df_2 <- ddply(df_, .(old), foooo)
	df_2_merge <- merge(df_, df_2, by="old")
	df_3 <- droplevels(df_2_merge[df_2_merge$V1 == 1, ])
	
	genusmatch <- function(x){
		one <- str_split(x$old, "_")[[1]]
		two <- str_split(x$new, "_")[[1]]
		ifelse(one[[1]] %in% two[[1]], as.character(x$new), as.character(x$old))
	}
	asdf <- adply(df_3, 1, genusmatch) 
	
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
	
	if(class(obj)=="phylo"){
		obj$tip.label <- tips_ # assign new names to tip.labels on phylo object
		temppp <- drop.tip(obj, tip="NULL NULL") # drop tips with "NULL NULL"
		
		if(writefile){
			trees_cleaned <- temppp
			save(trees_cleaned, file=paste(phylo,"_cleaned.rda",sep=""))
		} else
			{ return(obj) }
	} else
		if(class(obj)=="multiPhylo"){
			attr(obj,"TipLabel") <- tips_
			obj <- .uncompressTipLabel(obj)
			
			temppp <- lapply(obj, function(x) drop.tip(x, tip="NULL NULL") )
			class(temppp) <- "multiPhylo"
			
			if(writefile){
				trees_cleaned <- temppp
				save(trees_cleaned, file=paste(phylo,"_cleaned.rda",sep=""))
			} else
			{ return(temppp) }
		} else
	{
		charvector <- tips_
		return(charvector)
	}
}

# Check species names against Taxosaurus
# Run and write trees to directory with new file name, just appending "_new" to the end
# Run the function replacenames across all trees
# install_github("taxize_", "ropensci", ref="release")
library(taxize); library(ape); library(stringr)
# treefilenames <- c("AlfaroEtAl2009.rda","BergmannEtAl2012.rda","BinindaEmondsEtAl2007.rda",
# 		"EastmanEtAlUnpublished.rda","HeathEtAl2012.rda","JaffeEtAl2011.rda","megaTree.rda",
# 		"Oaks2011.rda","PyronWiens2011.rda","RaboskyEtAlUnpublished.rda","SantiniEtAl2009.rda",
# 		"SmithEtAl2011.rda","TreeBase.rda","ZhangandWake2009.rda")
treefilenames <- dir("/Users/scottmac2/phyloorchard/pkg/data")

# Load all trees into the workspace
# CHANGE PATH!!!!
l_ply(treefilenames, function(x) load(paste("/Users/scottmac2/phyloorchard/pkg/data/",x,sep=""), .GlobalEnv))

# Get tree names in the workspace
trees <- sapply(treefilenames, function(x) str_replace(x, ".rda", ""), USE.NAMES=F)
trees <- trees[-c(1,10)]

# Run all, with the same arguments to TNRS, change arguments per tree in the 
# future based on taxonomic group
checknames_safe <- plyr::failwith(NULL, checknames)
l_ply(trees, checknames_safe, source_="NCBI", splitby=100, writefile=TRUE)

# Check names that didn't match
# nonmatching_by_row()