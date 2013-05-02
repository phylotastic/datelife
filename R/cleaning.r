# Including TNRS in Datelife

#' Find only rows in a data.frame where the values don't match in a row
#' 
#' @import plyr
#' @param df A data.frame
#' @examples
#' dff <- data.frame(FabreEtAl2009$tip.label, out$tip.label)
#' nonmatching_by_row(dff)
nonmatching_by_row <- function(df){
	does_match <- function(y){ ifelse(y[,1] %in% y[,2], 0, 1) }
# 	df2 <- data.frame(old=sort(df[,1]), new=sort(df[,2]))
	temp <- adply(df, 1, does_match)
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
slice <- function(input, by = 2, equalpiecesof = NULL) {
	if(is.null(equalpiecesof)){	
		starts <- seq(1, length(input), by)
		tt <- lapply(starts, function(y) input[y:(y + (by - 1))])
		llply(tt, function(x) x[!is.na(x)])
	} else
	{
		splitby <- round(length(input)/equalpiecesof)+1
		starts <- seq(1, length(input), splitby)
		tt <- lapply(starts, function(y) input[y:(y + (splitby - 1))])
		llply(tt, function(x) x[!is.na(x)])
	}
}

#' Function to detect messy names and remove them, and remove sp.'s
#' 
#' @import stringr
#' @param taxon Taxonomic name, or a vector of taxonomic names.
#' @examples \dontrun{
#' drop_epithets(taxon = Oaks2011$tip.label)
#' }
drop_epithets <- function(taxon){
	# replace underscores
	taxon2 <- str_replace_all(taxon, "_", " ")
	
	# remove sp.'s, etc. 
	remove_sps <- gsub(" sp[^aeiouh]+| sp\\.+| sp$", "", taxon2)
	
	# remove epithets with just one letter
	remove_onechar_epithet <- gsub(" [A-Za-z0-9]{1}$", "", remove_sps)

	# put underscores back in
# 	return(str_replace_all(remove_onechar_epithet, " ", "_"))
	return(remove_onechar_epithet)
}

#' Function to check names on a species list or phylogeny.
#' 
#' @import taxize
#' @param phylo Phylogeny object, class phylo
#' @param charvector Character vector.
#' @param source_ Source to match names to. One of ncbi, iPlant
#' @param splitby Length of chunks by which to split species list
#' @param writefile Write file to directory or not. Remember to change directory to 
#' 		write to in the function. Right now it is set to 
#' 		"/Library/WebServer/Sites/datelife.org/datelife/data/"
#' @param writedir A directory to write file to (defaults to "~/"); ignore if writefile = FALSE.
#' @param byfilename Defaults to TRUE - if TRUE, then interprets the phylo argument input 
#' 		using eval() so executes it (used in name fixing on stored trees). If FALSE, 
#' 		then you can specify phylo as a phylo object instead of the name of that object.
#' @examples \dontrun{
#' # Cleaning stored trees
#' # A phylogeny as input
#' library(doMC); library(ape)
#' out <- checknames(phylo="BergmannEtAl2012", source_="NCBI", splitby=200)
#' out <- checknames(phylo="AlfaroEtAl2009", source_="NCBI", splitby=200)
#' out <- checknames(phylo="EastmanEtAlUnpublished_tree", source_="NCBI", splitby=100)
#' out <- checknames(phylo="HeathEtAl2012_tree", source_="NCBI")
#' out <- checknames(phylo="JaffeEtAl2011", source_="NCBI", splitby=100)
#' out <- checknames(phylo="Chaetodontidae2011", source_="NCBI")
#' out <- checknames(phylo="Drosophila2012_large", source_="NCBI", splitby=50, writefile=TRUE)
#' out <- checknames(phylo="Drosophila2012_small", source_="NCBI", splitby=50)
#' 
#' out <- checknames(phylo="Apogonidae2011", source_="NCBI")
#' out <- checknames(phylo="BinindaEmondsEtAl2007", source_="NCBI", splitby=500)
#' out <- checknames(phylo="HardyCook2012", source_="NCBI")
#' out <- checknames(phylo="FabreEtAl2009", source_="NCBI", splitby=50)
#' out <- checknames(phylo="PyronWiens2011", source_="NCBI", splitby=500, writefile=TRUE)
#' out <- checknames(phylo="Unpub3", source_="NCBI", splitby=50, writefile=TRUE)
#' 
#' # Cleaning user input trees
#' # A character vector as input (of species names)
#' out <- checknames(charvector=AlfaroEtAl2009$tip.label[1:10], source_="NCBI")
#' 
#' # A phylo object - user suppressMessages() around the fxn
#' mmm<-drop.tip(AlfaroEtAl2009, 1:190)
#' out <- suppressMessages(checknames(phylo=phy, source_="NCBI", byfilename=FALSE))
#' }
checknames <- function(phylo=NULL, charvector=NULL, source_ = "NCBI",
	writefile = FALSE, writedir="~/", byfilename = TRUE)
{	
	if(!is.null(phylo)){
		if(byfilename){
			obj <- eval(parse( text=phylo ))
		} else { obj <- phylo }
		
		# if multiphylo compress the multiphylo so just one set of tip labels
		if(class(obj)=="multiPhylo") {
			# Sample max of 100 trees from a multiPhylo object
			max.trees.per.study <- 100 #some studies will have a huge number of trees. They will be consistently sampled up to this number
			obj <- obj[unique(round(seq(from=1, to=length(obj), length.out=max.trees.per.study)))]
			
			# Compress tip labels
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
			obj <- NULL
			tree_tips <- str_replace_all(charvector, "_", " ") # replace underscores
			orig_tips <- tree_tips
			message(paste("Cleaning names on your vector ", " - ", length(tree_tips), " taxa", sep=""))
		} else stop("please provide an object in either phylo or charvector")
	
	# Run tnrs function from taxize to clean names
	if(length(tree_tips) < 500){
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
		registerDoMC(cores=8)
		out <- ldply(slice(tree_tips, equalpiecesof=8), function(x) tnrs(x, getpost="POST", source_=source_)[,1:4], .parallel=TRUE) # get TNRS data
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
	registerDoMC(8)
	temp2 <- ddply(temp, .(submittedName), names_match, .parallel=TRUE)
	
	nosourcematch <- unique(out$submittedName[!out$submittedName %in% temp$submittedName]) # no match to source
	temp22 <- rbind(temp2, data.frame(submittedName=nosourcematch, V1=nosourcematch)) # add in species for which there was no match
	
	# replace spaces with underscores in both columns
# 	temp22$submittedName <- str_replace_all(temp22$submittedName, " ", "_")
# 	temp22$V1 <- str_replace_all(temp22$V1, " ", "_")
	
	notnrsmatch <- orig_tips[!orig_tips %in% as.character(temp22$submittedName)] # no match at all
  
	# replace spaces with underscores in both columns
	temp22$submittedName <- str_replace_all(temp22$submittedName, " ", "_")
	temp22$V1 <- str_replace_all(temp22$V1, " ", "_")
  
	temp33 <- rbind(temp22, data.frame(submittedName=notnrsmatch, V1=notnrsmatch)) # add notnrsmatches to data.frame
	order_ <- sapply(str_replace(orig_tips," ","_"), function(x) match(x, temp33$submittedName), USE.NAMES=F)
	temp3 <- temp33[order_,] # reorder data.frame to order of tip.labels
	
	# Keep names that match genera, but not those that don't match genera
	names_match_bin <- function(x){ 
		if(as.character(x$old) %in% as.character(x$new)) { 0 } else { 1 }
	}
	df_ <- data.frame(old=orig_tips, new=temp3$V1)
# 	df_2 <- ddply(df_, .(old), names_match_bin)
	df_2 <- adply(df_, 1, names_match_bin)[,-2]
	df_2_merge <- merge(df_, df_2, by="old")
	df_3 <- droplevels(df_2_merge[df_2_merge$V1 == 1, ])
	
	genus_match <- function(x){
	  one <- str_split(str_replace(x$old, " ", "_"), "_")[[1]][[1]]
	  two <- str_split(x$new, "_")[[1]][[1]]
	  ifelse(one[[1]] %in% two[[1]], as.character(x$new), as.character(x$old))
	}
	asdf <- adply(df_3, 1, genus_match)
	
	if(!nrow(asdf)==0){
	  asdf$old <- str_replace_all(asdf$old, " ", "_")
		order_2 <- sapply(temp3[temp3$submittedName %in% asdf$old, "submittedName"], function(x) match(x, as.character(asdf$old)), USE.NAMES=F)
		asdf2 <- asdf[order_2,]
		temp3[temp3$submittedName %in% asdf$old, "V1"] <- asdf2$V1
	} else
	{
		temp3 <- temp3
	}	
	
	# clean names, remove unid. species, etc., 
	# when taxa match after cleaning, keep only the one taxon with longest branch
	# COMMENTED OUT FOR NOW, NEED TO WORK OUT SOME ISSUES WITH CLEANING NAMES - 
	# 	DO WE DROP SPECIES WITH EPTITHETS LIKE "sp."
# 	tips_ <- detect_unid(temp3$V1)
	
	# replace spaces with underscores
	tips_2 <- str_replace_all(temp3$V1, " ", "_")
	
	if(class(obj)=="phylo"){
		obj$tip.label <- tips_2 # assign new names to tip.labels on phylo object
# 		temppp <- drop.tip(obj, tip="NULL NULL") # drop tips with "NULL NULL"
		
		if(writefile){
			trees_cleaned <- obj
			save(trees_cleaned, file=paste(phylo,".rda",sep=""))
		} else
		{ return(obj) }
	} else
		if(class(obj)=="multiPhylo"){
			attr(obj,"TipLabel") <- tips_2
			obj <- .uncompressTipLabel(obj)
			
# 			temppp <- lapply(obj, function(x) drop.tip(x, tip="NULL NULL") )
			class(obj) <- "multiPhylo"
			
			if(writefile){
				trees_cleaned <- obj
				save(trees_cleaned, file=paste(writedir,phylo,".rda",sep=""))
# 				save(trees_cleaned, file=paste(phylo,"_cleaned.rda",sep=""))
			} else
			{ return(obj) }
		} else
		{
			return(tips_2)
		}
}

#' Test timing of name cleaning on datelife trees
#' 
#' @import plyr ggplot2 taxize stringr
#' @param listoftrees List of trees, each tree in phylo format, not multiPhylo
#' 		format. 
#' @param plotresults Plot results or not.
#' @examples \dontrun{
#' treefilenames <- dir("/Users/scottmac2/phyloorchard/pkg/data")
#' l_ply(treefilenames, function(x) load(paste("/Users/scottmac2/phyloorchard/pkg/data/",x,sep=""), .GlobalEnv))
#' trees <- sapply(treefilenames, function(x) str_replace(x, ".rda", ""), USE.NAMES=F)
#' binemonds <- BinindaEmondsEtAl2007[[1]]
#'  makesmallsubtrees <- function(phylo){
#'  		drops <- c(2000,2500,3000,3500,4000,4100,4200,4300,4350,4400,4420,4450,4470,4500,4505)
#'  		llply(rev(drops), function(x) drop.tip(phylo, tip=1:x))
#' }
#' mytreelist <- makesmallsubtrees(phylo = binemonds)
#' library(ggplot2); library(plyr); library(lubridate)
#' out <- time_cleaning(listoftrees = mytreelist, plotresults=TRUE)
#' out[[1]] # data.frame
#' out[[2]] # plot
#' }
time_cleaning <- function(listoftrees, plotresults=FALSE){
	checknames_safe <- plyr::failwith(NULL, checknames)
	temp <- llply(listoftrees, function(tree){
		start <- Sys.time()
		result <- suppressMessages(checknames_safe(phylo=tree, source_="MSW3", byfilename=FALSE))
		stop_ <- Sys.time()
		elapsed <- as.duration(stop_-start)
		as.numeric(elapsed)
	}, .progress="text"
	)
	
	names(temp) <- sapply(listoftrees, Ntip, USE.NAMES=F)
	resultsout <- ldply(temp)
	names(resultsout) <- c("num_tips","elapsed")
	resultsout$num_tips <- as.numeric(resultsout$num_tips)
	
	if(plotresults){
		p <- ggplot(resultsout, aes(num_tips, elapsed)) + 
			theme_bw(base_size=18) +
			geom_point(size = 4) + 
			labs(y = "Elapsed time", x = "Number of tips")
		list(resultsout, p)
	} else
	{ resultsout }
}