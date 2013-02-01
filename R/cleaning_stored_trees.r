# Check species names against Taxosaurus
# Run the function replacenames across all trees
library(taxize); library(ape); library(stringr); library(doMC)

# treefilenames <- dir("/Users/scottmac2/phyloorchard/pkg/data")
treefilenames <- dir("/Library/WebServer/Sites/datelife.org/datelife/data")

# Load all trees into the workspace
# CHANGE PATH!!!!
# l_ply(treefilenames, function(x) load(paste("/Users/scottmac2/phyloorchard/pkg/data/",x,sep=""), .GlobalEnv))
l_ply(treefilenames, function(x) load(paste("/Library/WebServer/Sites/datelife.org/datelife/data/",x,sep=""), .GlobalEnv))

# Get tree names in the workspace
trees <- sapply(treefilenames, function(x) str_replace(x, ".rda", ""), USE.NAMES=F)

# Taxonomic group assignment for trees  - Add to this as needed
treegrouplist <- list(
	list("AlfaroEtAl2009", "verts", "NCBI"),
	list("Apogonidae2011", "fishes", "NCBI"),
	list("BergmannEtAl2012", "reptiles", "NCBI"),
	list("BinindaEmondsEtAl2007", "mammals", "MSW3"),
	list("Conifers2012", "plants", "iPlant_TNRS"),
	list("Chaetodontidae2011", "fishes", "NCBI"),
	list("Labridae2011", "fishes", "NCBI"),
	list("Pomacentridae2011", "fishes", "NCBI"),
	list("Drosophila2012_large", "flies", "NCBI"),
	list("Drosophila2012_small", "flies", "NCBI"),
	list("FabreEtAl2009", "mammals", "MSW3"),
	list("HardyCook2012", "plants", "iPlant_TNRS"),
	list("HeathEtAl2012", "mammals", "MSW3"),
	list("JaffeEtAl2011", "reptiles", "NCBI"),
	list("LohmannEtAl2013", "plants", "iPlant_TNRS"),
	list("LopezFernandezEtAl2013", "fishes", "NCBI"),
	list("NauheimerEtAl2012", "plants", "iPlant_TNRS"),
	list("Oaks2011", "crocodiles", "NCBI"),
	list("Pyron2011", "amphibians", "NCBI"),
	list("PyronWiens2011", "amphibians", "NCBI"),
	list("SantiniEtAl2009", "fishes", "NCBI"),
	list("Unpub1", "amphibians", "NCBI"),
	list("Unpub2", "fishes", "NCBI"),
	list("Unpub3", "verts", "NCBI"),
	list("ZhangandWake2009", "amphibians", "NCBI")
)

# Create a safe version of checknames so that the below lapply doesn't stop if one fails 
checknames_safe <- plyr::failwith(NULL, checknames)

# Wrapper to checknames to be able to specify taxonomic group to speed up function
check_wrapper <- function(list_, ...){
	checknames_safe(phylo=list_[[1]], source_=list_[[3]], ...)
}

# Run all, with the same arguments to TNRS
l_ply(treegrouplist[1:2], check_wrapper, splitby=200, writefile=TRUE, writedir="/Users/schamber/datelife_testout/")
# l_ply(treegrouplist[1:2], check_wrapper, splitby=200, writefile=TRUE)