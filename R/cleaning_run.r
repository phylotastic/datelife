# Check species names against Taxosaurus
# Run and write trees to directory with new file name, just appending "_new" to the end
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

# taxonomic group assignment for trees
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
	list("Oaks2012", "plants", "iPlant_TNRS"),
	list("Pyron2011", "amphibians", "NCBI"),
	list("PyronWiens2011", "amphibians", "NCBI"),
	list("SantiniEtAl2009", "fishes", "NCBI"),
	list("Unpub1", "amphibians", "NCBI"),
	list("Unpub2", "fishes", "NCBI"),
	list("Unpub3", "verts", "NCBI"),
	list("ZhangandWake2009", "amphibians", "NCBI")
)

# Wrapper to checknames to be able to specify taxonomic group to speed up function
check_wrapper <- function(list_, ...){
	checknames_safe(phylo=list_[[1]], source_=list_[[3]], ...)
}

# Run all, with the same arguments to TNRS, change arguments per tree in the 
# future based on taxonomic group
checknames_safe <- plyr::failwith(NULL, checknames_safe)
l_ply(treegrouplist[1:2], check_wrapper, splitby=100, writefile=FALSE) # run the cleaning algorithm on all trees