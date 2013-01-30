\name{checknames}
\alias{checknames}
\title{Function to check names on a species list or phylogeny.}
\usage{
  checknames(phylo = NULL, charvector = NULL,
    source_ = "NCBI", splitby = NULL, writefile = FALSE)
}
\arguments{
  \item{phylo}{Phylogeny object, class phylo}

  \item{source_}{Source to match names to. One of ncbi,
  iPlant}

  \item{splitby}{Length of chunks by which to split species
  list}
}
\description{
  Function to check names on a species list or phylogeny.
}
\examples{
\dontrun{
# A phylogeny as input
library(doMC)
out <- checknames(phylo="BergmannEtAl2012", source_="NCBI", splitby=200)
out <- checknames(phylo=A"lfaroEtAl2009_tree", source_="NCBI", splitby=100)
out <- checknames(phylo="EastmanEtAlUnpublished_tree", source_="NCBI", splitby=100)
out <- checknames(phylo="HeathEtAl2012_tree", source_="NCBI")
out <- checknames(phylo="JaffeEtAl2011", source_="NCBI", splitby=100)
out <- checknames(phylo="Chaetodontidae2011", source_="NCBI")
out <- checknames(phylo="Drosophila2012_large", source_="NCBI", splitby=50, writefile=TRUE)
out <- checknames(phylo="Drosophila2012_small", source_="NCBI", splitby=50)

out <- checknames(phylo=Apogonidae2011, source_="NCBI")
out <- checknames(phylo="BinindaEmondsEtAl2007", source_="NCBI", splitby=500)
out <- checknames(phylo=FabreEtAl2009, source_="NCBI", splitby=50)
out <- checknames(phylo=HardyCook2012, source_="NCBI", writefile=TRUE)
out <- checknames(phylo="FabreEtAl2009", source_="NCBI", splitby=50, writefile=TRUE)

# A character vector as input (of species names)
out <- checknames(charvector=AlfaroEtAl2009_tree$tip.label, source_="NCBI", splitby=50)
out <- checknames(charvector=AlfaroEtAl2009_tree$tip.label[1:20], source_="NCBI")
}
}

