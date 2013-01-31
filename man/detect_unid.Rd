\name{detect_unid}
\alias{detect_unid}
\title{Function to detect messy names and remove them, and remove sp.'s}
\usage{
  detect_unid(taxon)
}
\arguments{
  \item{taxon}{Taxonomic name, or a vector of taxonomic
  names.}
}
\description{
  Function to detect messy names and remove them, and
  remove sp.'s
}
\examples{
\dontrun{
detect_unid(taxon = Drosophila2012_large[[1]]$tip.label)
}
}

