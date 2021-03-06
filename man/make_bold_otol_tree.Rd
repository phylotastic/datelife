% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bold_tree.R
\name{make_bold_otol_tree}
\alias{make_bold_otol_tree}
\title{Use Barcode of Life Data (BOLD) to get branch lengths on a tree topology.
If input is not a tree, then it gets and induced OToL subtree of the set of input taxa.}
\usage{
make_bold_otol_tree(
  input = c("Rhea americana", "Struthio camelus", "Gallus gallus"),
  marker = "COI",
  otol_version = "v3",
  chronogram = TRUE,
  doML = FALSE,
  verbose = FALSE,
  aligner = "muscle",
  ...
)
}
\arguments{
\item{input}{A character vector of taxon names, a tree as a 'phylo' object or a newick character string, or a 'datelifeQuery' object from make_datelife_query function.}

\item{marker}{A character vector indicating the gene from BOLD system to be used for branch length estimation.}

\item{otol_version}{Version of OToL to use}

\item{chronogram}{Boolean. If TRUE (default), branch lengths returned are estimated with ape::chronoMPL. If FALSE, branch lengths returned are estimated with phangorn::acctran and represent relative substitution rates .}

\item{doML}{Boolean; only relevant if chronogram = TRUE. If TRUE, it does ML branch length optimization with phangorn::optim.pml}

\item{verbose}{Boolean. If TRUE, it gives printed updates to the user.}

\item{aligner}{A character vector indicating whether to use MAFFT or MUSCLE to align BOLD sequences. It is not case sensitive. Default to MUSCLE.}

\item{...}{
  Arguments passed on to \code{\link[=get_otol_synthetic_tree]{get_otol_synthetic_tree}}
  \describe{
    \item{\code{resolve}}{boolean default to TRUE. Whether to resolve the tree at random or not.}
    \item{\code{ott_ids}}{If not NULL, it takes this argument and ignores input. A numeric vector of ott ids obtained with rotl::taxonomy_taxon_info or rolt::tnrs_match_names or datelife::tnrs_match}
  }}
}
\value{
A phylogeny with branch lengths proportional to relative substitution rate.
}
\description{
Use Barcode of Life Data (BOLD) to get branch lengths on a tree topology.
If input is not a tree, then it gets and induced OToL subtree of the set of input taxa.
}
\details{
If input is a phylo object or a newick string, it is used as backbone topology.
If input is a character vector of taxon names, an induced OToL tree is used as backbone.
If there are not enough sequences to return a tree with branch lengths, it returns
either the original input topology or the OToL tree obtained for the set of input taxon names.
}
