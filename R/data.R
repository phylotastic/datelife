# #' Description of data_object
# #'
# #' @name data_object
# #' @docType data OR package OR method OR class OR anything else
# #' @author Luna L. Sanchez-Reyes \email{lsanche7@utk.edu}
# #' Brian O'Meara \email{bomeara@utk.edu}
# #' @references \url{https://data_references.org}
# #' @source \url{http://data_source.org}
# #' @format A list
# #' \describe{
# #'   \item{name1}{Character vector of }
# #'   \item{name2}{Numeric vector of }
# #'   \item{name3}{Integer with }
# #' }
# #' @keywords
# #' @details
# #' How was this data obtained or generated
# "data_object"

#' Chronograms in Open Tree of Life and other related data
#'
#' Now storing >800 chronograms from OToL
#'
#' @name opentree_chronograms
#' @docType data
#' @format A list of four elements, containing data on OToL chronograms
#' \describe{
#'   \item{authors}{List of lists of authors of the included studies}
#'   \item{curators}{List of lists of curators of the included studies}
#'   \item{studies}{List of study identifiers}
#'   \item{trees}{List storing the chronograms from OpenTree}
#' }
#' @source \url{http://opentreeoflife.org}
#' @keywords otol tree chronogram
#' @details
#' Generated with update_datelife_cache() via get_otol_chronograms()
"opentree_chronograms"


#' Authors from studies with chronograms in Open tree of Life
#'
#' @name author.pretty
#' @docType data
#' @format A character vector with the author names from studies with chronograms that are in OToL
#' @source \url{http://opentreeoflife.org}
#' @keywords otol study studies tree chronogram author
#' @details
#' Generated with make_contributor_cache()
"author.pretty"


#' Information on authors, study ids and clades from chronograms in Open Tree of Life
#'
#'
#' @name author.results
#' @docType data
#' @format A dataframe with three variables
#' \describe{
#'   \item{person}{Character vector of authors of the included studies}
#'   \item{study}{Character vector of study identifiers}
#'   \item{clade}{Character vector with names of clades with chronograms from the included studies}
#' }
#' @source \url{http://opentreeoflife.org}
#' @keywords otol study studies tree chronogram id author clade
#' @details
#' Generated with make_contributor_cache()
"author.results"


#' Curators fof chronograms in Open tree of Life
#'
#' Curators are the ones uploading the chronograms to Open Tree of Life
#'
#' @name curator.pretty
#' @docType data
#' @format A character vector with the names of curators of chronograms that are in OToL
#' @source \url{http://opentreeoflife.org}
#' @keywords otol study studies tree chronogram curator
#' @details
#' Generated with make_contributor_cache()
"curator.pretty"

#' Information on curators, study ids and clades from chronograms in Open Tree of Life
#'
#' Curators are the ones uploading the chronograms to Open Tree of Life
#'
#' @name curator.results
#' @docType data
#' @format A dataframe with three variables
#' \describe{
#'   \item{person}{Character vector of curators of the included studies}
#'   \item{study}{Character vector of study identifiers}
#'   \item{clade}{Character vector with names of clades with chronograms from the included studies}
#' }
#' @source \url{http://opentreeoflife.org}
#' @keywords otol study studies tree chronogram curator id clade
#' @details
#' Generated with make_contributor_cache()
"curator.results"


#' Study ids of chronograms in Open Tree of Life not included in authors.
#'
#' There are different reasons why a doi cannot be found online, so author data
#' cannot be retrieved for that particular study.
#' This vector contains all the studies with chrongrams in OToL whose information
#' was not included automatically in author.pretty and author.results
#' @name missed_doi
#' @docType data
#' @format A character vector with study ids
#' @source \url{http://opentreeoflife.org}
#' @keywords otol study studies tree chronogram id
#' @details
#' Generated with make_contributor_cache()
"missed_doi"


#' Some plants chronogram
#'
#'
#' @name plant_bold_otol_tree
#' @docType data
#' @format A phylo object with 6 tips and 5 internal nodes
#' \describe{
#'   \item{edge}{Integer vector with edge (branch) numbers}
#'   \item{tip.label}{Character vector with species names of plants}
#'   \item{Nnode}{Integer vector with the number of nodes}
#'   \item{node.label}{Character vector wit node names}
#'   \item{edge.length}{Numeric vector with edge (branch) lengths}
#' }
#' @author Luna L. Sanchez-Reyes \email{lsanche7@utk.edu}
#' @author Brian O'Meara \email{bomeara@utk.edu}
#' @source \url{http://opentreeoflife.org}
#' @source \url{http://www.boldsystems.org}
#' @keywords otol plant tree chronogram bold
#' @details
#' Generated with make_bold_otol_tree(input = "((Zea mays,Oryza sativa),((Arabidopsis thaliana,(Glycine max,Medicago sativa)),Solanum lycopersicum)Pentapetalae);")
"plant_bold_otol_tree"


#' Authors from studies with chronograms in TreeBASE
#'
#'
#' @name tb.author.pretty
#' @docType data
#' @format A data.frame with two elements
#' \describe{
#'   \item{person}{Character vector of authors of the included studies}
#'   \item{study.count}{Numeric vector of studies authored by each person}
#' }
#' @source \url{http://treebase.org}
#' @keywords treebase tree chronogram author study count
#' @details
#' Generated with make_treebase_cache()
"tb.author.pretty"


#' Information on authors and study ids from chronograms in TreeBASE
#'
#'
#' @name tb.author.results
#' @docType data
#' @format A dataframe with three variables
#' \describe{
#'   \item{person}{Named character vector of authors of the included studies}
#'   \item{study}{Named character vector of study identifiers}
#' }
#' @source \url{http://treebase.org}
#' @keywords treebase tree chronogram study studies id
#' @details
#' Generated with make_treebase_cache()
"tb.author.results"
