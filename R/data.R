# #' Description of data_object
# #'
# #' @name data_object
# #' @docType data OR package OR method OR class OR anything else
# #' @author Luna L. Sanchez-Reyes \email{lsanche7@utk.edu}
# #' Brian O'Meara \email{bomeara@utk.edu}
# #' @source \url{http://data_source.org}
# #' @format A list
# #' \describe{
# #'   \item{name1}{Character vector of }
# #'   \item{name2}{Numeric vector of }
# #'   \item{name3}{Integer with }
# #' }
# #' @keywords
# #' @details
# #' How were the data obtained or generated
# "data_object"

#' Information on contributors, authors, study ids and clades from studies with
#' chronograms in Open Tree of Life (Open Tree)
#'
#' @name contributor_cache
#' @docType data
#' @format A list of five data sets.
#' \describe{
#'   \item{author.pretty}{A character vector with the author names from studies
#'    with chronograms that are in Open Tree.}
#'   \item{author.results}{A dataframe with three variables: authors, study ids and clades.}
#'   \item{curator.pretty}{A character vector with the names of curators of
#'    chronograms that are in Open Tree.}
#'   \item{curator.results}{A `data.frame` with three variables: curators, study ids and clades.}
#'   \item{missed_doi}{A character vector with study ids whose "doi" could not be retrieved.}
#' }
#' @source \url{http://opentreeoflife.org}
#' @keywords otol study studies tree chronogram author
#' @details
#' Generated with [make_contributor_cache()].
"contributor_cache"

# #' Authors from studies with chronograms in Open Tree of Life
# #'
# #' @name author.pretty
# #' @docType data
# #' @format A character vector with the author names from studies with chronograms that are in Open Tree
# #' @source \url{http://opentreeoflife.org}
# #' @keywords otol study studies tree chronogram author
# #' @details
# #'
# #' Generated with make_contributor_cache()
# "author.pretty"


# #' Information on authors, study ids and clades from chronograms in Open Tree of Life
# #'
# #'
# #' @name author.results
# #' @docType data
# #' @format A dataframe with three variables
# #' \describe{
# #'   \item{person}{Character vector of authors of the included studies}
# #'   \item{study}{Character vector of study identifiers}
# #'   \item{clade}{Character vector with names of clades with chronograms from the included studies}
# #' }
# #' @source \url{http://opentreeoflife.org}
# #' @keywords otol study studies tree chronogram id author clade
# #' @details
# #' Generated with make_contributor_cache()
# "author.results"


# #' Curators fof chronograms in Open tree of Life
# #'
# #' Curators are the ones uploading the chronograms to Open Tree of Life
# #'
# #' @name curator.pretty
# #' @docType data
# #' @format A character vector with the names of curators of chronograms that are in OToL
# #' @source \url{http://opentreeoflife.org}
# #' @keywords otol study studies tree chronogram curator
# #' @details
# #' Generated with make_contributor_cache()
# "curator.pretty"

# #' Information on curators, study ids and clades from chronograms in Open Tree of Life
# #'
# #' Curators are the ones uploading the chronograms to Open Tree of Life
# #'
# #' @name curator.results
# #' @docType data
# #' @format A dataframe with three variables
# #' \describe{
# #'   \item{person}{Character vector of curators of the included studies}
# #'   \item{study}{Character vector of study identifiers}
# #'   \item{clade}{Character vector with names of clades with chronograms from the included studies}
# #' }
# #' @source \url{http://opentreeoflife.org}
# #' @keywords otol study studies tree chronogram curator id clade
# #' @details
# #' Generated with make_contributor_cache()
# "curator.results"


# #' Study ids of chronograms in Open Tree of Life not included in authors.
# #'
# #' There are different reasons why a doi cannot be found online, so author data
# #' cannot be retrieved for that particular study.
# #' This vector contains all the studies of chrongrams in OToL whose information
# #' was not included automatically in author.pretty and author.results
# #' @name missed_doi
# #' @docType data
# #' @format A character vector with study ids
# #' @source \url{http://opentreeoflife.org}
# #' @keywords otol study studies tree chronogram id
# #' @details
# #' Generated with make_contributor_cache()
# "missed_doi"


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
#'   \item{node.label}{Character vector with node names}
#'   \item{edge.length}{Numeric vector with edge (branch) lengths}
#' }
#' @author Luna L. Sanchez-Reyes \email{lsanche7@utk.edu}
#' @author Brian O'Meara \email{bomeara@utk.edu}
#' @source \url{http://opentreeoflife.org}
#' @source \url{http://www.boldsystems.org}
#' @keywords otol plant tree chronogram bold
#' @details
#' Generated with make_bold_otol_tree(input = "((Zea mays,Oryza sativa),((Arabidopsis thaliana,(Glycine max,Medicago sativa)),Solanum lycopersicum)Pentapetalae);")
#' usethis::use_data(plant_bold_otol_tree)
"plant_bold_otol_tree"


#' Information on contributors, authors, study ids and clades from studies with chronograms in Open tree of Life
#'
#' @name treebase_cache
#' @docType data
#' @format A list of five data sets
#' \describe{
#'   \item{tb.author.pretty}{A dataframe with two elements: author names and number of studies in TreeBase authored by each}
#'   \item{tb.author.results}{A dataframe with two elements: author names and study identifiers}
#' }
#' @source TreeBASE database, no longer available online \url{https://en.wikipedia.org/wiki/TreeBASE}
#' @keywords treebase tree chronogram author study count id
#' @details
#'
#' Generated with make_treebase_cache()
"treebase_cache"

# #' Authors from studies with chronograms in TreeBASE
# #'
# #'
# #' @name tb.author.pretty
# #' @docType data
# #' @format A data frame with two elements
# #' \describe{
# #'   \item{person}{Character vector of authors of the included studies}
# #'   \item{study.count}{Numeric vector of studies authored by each person}
# #' }
# #' @source TreeBASE database, no longer available online \url{https://en.wikipedia.org/wiki/TreeBASE}
# #' @keywords treebase tree chronogram author study count
# #' @details
# #' Generated with make_treebase_cache()
# "tb.author.pretty"


# #' Information on authors and study ids from chronograms in TreeBASE
# #'
# #'
# #' @name tb.author.results
# #' @docType data
# #' @format A dataframe with two variables
# #' \describe{
# #'   \item{person}{Named character vector of authors of the included studies}
# #'   \item{study}{Named character vector of study identifiers}
# #' }
# #' @source TreeBASE database, no longer available online \url{https://en.wikipedia.org/wiki/TreeBASE}
# #' @keywords treebase tree chronogram study studies id
# #' @details
# #' Generated with make_treebase_cache()
# "tb.author.results"

#' datelifeSummary of a datelifeResult object of all Felidae species.
#'
#' @name felid_gdr_phylo_all
#' @docType data
#' @format A list of three elements, containing the summary of a datelifeResult object
#' \describe{
#'   \item{phylo_all}{List of subset chronograms in phylo format}
#'   \item{taxon_distribution}{A data frame with taxon presence across subset chronograms}
#'   \item{absent_taxa}{A dataframe with names of taxon not found in any chronogram}
#' }
#' @source \url{http://opentreeoflife.org}
#' @keywords otol tree subset chronogram felidae
#' @details
#' Generated with:
#' felid_spp <- make_datelife_query(input = "felidae", get_spp_from_taxon = TRUE)
#' felid_gdr <- get_datelife_result(input = felid_spp, get_spp_from_taxon = TRUE)
#' felid_gdr_phylo_all <- summarize_datelife_result(datelife_result = felid_gdr, taxon_summary = "summary", summary_format = "phylo_all", datelife_query = felid_spp)
#' usethis::use_data(felid_gdr_phylo_all)
"felid_gdr_phylo_all"

#' SDM tree of a datelifeResult object of all Felidae species.
#'
#' @name felid_sdm
#' @docType data
#' @format A list of two elements, containing the summary of a datelifeResult object
#' \describe{
#'   \item{phy}{An ultrametric phylo object with the SDM tree.}
#'   \item{data}{A datelifeResult object with data used to construct phy}
#' }
#' @source \url{http://opentreeoflife.org}
#' @keywords otol tree subset chronogram felidae supertree sdm
#' @details
#' Generated with:
#' felid_spp <- make_datelife_query(input = "felidae", get_spp_from_taxon = TRUE)
#' felid_gdr <- get_datelife_result(input = felid_spp, get_spp_from_taxon = TRUE)
#' felid_sdm <- datelife_result_sdm_phylo(felid_gdr)
#' usethis::use_data(felid_sdm)
"felid_sdm"

#' datelifeResult object of some ants
#'
#' @name some_ants_datelife_result
#' @docType data
#' @format A list of one element, containing a named patristic matrix
#' @source \url{http://opentreeoflife.org}
#' @keywords otol tree subset chronogram ants datelife
#' @details
#' Generated with:
#' some_ants_input <- "(Aulacopone_relicta,(((Myrmecia_gulosa,(Aneuretus_simoni,Dolichoderus_mariae)),((Ectatomma_ruidum,Huberia_brounii),Formica_rufa)),Apomyrma_stygia),Martialis_heureka)Formicidae;"
#' some_ants_datelife_query <- make_datelife_query(input = some_ants_input)
#' some_ants_datelife_result <- get_datelife_result(input = some_ants_datelife_query)
#' usethis::use_data(some_ants_datelife_result)
"some_ants_datelife_result"

#' Long list of >2.7k virus, bacteria, plant and animal taxon names
#'
#' @name subset2_taxa
#' @docType data
#' @format A character vector of length 2778
#' @source \url{https://github.com/phylotastic/rphylotastic/tree/master/tests/testthat}
#' @keywords taxon names subset virus bacteria plant animal
#' @details
#' Generated with:
#' subset2_taxa <- rphylotastic::url_get_scientific_names("https://github.com/phylotastic/rphylotastic/blob/master/tests/testthat/subset2.txt")
#' usethis::use_data(subset2_taxa)
"subset2_taxa"

#' A list with datelifeQuery and datelifeResult objects from a search of taxon names from subset2_taxa
#'
#' @name subset2_search
#' @docType data
#' @format A list with two named elements. datelifeResult object with 24 patristic matrices
#' \describe{
#'   \item{datelife_query}{A datelifeQuery object using names_subset 2 as input.}
#'   \item{datelife_result}{A datelifeResult object resulting from a search of names in datelifeQuery}
#' }
#' @keywords taxon names subset2 datelifeResult
#' @details
#' Generated with:
#' datelife_query <- make_datelife_query(subset2_taxa)
#' datelife_result <- get_datelife_result(datelife_query)
#' subset2_search <- list(query = datelife_query, result = datelife_result)
#' usethis::use_data(subset2_search, overwrite = TRUE)
"subset2_search"

#' A multiPhylo object with trees resulting from a datelife search of some birds and cats species
#'
#' @name birds_and_cats
#' @docType data
#' @format A multiPhylo object
#' @keywords datelifeResult
#' @details
#' Generated with:
#' taxa <- c("Rhea americana", "Pterocnemia pennata", "Struthio camelus", "Gallus", "Felis")
#' birds_and_cats <- datelife_search(input = taxa, summary_format = "phylo_all", get_spp_from_taxon = TRUE)
#' usethis::use_data(birds_and_cats)
"birds_and_cats"
