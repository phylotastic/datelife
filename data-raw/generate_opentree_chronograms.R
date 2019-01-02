#' Get taxonomy for all species within opentree_chronograms object
#' @inheritParams phylo_check
#' @return A cleaned up phylo object
#' @export
#' @examples
#' data(opentree_chronograms)
#' names(opentree_chronograms)
#' get_them_all(opentree_chronograms$trees)
get_them_all <- function(trees = opentree_chronograms$trees){

	all <- unique(unlist(sapply(trees, "[", "ott_ids")))
	all[which(is.na(as.numeric(all)))]
	# problems <- opentree_chronograms$trees[sapply(sapply(trees, "[", "ott_ids"), function(x) any(grepl("tip", x)))]
	problems <- opentree_chronograms$trees[sapply(sapply(trees, "[", "tip.label"), function(x) any(grepl("not mapped", x)))]
	length(problems)
	names(problems)
	# [1] "Upham, N.S. & B.D. Patterson. 2015. Evolution of the caviomorph rodents: a complete phylogeny and timetree of living genera. Pp. 63–120 In: Biology of caviomorph rodents: diversity and evolution (A. I. Vassallo & D. Antenucci, eds.). SAREM Series A, Buenos Aires. "
	# [2] "Prashant P. Sharma, Caitlin M. Baker, Julia G. Cosgrove, Joanne E. Johnson, Jill T. Oberski, Robert J. Raven, Mark S. Harvey, Sarah L. Boyer, Gonzalo Giribet, 2018, 'A revised dated phylogeny of scorpions: Phylogenomic support for ancient divergence of the temperate Gondwanan family Bothriuridae', Molecular Phylogenetics and Evolution, vol. 122, pp. 37-45"
	# [3] "Yun-Yun Lv, Kai He, Sebastian Klaus, Rafe M. Brown, Jia-Tang Li, 2018, 'A comprehensive phylogeny of the genus  Kurixalus  (Rhacophoridae, Anura) sheds light on the geographical range evolution of frilled swamp treefrogs', Molecular Phylogenetics and Evolution, vol. 121, pp. 224-232"
	# [4] "Andressa Paladini, Daniela M. Takiya, Julie M. Urban, Jason R. Cryan, 2018, 'New World spittlebugs (Hemiptera: Cercopidae: Ischnorhininae): Dated molecular phylogeny, classification, and evolution of aposematic coloration', Molecular Phylogenetics and Evolution, vol. 120, pp. 321-334"
	# [5] "Renske E. Onstein, H. Peter Linder, 2016, 'Beyond climate: convergence in fast evolving sclerophylls in Cape and Australian Rhamnaceae predates the mediterranean climate', Journal of Ecology, pp. n/a-n/a"	sapply(problems,names)
	clean_ott_chronogram(phy = problems[[5]])
	sapply(problems, "[", "ott_ids")
	tips <- unique(unlist(sapply(problems, "[", "tip.label")))
	tips2 <- gsub(".* - ", "", tips)
	grep(" ", tips2)
	tips2[grepl("[.*[:blank:].*]{2,}", tips2)]
	# I'm useless with grep, so Brian suggested stringr
	tips2[stringr::str_count(tips2, " ")>=2]
	tips_tnrs <- tnrs_match(tips2)
	tips_tnrs$unique_name
	sum(is.na(tips_tnrs$unique_name))
	names(tips_tnrs)
	tips_tnrs$search_string[is.na(tips_tnrs$unique)]
	length(tips_tnrs$unique_name)
	datelife_search(input = gsub(".* - ", "", sample(tips, 50)))
	all_got <- get_ott_lineage(ott_id = all)
}
