#' Gets Open Tree of Life synthetic tree of a set of lineages.
#' @inheritParams check_ott_input
#' @param otol_version Version of OToL to use
#' @inheritDotParams make_datelife_query -input
#' @return A phylo object
#' @export
get_otol_synthetic_tree <- function(input = NULL, ott_ids = NULL, otol_version = "v2", ...){
	# input <- birds_yellowstone
	# input <- birds_wiki
	# input <- c(birds_wiki, "ttttt")
	input <- unique(input)
	input_ott_match <- suppressMessages(check_ott_input(input, ott_ids, ...))
	if(length(input_ott_match) < 2){
		message("At least two valid names or numeric ott_ids are needed to get a tree")
		return(NA)
	}
	# enhance: we might need a check of ott_id elements, are they all numeric, are there no NAs, etc.
	# also, another check here, are all ott_ids from valid taxa? this is checked with get_ott_children, but from other functions we should check This
	# enhance: add a class to get_ott_children outputs so its easier to check here if all ott ids are valid taxa, indicate if it has been cleaned, maybe within an atrribute
	# system.time({sapply(rotl::taxonomy_taxon_info(df$ott_id), "[", "flags")})
	# system.time({tnrs_match(rownames(df))}) # this one is faster
	phy <- tryCatch(suppressWarnings(rotl::tol_induced_subtree(ott_ids = input_ott_match, label_format = "name",  otl_v = otol_version)),
					error = function(e){
						message("Some or all input taxa are absent from OToL synthetic tree, look for invalid taxa and clean them from input")
						# this will happen if there are some extinct taxa, barren or any invalid taxa in taxonomy
						# enhance: to avoid it, clean invalid taxa at the beginning
						NA })
	if(length(phy) == 1){
		return(phy)
	}
	if(!ape::is.binary(phy)){
		message(paste0("OToL synthetic tree of input taxa is not fully resolved (",
		phy$Nnode, " nodes/", length(phy$tip.label), " tips)."))
		message("Resolving at random...")
		phy <- ape::multi2di(phy)
	}
	# example of weird behaviour on tip labeling from otol:
	# tnrs <- rotl::tnrs_match_names(c("Staphylococcus aureus", "Bacillus subtilis", "Neisseria meningitidis"))
	# tol_sub <- rotl::tol_induced_subtree(ott_ids = tnrs$ott_id)
	# curl -X POST https://api.opentreeoflife.org/v3/tree_of_life/induced_subtree -H "content-type:application/json" -d '{"ott_ids":[1090496, 1084928, 611778]}'
	# tol_subpost <- ape::read.tree(text = "((((((((((((((((((((((((Bacillus_subtilis_ott1084928)mrcaott6603ott41501)mrcaott6603ott219996)mrcaott6603ott51434)mrcaott6603ott7907)mrcaott1417ott6603)mrcaott1417ott130524)mrcaott609ott1417)mrcaott609ott695050)mrcaott609ott3174)mrcaott218ott609,(((((((mrcaott4291ott4292)mrcaott3204ott3205)mrcaott3204ott109201)mrcaott3204ott935729)mrcaott3204ott460129)Staphylococcus_ott720488)Staphylococcaceae_ott949800)mrcaott1420ott85485)mrcaott218ott1420)mrcaott218ott9340)mrcaott218ott1491)mrcaott217ott218)mrcaott217ott1456)mrcaott47ott53)mrcaott47ott63791)mrcaott47ott10726)mrcaott47ott184124)mrcaott47ott64165)mrcaott47ott2035)mrcaott47ott19087)mrcaott47ott63363,((((((((((((((((((((((((((((Neisseria_meningitidis_ott611778)mrcaott22944ott239089)mrcaott22944ott237408)mrcaott22944ott239091)mrcaott22944ott469870)mrcaott22944ott233491)mrcaott22944ott279619)mrcaott22944ott279613)mrcaott22944ott279628)mrcaott22944ott279610)mrcaott14134ott22944)mrcaott14134ott67965)mrcaott14134ott1011641)Neisseria_ott611812)mrcaott5074ott139204)Neisseriaceae_ott286853)Neisseriales_ott779197)mrcaott90ott5074)mrcaott90ott103)mrcaott90ott11872)mrcaott90ott191429)mrcaott89ott90)mrcaott89ott3892)mrcaott50ott89)mrcaott50ott21523)mrcaott50ott6117)mrcaott50ott107113)mrcaott50ott1100)mrcaott50ott73)mrcaott47ott50;")

	# now include ott_ids as phy element:
	if(is.null(phy$ott_ids)){
		phy$ott_ids <- input_ott_match[match(phy$tip.label, gsub(" ", "_",
		names(input_ott_match)))] # will give NAs where there are mrcaotts, but fixed after
	}
	mrca_index <- grep("mrcaott", phy$tip.label)
	if(length(mrca_index) > 0){
		mrcaott_names <- unlist(lapply(phy$tip.label[mrca_index], recover_mrcaott))
		phy$tip.label[mrca_index] <- names(mrcaott_names)
		phy$ott_ids[mrca_index] <- mrcaott_names
	}
	phy$ott_ids <- as.numeric(phy$ott_ids)
	return(phy)
}

#' Get an otol induced dated subtree from your set of queried taxa
#' @inheritParams check_ott_input
#' @inheritDotParams make_datelife_query -input
#' @return A phylo object with edge length proportional to time in Myrs. It will return NA if any ott_id is invalid.
#' @export
#' @details otol dated tree from Stephen Smith's otol scaling service.
# # ' @examples
#' # if you want to make an ltt plot of a dated OToL tree you'll need to get rid of singleton nodes with ape::collapse.singles
#' # and also probably do phytools::force.ultrametric
get_dated_otol_induced_subtree <- function(input = NULL, ott_ids = NULL, ...){
	# for debugging:
	# utils::data(threebirds.rda)
	# input <- threebirds_query$cleaned_names
	# ott_ids <- threebirds_query$ott_ids
	input_ott_match <- suppressMessages(check_ott_input(input, ott_ids, ...))
	if(length(input_ott_match) < 2){
		message("At least two valid names or numeric ott_ids are needed to get a tree")
		return(NA)
	}
    pp <- tryCatch(httr::POST("http://141.211.236.35:10999/induced_subtree",
						body = list(ott_ids = input_ott_match),
						encode = "json", httr::timeout(10)), error = function(e) NA)
	if(length(pp) > 1){ # this means it retrieved a tree succesfully
		pp <- httr::content(pp)
		rr <- httr::POST("http://141.211.236.35:10999/rename_tree",
		          body = list(newick = pp$newick),
		          encode = "json", httr::timeout(10))
		rr <- httr::content(rr)
		rr <- ape::read.tree(text = rr$newick)
		rr$ott_ids <- ape::read.tree(text = pp$newick)$tip.label
	  return(rr)
	}	else { # this means it errored and we return NA
		message("There was probably a connection error with OToL's tree-dating service; try again later.")
		return(NA)
	}
}
