# helper data sets for testing

threebirds <- c("Rhea americana", "Pterocnemia pennata", "Struthio camelus")
threebirds_query <- datelife::make_datelife_query(threebirds)
threebirds_result <- datelife::get_datelife_result(threebirds_query)
threebirds_median <- datelife::summarize_datelife_result(threebirds_result, threebirds_query,
    summary_format = "phylo_median")
threebirds_nbl <- threebirds_median
threebirds_nbl$edge.length <- NULL
threebirds_all <- datelife::summarize_datelife_result(threebirds_result, threebirds_query,
        summary_format = "phylo_all")

utils::data(subset2_taxa)
subset2_query <- datelife::make_datelife_query(subset2_taxa)
subset2_result <- datelife::get_datelife_result(subset2_query)
# utils::data(subset2_search)
# subset2_result <- subset2_search$result
subset2_bestgrove <- datelife::get_best_grove(subset2_result, n=2)$best_grove
subset2_sdm_matrix <- datelife::get_sdm_matrix(subset2_bestgrove)

birds_yellowstone_otoltree <- datelife::get_otol_synthetic_tree(birds_yellowstone)
# names(birds_yellowstone_otoltree)
# data.frame(birds_yellowstone_otoltree$tip.label, birds_yellowstone_otoltree$ott_ids)
birds_yellowstone_query <- datelife::make_datelife_query(input = birds_yellowstone_otoltree)
birds_yellowstone_result <- datelife::get_datelife_result(birds_yellowstone_query)
birds_yellowstone_phyloall <- datelife::summarize_datelife_result(birds_yellowstone_result, birds_yellowstone_query, "phylo_all")
# unname(ape::is.ultrametric(phylo_all))
# unname(sapply(phylo_all, function(x) max(ape::branching.times(x))))
# unname(sapply(datelife_result, function(x) max(x)))
# plot_phylo_all(birds_yellowstone_phyloall)
birds_yellowstone_phylomedian <- datelife::summarize_datelife_result(birds_yellowstone_result,
    birds_yellowstone_query, "phylo_median", summary_print = "citations")
# names(birds_yellowstone_phylomedian)
# plot(birds_yellowstone_phylomedian, cex = 0.25)
# ape::axisPhylo(cex = 0.5)
# birds_yellowstone_phylosdm <- summarize_datelife_result(birds_yellowstone_result,
#     birds_yellowstone_query, "phylo_sdm", summary_print = "citations")
# names(birds_yellowstone_phylosdm)

birdswiki_query <- datelife::make_datelife_query(input = birds_wiki)
birdswiki_result <- datelife::get_datelife_result(birds_wiki)

cetacea_query <- datelife::make_datelife_query(input = "cetacea", get_spp_from_taxon = TRUE)
cetacea_result <- datelife::get_datelife_result(cetacea_query)
cetaceae_phyloall <- datelife::summarize_datelife_result(cetacea_result, cetacea_query, "phylo_all")

pan_query <- datelife::make_datelife_query(input = "Pan", get_spp_from_taxon = TRUE)
pan_result <- datelife::get_datelife_result(pan_query)

catsanddogs_query <- datelife::make_datelife_query(input = c("felis", "canidae"), get_spp_from_taxon = TRUE)
catsanddogs_results <- datelife::get_datelife_result(input = catsanddogs_query)
catsanddogs_phyloall <- datelife::summarize_datelife_result(catsanddogs_results, catsanddogs_query, "phylo_all")
catsanddogs_calibrations <- datelife::extract_calibrations_phylo(catsanddogs_phyloall)
# try(catsanddogs_bold <- make_bold_otol_tree(input = catsanddogs_query, chronogram = FALSE))
# try(ape::is.binary(catsanddogs_bold))
# try(catsanddogs_bladj <- use_calibrations_bladj(phy = catsanddogs_bold, calibrations = catsanddogs_calibrations))
# try(catsanddogs_treepl <- use_calibrations_treePL(phy = catsanddogs_bold, calibrations = catsanddogs_calibrations))
# plot(catsanddogs_treepl$phy)
# catsanddogs_treepl$phy$edge.length



utils::data(felid_sdm)
utils::data(felid_gdr_phylo_all)
utils::data(some_ants_datelife_result)
try(utils::data(plant_bold_otol_tree, subset2_search))

canis <- rotl::tnrs_match_names("canis")  # get the ott id of the genus Canis (dogs)
canis_taxonomy <- rotl::taxonomy_subtree(canis$ott_id)
