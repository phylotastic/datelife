# helper data sets for testing

threebirds <- c("Rhea americana", "Pterocnemia pennata", "Struthio camelus")
threebirds_query <- make_datelife_query(threebirds)
threebirds_result <- get_datelife_result(threebirds_query)
threebirds_median <- summarize_datelife_result(threebirds_result, threebirds_query,
    summary_format = "phylo_median")
threebirds_nbl <- threebirds_median
threebirds_nbl$edge.length <- NULL
threebirds_all <- summarize_datelife_result(threebirds_result, threebirds_query,
        summary_format = "phylo_all")

utils::data(subset2_taxa)
subset2_query <- make_datelife_query(subset2_taxa)
subset2_result <- get_datelife_result(subset2_query)
# utils::data(subset2_search)
# subset2_result <- subset2_search$result
subset2_bestgrove <- get_best_grove(subset2_result, n=2)$best_grove
subset2_sdm_matrix <- get_sdm_matrix(subset2_bestgrove)

birds_yellowstone_otoltree <- get_otol_synthetic_tree(birds_yellowstone)
# names(birds_yellowstone_otoltree)
# data.frame(birds_yellowstone_otoltree$tip.label, birds_yellowstone_otoltree$ott_ids)
birds_yellowstone_query <- make_datelife_query(input = birds_yellowstone_otoltree)
birds_yellowstone_result <- get_datelife_result(birds_yellowstone_query)
birds_yellowstone_phyloall <- summarize_datelife_result(birds_yellowstone_result, birds_yellowstone_query, "phylo_all")
# unname(ape::is.ultrametric(phylo_all))
# unname(sapply(phylo_all, function(x) max(ape::branching.times(x))))
# unname(sapply(datelife_result, function(x) max(x)))
# plot_phylo_all(birds_yellowstone_phyloall)
birds_yellowstone_phylomedian <- summarize_datelife_result(birds_yellowstone_result,
    birds_yellowstone_query, "phylo_median", summary_print = "citations")
# names(birds_yellowstone_phylomedian)
# plot(birds_yellowstone_phylomedian, cex = 0.25)
# ape::axisPhylo(cex = 0.5)
# birds_yellowstone_phylosdm <- summarize_datelife_result(birds_yellowstone_result,
#     birds_yellowstone_query, "phylo_sdm", summary_print = "citations")
# names(birds_yellowstone_phylosdm)

birdswiki_query <- make_datelife_query(input = birds_wiki)
birdswiki_result <- get_datelife_result(birds_wiki)

cetacea_query <- make_datelife_query(input = "cetacea", get_spp_from_taxon = TRUE)
cetacea_result <- get_datelife_result(cetacea_query)
cetaceae_phyloall <- summarize_datelife_result(cetacea_result, cetacea_query, "phylo_all")

pan_query <- make_datelife_query(input = "Pan", get_spp_from_taxon = TRUE)
pan_result <- get_datelife_result(pan_query)

catsanddogs_query <- make_datelife_query(input = c("felis", "canidae"), get_spp_from_taxon = TRUE)
catsanddogs_results <- get_datelife_result(input = catsanddogs_query)
catsanddogs_phyloall <- summarize_datelife_result(catsanddogs_results, catsanddogs_query, "phylo_all")
catsanddogs_calibrations <- get_all_calibrations(input = catsanddogs_phyloall)
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
