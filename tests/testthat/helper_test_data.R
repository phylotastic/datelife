# helper data sets for testing

threebirds <- c("Rhea americana", "Pterocnemia pennata", "Struthio camelus")
threebirds_query <- make_datelife_query(threebirds)
threebirds_result <- get_datelife_result(threebirds_query)
threebirds_median <- summarize_datelife_result(threebirds_result, threebirds_query,
    summary_format = "phylo_median")
threebirds_all <- summarize_datelife_result(threebirds_result, threebirds_query,
        summary_format = "phylo_all")

utils::data(subset2_taxa)
subset2_query <- make_datelife_query(subset2_taxa)
subset2_result <- get_datelife_result(subset2_query)
# utils::data(subset2_search)
# subset2_result <- subset2_search$result
subset2_bestgrove <- get_best_grove(subset2_result)$best_grove
# subset2_bestgrove_phyloall <- summarize_datelife_result(subset2_bestgrove, summary_format = "phylo_all")
subset2_sdm_matrix <- get_sdm_matrix(subset2_bestgrove)
xx_bionjs <- ape::bionjs(subset2_sdm_matrix)
xx_triang <- ape::triangMtds(subset2_sdm_matrix)
plot(ape::ladderize(xx_triang), cex = 0.8)
xx_mvrs <- ape::mvrs(subset2_sdm_matrix)
plot(ape::ladderize(xx_bionjs))
subset2_sdm_matrix_filled1 <- ape::additive(subset2_sdm_matrix)
dimnames(subset2_sdm_matrix_filled1) <- list(rownames(subset2_sdm_matrix), colnames(subset2_sdm_matrix))
subset2_sdm_matrix_filled2 <- ape::ultrametric(subset2_sdm_matrix)
dimnames(subset2_sdm_matrix_filled2) <- list(rownames(subset2_sdm_matrix), colnames(subset2_sdm_matrix))
# source("tests/testthat/helper_birdnames.R")
xx <- cluster_patristicmatrix(subset2_sdm_matrix_filled1)
xx <- cluster_patristicmatrix(subset2_sdm_matrix_filled2)
plot(ape::ladderize(xx$upgma), cex =0.8)
plot(ape::ladderize(xx$nj), cex =0.8)
ape::axisPhylo()
ape::is.ultrametric(xx$nj)
xx3 <- ape::bionj(subset2_sdm_matrix_filled1)
# Error in ape::bionj(subset2_sdm_matrix_filled2) :
#   at least one distance was greater than 100

otol <- get_otol_synthetic_tree(xx$upgma$tip.label)
plot(ape::ladderize(otol), cex = 0.8)
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

pan_query <- make_datelife_query(input = "Pan", get_spp_from_taxon = TRUE)
pan_result <- get_datelife_result(pan_query)

catsanddogs_query <- make_datelife_query(input = c("felis", "canidae"), get_spp_from_taxon = TRUE)
catsanddogs_results <- get_datelife_result(input = catsanddogs_query)


utils::data(felid_sdm)
utils::data(felid_gdr_phylo_all)
utils::data(some_ants_datelife_result)
utils::data(plant_bold_otol_tree, subset2_search)
