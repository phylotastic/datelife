# helper data sets for testing

threebirds <- c("Rhea americana", "Pterocnemia pennata", "Struthio camelus")
threebirds_query <- make_datelife_query(threebirds)
threebirds_result <- get_datelife_result(threebirds_query)
threebirds_median <- summarize_datelife_result(threebirds_result, threebirds_query,
    summary_format = "phylo_median")
threebirds_all <- summarize_datelife_result(threebirds_result, threebirds_query,
        summary_format = "phylo_all")

names(threebirds_all)
utils::data(subset2_taxa)
subset2_query <- make_datelife_query(subset2_taxa)
subset2_result <- get_datelife_result(subset2_query)
utils::data(subset2_search)
subset2_result <- subset2_search$result
subset2_bestgrove <- get_best_grove(subset2_result)$best_grove
# subset2_bestgrove_phyloall <- summarize_datelife_result(subset2_bestgrove, summary_format = "phylo_all")
subset2_sdm_matrix <- get_sdm_matrix(subset2_bestgrove)

# source("tests/testthat/helper_birdnames.R")

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
