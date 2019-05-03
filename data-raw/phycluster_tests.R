utils::data(subset2_taxa)
subset2_query <- make_datelife_query(subset2_taxa)
subset2_result <- get_datelife_result(subset2_query)
# utils::data(subset2_search)
# subset2_result <- subset2_search$result
subset2_bestgrove <- get_best_grove(subset2_result)$best_grove
# subset2_bestgrove_phyloall <- summarize_datelife_result(subset2_bestgrove, summary_format = "phylo_all")
subset2_sdm_matrix <- get_sdm_matrix(subset2_bestgrove)


class(subset2_sdm_matrix) <- c("matrix", "datelifeSummaryMatrix")
subset2_sdm_matrix_filled1 <- ape::additive(subset2_sdm_matrix)
dimnames(subset2_sdm_matrix_filled1) <- list(rownames(subset2_sdm_matrix), colnames(subset2_sdm_matrix))
subset2_sdm_matrix_filled2 <- ape::ultrametric(subset2_sdm_matrix)
dimnames(subset2_sdm_matrix_filled2) <- list(rownames(subset2_sdm_matrix), colnames(subset2_sdm_matrix))
# source("tests/testthat/helper_birdnames.R")
xx <- cluster_patristicmatrix(subset2_sdm_matrix_filled1)
phycluster <- xx
# class(xx) <- "datelifeCluster"
choose_cluster(xx)
vv <- datelife_result_variance_matrix(subset2_result)
yy <- ape::mvrs(subset2_sdm_matrix, vv)
plot(yy, cex = 0.8)
ape::is.ultrametric(yy, 2)
vv[is.na(vv)] <- 0
yy <- ape::mvr(subset2_sdm_matrix_filled1, vv)
xx <- cluster_patristicmatrix(subset2_sdm_matrix_filled2)
lapply(xx, function(x) {if(!inherits(x, "phylo")) return(NA)
    ape::is.ultrametric(x)})
for (i in seq(xx)){
    if(!inherits(xx[[i]], "phylo")) next
    plot(ape::ladderize(xx[[i]]), cex = 0.8)
    ape::axisPhylo()
    mtext(names(xx)[i])
}
plot(ape::ladderize(xx$upgma), cex =0.8)
plot(ape::ladderize(xx$nj), cex =0.8)
ape::axisPhylo()
ape::is.ultrametric(xx$nj)
xx3 <- ape::bionj(subset2_sdm_matrix_filled1)
# Error in ape::bionj(subset2_sdm_matrix_filled2) :
#   at least one distance was greater than 100

otol <- get_otol_synthetic_tree(xx$upgma$tip.label)
plot(ape::ladderize(otol), cex = 0.8)
mtext("OTOL")

ss1 <- summary_matrix_to_phylo(subset2_sdm_matrix)
