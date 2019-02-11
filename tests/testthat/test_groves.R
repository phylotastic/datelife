# test that groves functions are doing what they're supposed to do

test_that("test groves", {
    utils::data(names_subset2)
    subset2_query <- make_datelife_query(names_subset2)
    subset2_dl <- get_datelife_result(subset2_query)
    grove_taxa <- filter_for_grove(subset2_dl, criterion = "taxa")
    grove_tree <- filter_for_grove(subset2_dl, criterion = "tree")
})

# test that best grove obtained for median method in datelife_summary function
# is the best one that can be used for sdm (maybe groves with smaller n work for
# sdm and not for median and visceversa)


test_that("median and sdm give a tree when source trees have different degrees of overlapping names", {
    utils::data(names_subset2)
    spp_query <- make_datelife_query(names_subset2)
    spp_dl_result <- get_datelife_result(spp_query)
    # median.result <- NULL
		# overlap <- 2
		# while(!inherits(median.result, "phylo")){
    #   print(overlap)
		#   best_grove <- datelife::filter_for_grove(spp_dl_result, criterion = "taxa", n = overlap)
		#   median.result <- tryCatch(datelife_result_median(best_grove), error = function(e) NULL)
		# 	# sometimes max(branching.times) is off (too big or too small), so we could
		# 	# standardize by real median of original data (max(mrcas)).
		# 	# median.phylo$edge.length <- median.phylo$edge.length * stats::median(mrcas)/max(ape::branching.times(median.phylo))
		#   overlap <- overlap + 1
		# }
    # patristic.array <- patristic_matrix_list_to_array(best_grove)
  	# median.matrix <- summary_patristic_matrix_array(patristic.array)
    # xx <- patristic_matrix_to_phylo(patristic_matrix = median.matrix)
    # yy <- cluster_patristicmatrix(patristic_matrix = median.matrix)
    # ape::is.ultrametric(yy[[1]])
    # any(yy[[1]]$edge.length < 0)
    # min(yy[[1]]$edge.length)
    # yy2 <- phytools::force.ultrametric(yy[[1]])
    # any(yy2$edge.length < 0)
    # min(yy2$edge.length)
    # plot(yy2)
    # xx$clust
    # ape::plot.phylo(xx, cex = 0.5)
    # ape::axisPhylo()
    # ape::is.ultrametric(xx)
    expect_true(inherits(summarize_datelife_result(datelife_query = spp_query,
      datelife_result = spp_dl_result, summary_format = "phylo_median"), "phylo"))
    expect_true(inherits(summarize_datelife_result(datelife_query = spp_query,
      datelife_result = spp_dl_result, summary_format = "phylo_sdm"), "phylo"))
    expect_equal(class(summarize_datelife_result(datelife_query = spp_query,
      datelife_result = spp_dl_result, summary_format = "newick_median")), "character")
    expect_equal(class(summarize_datelife_result(datelife_query = spp_query,
      datelife_result = spp_dl_result, summary_format = "newick_sdm")), "character")
})
