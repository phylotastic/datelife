test_that("felis/canidae divergence is accurate", {
    query <- make_datelife_query(input = c("felis", "canidae"), get_spp_from_taxon = TRUE)
    cats_and_dogs_results <- get_datelife_result(input = query)
    matrix_max_ages <- sapply(cats_and_dogs_results, max)
    taxa <- c("felis", "canidae")
    cats_and_dogs <- datelife_search(input = taxa, get_spp_from_taxon = TRUE,
      summary_format = "phylo_all")
    phylo_max_ages <- sapply(cats_and_dogs, function(x) max(ape::branching.times(x)))
    expect_true(all(names(matrix_max_ages) == names(phylo_max_ages)))
    # names(matrix_max_ages) <- names(phylo_max_ages)<- NULL
    ns <- 20
    format(round(sort(matrix_max_ages/2)), nsmall = ns) == format(round(sort(phylo_max_ages)), nsmall = ns)
    # ages from our cache range from 54.9 to 70.9, this includes upper limit confidence interval chronograms
    # timetree study-derived ages range from 39.7 to 67.1. This excludes confidence intervals
    median_phy <- summarize_datelife_result(datelife_result = cats_and_dogs_results, datelife_query = query,
      summary_format = "phylo_median")
    sdm_phy <- summarize_datelife_result(datelife_result = cats_and_dogs_results, datelife_query = query,
        summary_format = "phylo_sdm")
})

# test_that("We get species names back from subspecies and var names",{
#   # rotl::tnrs_match_names("Ceratopogon slossonae")
#   NA
# })

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
