test_that("results_list_process", {
utils::data(opentree_chronograms)
taxa <- c("Rhea americana", "Pterocnemia pennata", "Struthio camelus")
results_list <- lapply(opentree_chronograms$trees,get_subset_array_dispatch, taxa=taxa, phy=NULL)
expect_gte(length(results_list_process(results_list, taxa, TRUE)), 1)
})

test_that("Summarize as mrca works correctly", {
  utils::data(opentree_chronograms)
  taxa <- c("Rhea americana", "Pterocnemia pennata", "Struthio camelus")
  results_list <- lapply(opentree_chronograms$trees,get_subset_array_dispatch, taxa=taxa, phy=NULL)
  datelife_result <- results_list_process(results_list, taxa, TRUE)
  mrca.vector <- summarize_datelife_result(datelife_result, summary_format="mrca", cache=opentree_chronograms)
  expect_equal(class(mrca.vector), "numeric")
  expect_gte(min(mrca.vector),5)
  expect_lte(max(mrca.vector),150)
})

test_that("Summarize as citations works correctly", {
  utils::data(opentree_chronograms)
  taxa <- c("Rhea americana", "Pterocnemia pennata", "Struthio camelus")
  results_list <- lapply(opentree_chronograms$trees,get_subset_array_dispatch, taxa=taxa, phy=NULL)
  datelife_result <- results_list_process(results_list, taxa, TRUE)
  citation.results <- summarize_datelife_result(datelife_result, summary_format="cit", cache=opentree_chronograms)
  expect_equal(class(citation.results), "character")
  expect_gte(sum(grepl("Prum", citation.results)),1)
 })

test_that("Summarize as newick.all works correctly", {
  utils::data(opentree_chronograms)
  taxa <- c("Rhea americana", "Pterocnemia pennata", "Struthio camelus")
  results_list <- lapply(opentree_chronograms$trees,get_subset_array_dispatch, taxa=taxa, phy=NULL)
  datelife_result <- results_list_process(results_list, taxa, TRUE)
  trees <- summarize_datelife_result(datelife_result, summary_format="newick.all", cache=opentree_chronograms)
  expect_equal(class(trees), "character")
  expect_false(anyNA(trees))
  expect_equal(class(ape::read.tree(text=trees[1])), "phylo")
})

test_that("Summarize as newick.median works correctly", {
  utils::data(opentree_chronograms)
  taxa <- c("Rhea americana", "Pterocnemia pennata", "Struthio camelus")
  results_list <- lapply(opentree_chronograms$trees,get_subset_array_dispatch, taxa=taxa, phy=NULL)
  datelife_result <- results_list_process(results_list, taxa, partial=FALSE)
  tree <- summarize_datelife_result(datelife_result, summary_format="newick.median", cache=opentree_chronograms)
  expect_equal(class(tree), "character")
  expect_false(anyNA(tree))
  expect_equal(class(ape::read.tree(text=tree)), "phylo")
})

test_that("Processing input newick", {
	skip_on_cran()
	skip_on_travis() #b/c no pathd8
  utils::data(opentree_chronograms)
  input.processed <- make_datelife_query(ape::write.tree(ape::rcoal(3, tip.label=c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"))), use_tnrs=FALSE, approximate_match=TRUE)
  expect_equal(class(input.processed$phy),"phylo")
})

# test_that("Processing complex newick works", {
# 	skip_on_cran()
# 	skip_on_travis()
#   utils::data(opentree_chronograms)
# 	expect_error(datelife_search("((((((Typha latifolia,(Phragmites australis,(Sporobolus alterniflorus,Sporobolus pumilus)Sporobolus)PACMAD clade)Poales,(((Hydrilla verticillata,Vallisneria americana)Hydrocharitaceae,Potamogeton perfoliatus),Zostera marina,Ruppia maritima)Alismatales),(Lythrum salicaria,Myriophyllum spicatum)),(Ulva,Caulerpa taxifolia))Chloroplastida,((Skeletonema,(Gomphonema,Didymosphenia geminata)Bacillariophyceae)Bacillariophytina,Prorocentrum)SAR),Microcystis)Eukaryota;", summary_format="phylo.all"), NA)
# })

test_that("datelife_search returns phylo.all", {
  skip_on_cran()
  utils::data(opentree_chronograms)
  phylo.results <- datelife_search(input=c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"), partial=TRUE, use_tnrs=FALSE, approximate_match=TRUE, cache=get("opentree_chronograms"), summary_format="phylo.all")
  expect_equal(class(phylo.results[[1]]),"phylo")
})

test_that("datelife_search returns phylo.sdm", {
  skip_on_cran()
  skip_on_travis() #b/c no pathd8
  utils::data(opentree_chronograms)
phylo.results <- datelife_search(input=c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"), partial=TRUE, use_tnrs=FALSE, approximate_match=TRUE, cache=get("opentree_chronograms"), summary_format="phylo.sdm")
expect_equal(class(phylo.results),"phylo")
expect_true(!is.null(phylo.results$edge.length))
})


test_that("datelife_search returns mrca", {
  skip_on_cran()
  utils::data(opentree_chronograms)
  phylo.results <- datelife_search(input=c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"), partial=TRUE, use_tnrs=FALSE, approximate_match=TRUE, cache=get("opentree_chronograms"), summary_format="mrca")
expect_equal(class(phylo.results),"numeric")
expect_gte(length(phylo.results), 2)
})

test_that("datelife_search returns mrca", {
  skip_on_cran()
  utils::data(opentree_chronograms)
  phylo.results <- datelife_search(input=c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"), partial=TRUE, use_tnrs=FALSE, approximate_match=TRUE, cache=get("opentree_chronograms"), summary_format="mrca")
expect_equal(class(phylo.results),"numeric")
expect_gte(length(phylo.results), 2)
})

test_that("get_datelife_result works", {
  skip_on_cran()
  skip_on_travis() #b/c no pathd8
  utils::data(opentree_chronograms)
  datelife_result.in <- get_datelife_result(input=c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"), partial=TRUE, use_tnrs=FALSE, approximate_match=TRUE, cache=get("opentree_chronograms"))
  # expect_equal(class(datelife_result.in), "list")
  expect_equal(class(datelife_result.in), "datelifeResult")
  expect_gte(length(datelife_result.in), 4) #as of Nov 4, 2016, had length 8
  expect_equal(class(datelife_result.in[[1]]), "matrix")
})

test_that("Processing input string", {
	skip_on_cran()
  utils::data(opentree_chronograms)
  input.processed <- make_datelife_query(c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"), use_tnrs=FALSE, approximate_match=TRUE)
  expect_equal(length(input.processed$cleaned_names),3)
})

test_that("Making OToL and BOLD tree works", {
  skip_on_cran()
  skip_on_travis() #b/c no pathd8
  utils::data(opentree_chronograms)
  phy.pars <- make_bold_otol_tree(input=c("Rhea americana",  "Struthio camelus","Gallus gallus", "Pterocnemia pennata"), marker="COI", otol_version="v2", doML=FALSE)
  phy.ml <- make_bold_otol_tree(input=c("Rhea americana",  "Struthio camelus","Gallus gallus", "Pterocnemia pennata"), marker="COI", otol_version="v2", doML=TRUE)
  expect_equal(class(phy.pars), "phylo")
  expect_equal(class(phy.ml), "phylo")
  expect_gte(max(phy.pars$edge.length), 1)
  expect_lte(min(phy.ml$edge.length), 1)
})

test_that("patristic_matrix_array_congruifyWorks", {
  skip_on_cran()
  skip_on_travis() #b/c no pathd8


#TODO add test here


})

test_that("Congruification works", {
	skip_on_cran()
	skip_on_travis() #b/c no pathd8
  utils::data(opentree_chronograms)
  datelife_result <- get_datelife_result(ape::write.tree(ape::rcoal(3, tip.label=c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"))), partial=TRUE, use_tnrs=FALSE, approximate_match=TRUE)
  # expect_gte(length(datelife_result), 2)
  expect_gte(length(datelife_result), 4)
})

test_that("SDM correctly returns tree", {
  skip_on_cran()
  skip_on_travis() #b/c no pathd8
  utils::data(opentree_chronograms)
  taxa <- c("Rhea americana", "Pterocnemia pennata", "Struthio camelus")
  results_list <- lapply(opentree_chronograms$trees, get_subset_array_dispatch, taxa = taxa, phy = NULL)
  datelife_result <- results_list_process(results_list, taxa, TRUE)
  result.tree <- datelife_result_sdm(datelife_result)$phy
  expect_equal(class(result.tree), "phylo")
})

test_that("use_all_calibrations actually works", {
  skip_on_cran()
  skip_on_travis() #b/c no pathd8
  utils::data(opentree_chronograms)
  results <- use_all_calibrations()
  expect_true(ape::is.ultrametric(results$phy, tol=0.00001, option = 2))
  expect_equal(class(results$phy), "phylo")
})

test_that("Crop plant taxa work", {
  utils::data(opentree_chronograms)
  taxa <- c("Zea mays", "Oryza sativa", "Arabidopsis thaliana", "Glycine max", "Medicago sativa", "Solanum lycopersicum")
  results <- datelife_search(input=taxa, summary_format="phylo.all")
  expect_equal(class(results), "list")
  expect_equal(class(results[[1]]), "phylo")
  expect_gte(length(results), 2)
})

test_that("Processing newick input works", {
  processed <- make_datelife_query("((Zea mays,Oryza sativa),((Arabidopsis thaliana,(Glycine max,Medicago sativa)),Solanum lycopersicum));")
  expect_equal(class(processed$phy), "phylo")
  expect_equal(ape::Ntip(processed$phy), 6)
  expect_equal(ape::Nnode(processed$phy), 5)
  expect_equal(length(processed$cleaned_names), 6)
})

test_that("Crop plant newick works", {
  skip_on_cran()
	skip_on_travis() #b/c no pathd8

trees <- datelife_search(input = "((Zea mays,Oryza sativa),((Arabidopsis thaliana,(Glycine max,Medicago sativa)),Solanum lycopersicum)Pentapetalae);", summary_format = "phylo.all", partial = TRUE, use_tnrs = FALSE, approximate_match = TRUE, cache = opentree_chronograms, method = "PATHd8")
expect_equal(class(trees[[1]]), "phylo")
})

# to test https://github.com/phylotastic/datelife/issues/11
test_that("We don't get negative brlen from pathd8", {
  skip_on_cran()
  skip_on_travis()

  tree <- datelife_search(input = "(((((((Homo sapiens,(Ara ararauna,Alligator mississippiensis)Archosauria)Amniota,Salamandra atra)Tetrapoda,Katsuwonus pelamis)Euteleostomi,Carcharodon carcharias)Gnathostomata,Asymmetron lucayanum)Chordata,(Echinus esculentus,Linckia columbiae)Eleutherozoa)Deuterostomia,(((((Procambarus alleni,Homarus americanus)Astacidea,Callinectes sapidus),(Bombus balteatus,Periplaneta americana)Neoptera)Pancrustacea,Latrodectus mactans)Arthropoda,((Lineus longissimus,(Octopus vulgaris,Helix aspersa)),Lumbricus terrestris))Protostomia);", summary_format = "phylo.median", partial = TRUE, use_tnrs = FALSE, approximate_match = TRUE, method = "PATHd8")
  expect_true(min(tree$edge.length)>=0)

})


# test_that("TNRS with approximate match works", {
	 # taxa <- c("Rhea_americana", "Pterocnemia pennato", "Strutho camelus")
 	# input.processed <- make_datelife_query(taxa, use_tnrs=TRUE, approximate_match=TRUE)
 	# expect_true(all.equal(c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"), input.processed$cleaned_names))
# })


# test_that("TNRS with unmatchable taxa works", {
	# taxa <- c("Rhea_americana", "Pterocnemia pennato", "Oscar the grouch", "Strutho camelus")
 	# input.processed <- make_datelife_query(taxa, use_tnrs=TRUE, approximate_match=TRUE)
 	# expect_true(all.equal(c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"), input.processed$cleaned_names))
# })
