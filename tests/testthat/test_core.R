# test_that("update_datelife_cache works", {
#     skip_on_cran()
# 	  skip_on_travis() #b/c super time consuming
#     xx <- update_datelife_cache(save = TRUE, file = "/tmp/opentree_chronograms_tmp.RData", verbose = TRUE)  # this works, opentree_chronograms_tmp is saved in tmp
#     expect_true(all(lapply(opentree_chronograms, length) == length(opentree_chronograms$trees)))  # all elements have the same length
# })  # this test takes around 20 min

test_that("results_list_process works", {
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
  mrca.vector <- summarize_datelife_result(datelife_result = datelife_result, summary_format="mrca", cache=opentree_chronograms)
  expect_equal(class(mrca.vector), "numeric")
  expect_gte(min(mrca.vector),5)
  expect_lte(max(mrca.vector),150)
})

test_that("Summarize as citations works correctly", {
  utils::data(opentree_chronograms)
  taxa <- c("Rhea americana", "Pterocnemia pennata", "Struthio camelus")
  results_list <- lapply(opentree_chronograms$trees,get_subset_array_dispatch, taxa=taxa, phy=NULL)
  datelife_result <- results_list_process(results_list, taxa, TRUE)
  citation.results <- summarize_datelife_result(datelife_result = datelife_result, summary_format="cit", cache=opentree_chronograms)
  expect_equal(class(citation.results), "character")
  expect_gte(sum(grepl("Prum", citation.results)),1)
 })

test_that("Summarize as newick_all works correctly", {
  utils::data(opentree_chronograms)
  taxa <- c("Rhea americana", "Pterocnemia pennata", "Struthio camelus")
  results_list <- lapply(opentree_chronograms$trees,get_subset_array_dispatch, taxa=taxa, phy=NULL)
  datelife_result <- results_list_process(results_list, taxa, TRUE)
  trees <- summarize_datelife_result(datelife_result = datelife_result, summary_format="newick_all", cache=opentree_chronograms)
  expect_equal(class(trees), "character")
  expect_false(anyNA(trees))
  expect_equal(class(ape::read.tree(text=trees[1])), "phylo")
})

test_that("add_taxon_distribution argument from summarize_datelife_result() works", {
  utils::data(opentree_chronograms)
  taxa <- c("Rhea americana", "Pterocnemia pennata", "Struthio camelus")
  results_list <- lapply(opentree_chronograms$trees,get_subset_array_dispatch, taxa=taxa, phy=NULL)
  datelife_result <- results_list_process(results_list, taxa, TRUE)
  trees <- summarize_datelife_result(datelife_result = datelife_result, summary_format="newick_all", cache=opentree_chronograms, add_taxon_distribution = "summary")
  # str(trees)
  # trees$taxon_distribution
  # trees$absent_taxa
  expect_gte(length(trees), 3)
  trees2 <- summarize_datelife_result(datelife_result = datelife_result, summary_format="newick_all", cache=opentree_chronograms, add_taxon_distribution = "matrix")
  expect_gte(length(trees2), 3)
  # str(trees2)
  # trees2$taxon_distribution
  # trees2$absent_taxa
})

test_that("datelife_search returns phylo_all", {
  skip_on_cran()
  utils::data(opentree_chronograms)
  datelife_phylo <- datelife_search(input=c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"), partial=TRUE, use_tnrs=FALSE, approximate_match=TRUE, cache=get("opentree_chronograms"), summary_format="phylo_all")
  expect_true(inherits(datelife_phylo,"multiPhylo"))
})


test_that("datelife_search returns phylo_biggest", {
  skip_on_cran()
  utils::data(opentree_chronograms)
  datelife_phylo <- datelife_search(input=c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"), partial=TRUE, use_tnrs=FALSE, approximate_match=TRUE, cache=get("opentree_chronograms"), summary_format="phylo_biggest")
  expect_true(inherits(datelife_phylo,"phylo"))
})


test_that("datelife_search returns phylo_sdm", {
  skip_on_cran()
  skip_on_os("linux") #b/c no pathd8 on travis linux
  utils::data(opentree_chronograms)
  datelife_phylo <- datelife_search(input=c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"), partial=TRUE, use_tnrs=FALSE, approximate_match=TRUE, cache=get("opentree_chronograms"), summary_format="phylo_sdm")
  expect_true(inherits(datelife_phylo,"phylo"))
  expect_true(!is.null(datelife_phylo$edge.length))
})


test_that("datelife_search returns mrca", {
  skip_on_cran()
  utils::data(opentree_chronograms)
  datelife_phylo <- datelife_search(input=c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"), partial=TRUE, use_tnrs=FALSE, approximate_match=TRUE, cache=get("opentree_chronograms"), summary_format="mrca")
expect_equal(class(datelife_phylo),"numeric")
expect_gte(length(datelife_phylo), 2)
})

test_that("datelife_search returns mrca", {
  skip_on_cran()
  utils::data(opentree_chronograms)
  datelife_phylo <- datelife_search(input=c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"), partial=TRUE, use_tnrs=FALSE, approximate_match=TRUE, cache=get("opentree_chronograms"), summary_format="mrca")
expect_equal(class(datelife_phylo),"numeric")
expect_gte(length(datelife_phylo), 2)
})

test_that("get_datelife_result works", {
  skip_on_cran()
  skip_on_os("linux") #b/c no pathd8 on travis linux
  utils::data(opentree_chronograms)
  datelife_result.in <- get_datelife_result(input = c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"), partial=TRUE, use_tnrs=FALSE, approximate_match=TRUE, cache=get("opentree_chronograms"))
  expect_equal(typeof(datelife_result.in), "list")
  expect_s3_class(datelife_result.in, "datelifeResult")
  # expect_gte(length(datelife_result.in), 4) #as of Nov 4, 2016, had length 8
  expect_equal(class(datelife_result.in[[1]]), "matrix")
})

test_that("Making OToL and BOLD tree works", {
  skip_on_cran()
  skip_on_os("linux") #b/c no pathd8 on travis linux
  utils::data(opentree_chronograms)
  phy.pars <- make_bold_otol_tree(input=c("Rhea americana",  "Struthio camelus","Gallus gallus", "Pterocnemia pennata"), marker="COI", otol_version="v2", doML=FALSE)
  phy.ml <- make_bold_otol_tree(input=c("Rhea americana",  "Struthio camelus","Gallus gallus", "Pterocnemia pennata"), marker="COI", otol_version="v2", doML=TRUE)
  expect_equal(class(phy.pars), "phylo")
  expect_equal(class(phy.ml), "phylo")
  expect_gte(max(phy.pars$edge.length), 1)
  expect_lte(min(phy.ml$edge.length), 1)
})

test_that("patristic_matrix_array_congruify Works", {
  skip_on_cran()
  skip_on_os("linux") #b/c no pathd8 on travis linux


#TODO add test here


})

test_that("Congruification works", {
	skip_on_cran()
  skip_on_os("linux") #b/c no pathd8 on travis linux
  utils::data(opentree_chronograms)
  datelife_result <- get_datelife_result(ape::write.tree(ape::rcoal(3, tip.label=c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"))), partial=TRUE, use_tnrs=FALSE, approximate_match=TRUE)
  # expect_gte(length(datelife_result), 2)
  # length(datelife_result) is equal 6
  #TODO add a different test here, testing length is not efficient
})

test_that("Congruification works with treePL", {
    skip_on_cran()
    skip_on_travis() #b/c no treepl on travis
    utils::data(opentree_chronograms)
    datelife_result <- get_datelife_result(ape::write.tree(ape::rcoal(3, tip.label=c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"))), partial=TRUE, use_tnrs=FALSE, approximate_match=TRUE)
    # expect_gte(length(datelife_result), 2)
    #TODO add a different test here, testing length is not efficient
})

test_that("use_all_calibrations actually works", {
  skip_on_cran()
  skip_on_os("linux") #b/c no pathd8 on travis linux
  utils::data(opentree_chronograms)
  # results <- suppressWarnings(use_all_calibrations())
  # # expect_true(ape::is.ultrametric(results$phy, tol=0.0000))
  # expect_true(ape::is.ultrametric(results$phy, option = 2))
  # expect_s3_class(results$phy, "phylo")
})

test_that("Crop plant taxa work", {
  utils::data(opentree_chronograms)
  taxa <- c("Zea mays", "Oryza sativa", "Arabidopsis thaliana", "Glycine max", "Medicago sativa", "Solanum lycopersicum")
  results <- datelife_search(input=taxa, summary_format="phylo_all")
  expect_equal(class(results), "multiPhylo")
  expect_true(inherits(results[[1]], "phylo"))
  expect_gte(length(results), 2)
})



test_that("Crop plant newick works", {
  skip_on_cran()
  skip_on_os("linux") #b/c no pathd8 on travis linux
  trees <- datelife_search(input = "((Zea mays,Oryza sativa),((Arabidopsis thaliana,(Glycine max,Medicago sativa)),Solanum lycopersicum)Pentapetalae);", summary_format = "phylo_all", partial = TRUE, use_tnrs = FALSE, approximate_match = TRUE, cache = opentree_chronograms, dating_method = "PATHd8")
  # expect_s3_class(trees[[1]], "phylo")
})

# to test https://github.com/phylotastic/datelife/issues/11
test_that("That we don't get negative brlen from pathd8", {
  skip_on_cran()
  skip_on_os("linux") #b/c no pathd8 on travis linux
  # tree <- datelife_search(input = "(((((((Homo sapiens,(Ara ararauna,Alligator mississippiensis)Archosauria)Amniota,Salamandra atra)Tetrapoda,Katsuwonus pelamis)Euteleostomi,Carcharodon carcharias)Gnathostomata,Asymmetron lucayanum)Chordata,(Echinus esculentus,Linckia columbiae)Eleutherozoa)Deuterostomia,(((((Procambarus alleni,Homarus americanus)Astacidea,Callinectes sapidus),(Bombus balteatus,Periplaneta americana)Neoptera)Pancrustacea,Latrodectus mactans)Arthropoda,((Lineus longissimus,(Octopus vulgaris,Helix aspersa)),Lumbricus terrestris))Protostomia);", summary_format = "phylo_median", partial = TRUE, use_tnrs = FALSE, approximate_match = TRUE, dating_method = "PATHd8")
  # expect_true(min(tree$edge.length)>=0)
})

test_that("We can get trees of two taxa back", {
skip_on_cran()
res <- datelife_search(c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"), summary_format="data_frame")
two.taxon <- which(res$Ntax==2)[1]
tree <- ape::read.tree(text=as.character(res$Newick[two.taxon]))
expect_equal(ape::Ntip(tree), 2)
expect_gte(tree$edge.length[1], 10)
})


test_that("tree_fix_brlen works", {
    utils::data(plant_bold_otol_tree, dlsearch_subset2)
    x1 <- tree_fix_brlen(tree = plant_bold_otol_tree, fixing_criterion = "negative", fixing_method = 0)
    expect_true(ape::is.ultrametric(x1, option = 2))
    skip_on_cran()
    skip_on_os("linux") #b/c no mrbayes on travis linux
    install.packages("phylocomr", repos = "https://cloud.r-project.org")
    devtools::install_github("ropensci/phylocomr")
    x2 <- tree_fix_brlen(tree = plant_bold_otol_tree, fixing_criterion = "negative", fixing_method = "bladj")
    expect_true(ape::is.ultrametric(x2, option = 2))
    # mrbayes fix brlen test. It takes a while:
    wwdd <- getwd()
    setwd("~/")
    x3 <- tree_fix_brlen(tree = plant_bold_otol_tree, fixing_criterion = "negative", fixing_method = "mrbayes")
    setwd(wwdd)
    expect_true(ape::is.ultrametric(x3, option = 2))
    # prefer option = 2, using the variance to test ultrametricity, cf. E, Paradis'
    # comments on this post http://blog.phytools.org/2017/03/forceultrametric-method-for-ultrametric.html
    # fixing_criterion = "negative", fixing_method = 0
    x4 <- summarize_datelife_result(datelife_result = dlsearch_subset2$datelife_result, summary_format = "phylo_sdm")
    skip("we need to debug tree_fix_brlen for this example")
    expect_true(ape::is.ultrametric(x4, option = 2))
})

test_that("get_otol_synthetic_tree works", {
  otol_tree <- get_otol_synthetic_tree(input = c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"))
  expect_s3_class(otol_tree, "phylo") # output is phylo
  expect_gte(length(otol_tree), 4)  # it has no branch lengths
  # otol_tree <- get_otol_synthetic_tree(input = c("Struthio camelus"))
  # should it return all names always?
  # it should return the same number of names matched by tnrs_match
  # add such a test:
})

test_that("missing_taxa_check works", {
  utils::data(felid_gdr_phylo_all)
  utils::data(felid_sdm)
  mt1 <- missing_taxa_check(missing_taxa = felid_gdr_phylo_all$absent_taxa, dated_tree = felid_sdm$phy)
  mt2 <- missing_taxa_check(missing_taxa = NA, dated_tree = felid_sdm$phy)  # returns "NA"
  mt3 <- missing_taxa_check(missing_taxa = FALSE, dated_tree = felid_sdm$phy)  # returns "FALSE"
  mt4 <- missing_taxa_check(missing_taxa = NULL, dated_tree = felid_sdm$phy)  # does not return error, bc missing_taxa can be NULL
  expect_s3_class(mt1, "data.frame") # output is data.frame
  expect_equal(mt2, "NA") # output is a character vector
  expect_equal(mt3, "FALSE") # output is a character vector
  expect_true(length(mt4) == 0)
})

# until we figure out how to make mrbayes work, this testhat does not work
# test_that("generate_uncertainty gives a tree with different branch lengths and sample_trees works", {
#   skip_on_cran()
#   skip_on_os("linux")
#   utils::data(felid_sdm)
#   xx <- phylo_generate_uncertainty(felid_sdm$phy, uncertainty_method = "mrbayes", age_distribution = "uniform", size = 10)
#   expect_true("consensus_tree" %in% names(xx))
#   expect_true("trees" %in% names(xx))
#   expect_equal(class(xx$trees), "multiPhylo")
#   expect_true(length(xx$trees) == 10)
#   expect_true(all(length(felid_sdm$phy$edge.length) == sapply(xx$trees, function(x) length(x$edge.length))))
#   xx$trees1 <- sample_trees(trees_file = "felid_sdm_phy_mrbayes_uncertainty_uniform.nexus.t", burnin = 0.25, size = 1)
#   expect_equal(class(xx$trees1), "phylo")
#   expect_true(length(xx$trees1) == 4)  # it contains the 4 minimum elements of a phylo object
#   xx$trees1 <- ape::drop.tip(xx$trees1, "fake_outgroup")
#   expect_true(length(felid_sdm$phy$edge.length) == length(xx$trees1$edge.length))
#   expect_true(any(sort(felid_sdm$phy$edge.length) != sort(xx$trees1$edge.length))) #
# })

test_that("tree_add_outgroup and tree_get_singleton_outgroup work", {
  utils::data(felid_sdm)
  xx <- tree_add_outgroup(felid_sdm$phy)
  expect_true("outgroup" %in% xx$tip.label)
  expect_true("outgroup" %in% tree_get_singleton_outgroup(tree = xx))
  expect_true(ape::is.ultrametric(xx))
  expect_true(length(felid_sdm$phy$tip.label) != length(xx$tip.label))

  xx <- tree_add_outgroup(felid_sdm$phy, outgroup = "fake_outgroup")
  expect_true("fake_outgroup" %in% xx$tip.label)
  expect_true("fake_outgroup" %in% tree_get_singleton_outgroup(tree = xx))
  expect_true(ape::is.ultrametric(xx))
  expect_true(length(felid_sdm$phy$tip.label) != length(xx$tip.label))

  felid_sdm$phy$root.edge <- 10
  xx <- tree_add_outgroup(felid_sdm$phy)
  expect_true("outgroup" %in% xx$tip.label)
  expect_true("outgroup" %in% tree_get_singleton_outgroup(tree = xx))
  expect_true(ape::is.ultrametric(xx))
  expect_true(length(felid_sdm$phy$tip.label) != length(xx$tip.label))
  expect_true(is.null(xx$root.edge))

  felid_sdm$phy$edge.length <- NULL
  felid_sdm$phy$root.edge <- NULL
  xx <- tree_add_outgroup(felid_sdm$phy)
  expect_true("outgroup" %in% xx$tip.label)
  expect_true("outgroup" %in% tree_get_singleton_outgroup(tree = xx))
  # expect_false(ape::is.ultrametric(xx)) # this tree has no branch lengths, but if we leave a root edge, the function will run and give a false result
  expect_true(is.null(xx$edge.length))
  expect_true(length(felid_sdm$phy$tip.label) != length(xx$tip.label))

  new_tree <- "((Zea mays,Oryza sativa),((Arabidopsis thaliana,(Glycine max,Medicago sativa)),Solanum lycopersicum));"
  xx <- tree_add_outgroup(new_tree)
  expect_true("outgroup" %in% xx$tip.label)
  expect_true("outgroup" %in% tree_get_singleton_outgroup(tree = xx))
  expect_true(is.null(xx$edge.length))
  expect_true(length(felid_sdm$phy$tip.label) != length(xx$tip.label))

  new_tree <- "(((t4:19.93154846,t3:19.93154846):1.42784927,t2:21.35939773):4.956403277,(t5:7.102366565,t1:7.102366565):19.21343444):48.85405041;"
  xx <- tree_add_outgroup(new_tree, outgroup = "fake_outgroup")
  expect_true("fake_outgroup" %in% xx$tip.label)
  expect_true("fake_outgroup" %in% tree_get_singleton_outgroup(tree = xx))
  expect_true(ape::is.ultrametric(xx))
  expect_true(is.null(xx$root.edge))
})

test_that("patristic_matrix_to_phylo runs", {
    utils::data(some_ants_datelife_result)
    xx <- patristic_matrix_to_phylo(patristic_matrix = some_ants_datelife_result[[1]])
    expect_s3_class(xx, "phylo")
    expect_true(ape::is.ultrametric(xx))
    #make sure it works with missing data:
    withNaN <- some_ants_datelife_result[[1]]
    withNaN[9, 8] <- NaN
    withNaN[8, 9] <- NaN
    xx <- patristic_matrix_to_phylo(patristic_matrix = withNaN)
    expect_s3_class(xx, "phylo")
    expect_true(ape::is.ultrametric(xx))# because NaN, it uses njs, and gives a tree
  })

test_that("tree_add_dates works", {
    skip_on_cran()
    skip_on_os("linux") #b/c no mrbayes on travis linux
    utils::data(felid_sdm)
    y <- tree_add_dates(felid_sdm$phy, missing_taxa = letters[1:5])
  })




    # test_that("bold tree from datelife_search is the same as the one from make_bold_otol_tree", {
# 	tax2 <- c("Homo sapiens", "Macaca mulatta", "Melursus ursinus","Canis lupus pallipes", "Panthera pardus", "Panthera tigris", "Herpestes fuscus", "Elephas maximus", "Haliastur indus")
# 	other <- "(((((((Homo sapiens,(Ara ararauna,Alligator mississippiensis)Archosauria)Amniota,Salamandra atra)Tetrapoda,Katsuwonus pelamis)Euteleostomi,Carcharodon carcharias)Gnathostomata,Asymmetron lucayanum)Chordata,(Echinus esculentus,Linckia columbiae)Eleutherozoa)Deuterostomia,(((((Procambarus alleni,Homarus americanus)Astacidea,Callinectes sapidus),(Bombus balteatus,Periplaneta americana)Neoptera)Pancrustacea,Latrodectus mactans)Arthropoda,((Lineus longissimus,(Octopus vulgaris,Helix aspersa)),Lumbricus terrestris))Protostomia);"
# 	b1 <- make_bold_otol_tree(input = other)
# 	# nb1 <- length(b1$tiplabel)
# 	ed1 <- datelife_search(input = other, summary_format = "phylo_all", partial = TRUE, use_tnrs = FALSE, approximate_match = TRUE, method = "PATHd8", bold = TRUE)
# 	# ned1 <- length(ed1[[length(ed1)]]$tip.label)
# 	# expect_equal(nb1, ned1)
# 	expect_equal(b1$tiplabel, ed1[[length(ed1)]]$tip.label) # tests both trees have the same taxa, in the same number and order. It's ok, cause it should be the same tree
# 	# expect_identical(b1$tiplabel, ed1[[length(ed1)]]$tip.label)
# 	b2 <- make_bold_otol_tree(input = tax2)
# 	# nb1 <- length(b1$tiplabel)
# 	ed2 <- datelife_search(input = tax2, summary_format = "phylo_all", partial = TRUE, use_tnrs = FALSE, approximate_match = TRUE, method = "PATHd8", bold = TRUE)
# 	# ned1 <- length(ed1[[length(ed1)]]$tip.label)
# 	# expect_equal(nb1, ned1)
# 	expect_equal(b2$tip.label, ed2[[length(ed2)]]$tip.label)
# })

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
