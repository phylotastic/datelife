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

test_that("Summarize as newick.all works correctly", {
  utils::data(opentree_chronograms)
  taxa <- c("Rhea americana", "Pterocnemia pennata", "Struthio camelus")
  results_list <- lapply(opentree_chronograms$trees,get_subset_array_dispatch, taxa=taxa, phy=NULL)
  datelife_result <- results_list_process(results_list, taxa, TRUE)
  trees <- summarize_datelife_result(datelife_result = datelife_result, summary_format="newick.all", cache=opentree_chronograms)
  expect_equal(class(trees), "character")
  expect_false(anyNA(trees))
  expect_equal(class(ape::read.tree(text=trees[1])), "phylo")
})

test_that("add_taxon_distribution argument from summarize_datelife_result() works", {
  utils::data(opentree_chronograms)
  taxa <- c("Rhea americana", "Pterocnemia pennata", "Struthio camelus")
  results_list <- lapply(opentree_chronograms$trees,get_subset_array_dispatch, taxa=taxa, phy=NULL)
  datelife_result <- results_list_process(results_list, taxa, TRUE)
  trees <- summarize_datelife_result(datelife_result = datelife_result, summary_format="newick.all", cache=opentree_chronograms, add_taxon_distribution = "summary")
  # str(trees)
  # trees$taxon_distribution
  # trees$absent_taxa
  expect_gte(length(trees), 3)
  trees2 <- summarize_datelife_result(datelife_result = datelife_result, summary_format="newick.all", cache=opentree_chronograms, add_taxon_distribution = "matrix")
  expect_gte(length(trees2), 3)
  # str(trees2)
  # trees2$taxon_distribution
  # trees2$absent_taxa
})

test_that("Summarize as newick.median works correctly", {
  utils::data(opentree_chronograms)
  taxa <- c("Rhea americana", "Pterocnemia pennata", "Struthio camelus")
  results_list <- lapply(opentree_chronograms$trees,get_subset_array_dispatch, taxa=taxa, phy=NULL)
  datelife_result <- results_list_process(results_list, taxa, partial=FALSE)
  tree <- summarize_datelife_result(datelife_result = datelife_result, summary_format="newick.median", cache=opentree_chronograms)
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
# 	expect_error(datelife_search("((((((Typha latifolia,(Phragmites australis,(Sporobolus alterniflorus,Sporobolus pumilus)Sporobolus)PACMAD clade)Poales,(((Hydrilla verticillata,Vallisneria americana)Hydrocharitaceae,Potamogeton perfoliatus),Zostera marina,Ruppia maritima)Alismatales),(Lythrum salicaria,Myriophyllum spicatum)),(Ulva,Caulerpa taxifolia))Chloroplastida,((Skeletonema,(Gomphonema,Didymosphenia geminata)Bacillariophyceae)Bacillariophytina,Prorocentrum)SAR),Microcystis)Eukaryota;", summary_format="phylo_all"), NA)
# })

test_that("datelife_search returns phylo_all", {
  skip_on_cran()
  utils::data(opentree_chronograms)
  datelife_phylo <- datelife_search(input=c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"), partial=TRUE, use_tnrs=FALSE, approximate_match=TRUE, cache=get("opentree_chronograms"), summary_format="phylo_all")
  expect_equal(class(datelife_phylo[[1]]),"phylo")
})

test_that("datelife_search returns phylo.sdm", {
  skip_on_cran()
  skip_on_travis() #b/c no pathd8
  utils::data(opentree_chronograms)
datelife_phylo <- datelife_search(input=c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"), partial=TRUE, use_tnrs=FALSE, approximate_match=TRUE, cache=get("opentree_chronograms"), summary_format="phylo.sdm")
expect_equal(class(datelife_phylo),"phylo")
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
  skip_on_travis() #b/c no pathd8
  utils::data(opentree_chronograms)
  datelife_result.in <- get_datelife_result(input = c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"), partial=TRUE, use_tnrs=FALSE, approximate_match=TRUE, cache=get("opentree_chronograms"))
  expect_equal(typeof(datelife_result.in), "list")
  expect_s3_class(datelife_result.in, "datelifeResult")
  # expect_gte(length(datelife_result.in), 4) #as of Nov 4, 2016, had length 8
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
  # length(datelife_result) is equal 6
  #TODO add a different test here, testing length is not efficient
})

test_that("Congruification works with treePL", {
    skip_on_cran()
    skip_on_travis() #b/c no treepl
    utils::data(opentree_chronograms)
    datelife_result <- get_datelife_result(ape::write.tree(ape::rcoal(3, tip.label=c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"))), partial=TRUE, use_tnrs=FALSE, approximate_match=TRUE)
    # expect_gte(length(datelife_result), 2)
    #TODO add a different test here, testing length is not efficient
})

test_that("SDM correctly returns tree", {
  skip_on_cran()
  skip_on_travis() #b/c no pathd8
  utils::data(opentree_chronograms)
  taxa <- c("Rhea americana", "Pterocnemia pennata", "Struthio camelus")
  results_list <- lapply(opentree_chronograms$trees, get_subset_array_dispatch, taxa=taxa, phy=NULL)
  datelife_result <- results_list_process(results_list, taxa, TRUE)
  result.tree <- datelife_result_sdm(datelife_result)$phy
  expect_equal(class(result.tree), "phylo")
})

test_that("use_all_calibrations actually works", {
  skip_on_cran()
  skip_on_travis() #b/c no pathd8
  utils::data(opentree_chronograms)
  results <- use_all_calibrations()
  # expect_true(ape::is.ultrametric(results$phy, tol=0.0000))
  expect_true(ape::is.ultrametric(results$phy, option = 2))
  expect_s3_class(results$phy, "phylo")
})

test_that("Crop plant taxa work", {
  utils::data(opentree_chronograms)
  taxa <- c("Zea mays", "Oryza sativa", "Arabidopsis thaliana", "Glycine max", "Medicago sativa", "Solanum lycopersicum")
  results <- datelife_search(input=taxa, summary_format="phylo_all")
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
  trees <- datelife_search(input = "((Zea mays,Oryza sativa),((Arabidopsis thaliana,(Glycine max,Medicago sativa)),Solanum lycopersicum)Pentapetalae);", summary_format = "phylo_all", partial = TRUE, use_tnrs = FALSE, approximate_match = TRUE, cache = opentree_chronograms, dating_method = "PATHd8")
  expect_s3_class(trees[[1]], "phylo")
})

# to test https://github.com/phylotastic/datelife/issues/11
test_that("We don't get negative brlen from pathd8", {
  skip_on_cran()
  skip_on_travis()
  tree <- datelife_search(input = "(((((((Homo sapiens,(Ara ararauna,Alligator mississippiensis)Archosauria)Amniota,Salamandra atra)Tetrapoda,Katsuwonus pelamis)Euteleostomi,Carcharodon carcharias)Gnathostomata,Asymmetron lucayanum)Chordata,(Echinus esculentus,Linckia columbiae)Eleutherozoa)Deuterostomia,(((((Procambarus alleni,Homarus americanus)Astacidea,Callinectes sapidus),(Bombus balteatus,Periplaneta americana)Neoptera)Pancrustacea,Latrodectus mactans)Arthropoda,((Lineus longissimus,(Octopus vulgaris,Helix aspersa)),Lumbricus terrestris))Protostomia);", summary_format = "phylo.median", partial = TRUE, use_tnrs = FALSE, approximate_match = TRUE, dating_method = "PATHd8")
  expect_true(min(tree$edge.length)>=0)
})

test_that("We can get trees of two taxa back", {
skip_on_cran()
res <- datelife_search(c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"), summary_format="data.frame")
two.taxon <- which(res$Ntax==2)[1]
tree <- ape::read.tree(text=as.character(res$Newick[two.taxon]))
expect_equal(ape::Ntip(tree), 2)
expect_gte(tree$edge.length[1], 10)
})

test_that("input_process works", {
	new <- "(((((Pterois miles,Pterois volitans)Pteroinae)Teleostei)Chordata,Lymnaea))Metazoa;"
	phy <- ape::read.tree(text="((Zea mays,Oryza sativa),((Arabidopsis thaliana,(Glycine max,Medicago sativa)),Solanum lycopersicum)Pentapetalae);")
	notnew <- "a,b;"
	expect_error(input_process(c(new, new)), verbose=TRUE) #trying to process two phylogenies will give an error
	expect_error(input_process(c(phy, phy)), verbose=TRUE) #trying to process two phylogenies will give an error
	expect_output(x <- input_process(new, verbose=TRUE)) # when verbose=TRUE it will give a printed message
	expect_output(x <- input_process(phy, verbose=TRUE)) # idem
	expect_output(x <- input_process(notnew, verbose=TRUE)) # idem
	expect_output(x <- input_process(new, verbose=FALSE), NA) # when verbose=FALSE there is no printed message
	expect_output(x <- input_process(notnew, verbose=FALSE), NA) # idem
	expect_output(x <- input_process("purrr", verbose=FALSE), NA) # idem
	expect_type(x <- input_process(notnew, verbose=FALSE), "logical") # output is NA
	expect_type(x <- input_process("purrr", verbose=FALSE), "logical") # output is NA
	expect_s3_class(x <- input_process(new, verbose=FALSE), "phylo") # output is phylo
	expect_s3_class(x <- input_process(phy, verbose=FALSE), "phylo") # output is phylo
})

test_that("tree_fix_brlen works", {
    utils::data(plant_bold_otol_tree)
    x1 <- tree_fix_brlen(tree = plant_bold_otol_tree, fixing_criterion = "negative", fixing_method = 0)
    expect_true(ape::is.ultrametric(x1, option = 2))
    skip_on_cran()
    skip_on_travis() #bc no mrbayes
    install.packages("phylocomr", repos = "https://cloud.r-project.org")
    devtools::install_github("ropensci/phylocomr")
    x2 <- tree_fix_brlen(tree = plant_bold_otol_tree, fixing_criterion = "negative", fixing_method = "bladj")
    expect_true(ape::is.ultrametric(x2, option = 2))
    # mrbayes fix brlen test. It takes a while:
    wwdd <- getwd()
    setwd("~/")
    x3 <- tree_fix_brlen(tree = plant_bold_otol_tree, fixing_criterion = "negative", fixing_method = "mrbayes")
    setwd(wwdd)
    expect_true(ape::is.ultrametric(x3, option = 2))  # prefer option = 2, using the variance to test ultrametricity, cf. E, Paradis' comments on this post http://blog.phytools.org/2017/03/forceultrametric-method-for-ultrametric.html
})

test_that("get_otol_synthetic_tree works", {
  otol_tree <- get_otol_synthetic_tree(input = c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"))
  expect_s3_class(otol_tree, "phylo") # output is phylo
  expect_gte(length(otol_tree), 4)  # it has no branch lengths
  # should it return all names always?
  # it should return the same number of names matched by input_tnrs
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

test_that("generate_uncertainty gives a tree with different branch lengths", {
  utils::data(felid_sdm)
  xx <- phy_generate_uncertainty(felid_sdm$phy, n = 10)
  expect_equal(class(xx), "multiPhylo")
  expect_true(all(length(felid_sdm$phy$edge.length) == sapply(xx, function(x) length(x$edge.length))))
  xx <- phy_generate_uncertainty(felid_sdm$phy, n = 1)
  expect_equal(class(xx), "phylo")
  expect_true(length(felid_sdm$phy$edge.length) == length(xx$edge.length))
  expect_true(any(sort(felid_sdm$phy$edge.length) != sort(xx$edge.length))) #
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
