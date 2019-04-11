# test_that("update_datelife_cache works", {
#     skip_on_cran()
# 	  skip_on_travis() #b/c super time consuming
#     xx <- update_datelife_cache(save = TRUE, file = "/tmp/opentree_chronograms_tmp.RData", verbose = TRUE)  # this works, opentree_chronograms_tmp is saved in tmp
#     expect_true(all(lapply(opentree_chronograms, length) == length(opentree_chronograms$trees)))  # all elements have the same length
# })  # this test takes around 20 min

test_that("results_list_process works", {
taxa <- c("Rhea americana", "Pterocnemia pennata", "Struthio camelus")
results_list <- lapply(opentree_chronograms$trees,get_subset_array_dispatch, taxa=taxa, phy=NULL)
expect_gte(length(results_list_process(results_list, taxa, TRUE)), 1)
})

test_that("datelife_search returns phylo_all: three birds", {
  datelife_phylo <- datelife_search(input =
    c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"),
    summary_format="phylo_all")
  expect_true(inherits(datelife_phylo,"multiPhylo"))
})
test_that("datelife_search returns phylo_all: Crop plants, vector input", {
  taxa <- c("Zea mays", "Oryza sativa", "Arabidopsis thaliana", "Glycine max", "Medicago sativa", "Solanum lycopersicum")
  results <- datelife_search(input=taxa, summary_format="phylo_all")
  expect_equal(class(results), "multiPhylo")
  expect_true(inherits(results[[1]], "phylo"))
  expect_gte(length(results), 2)
})
test_that("datelife_search returns phylo_all: Crop plants, newick input", {
  trees <- datelife_search(input = "((Zea mays,Oryza sativa),((Arabidopsis thaliana,(Glycine max,Medicago sativa)),Solanum lycopersicum)Pentapetalae);",
  summary_format = "phylo_all", partial = TRUE, use_tnrs = FALSE, approximate_match = TRUE,
  cache = opentree_chronograms, dating_method = "PATHd8")
  # expect_s3_class(trees[[1]], "phylo")
})
test_that("datelife_search returns one phylo_biggest", {
    datelife_phylo <- datelife_search(input =
      c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"),
      summary_format="phylo_biggest")
  expect_true(inherits(datelife_phylo,"phylo"))
})


test_that("datelife_search returns phylo_sdm", {
  skip_on_cran()
  skip_on_os("linux") #b/c no pathd8 on travis linux
  datelife_phylo <- datelife_search(input =
    c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"),
    summary_format="phylo_sdm")
  expect_true(inherits(datelife_phylo,"phylo"))
  expect_true(!is.null(datelife_phylo$edge.length))
})


test_that("datelife_search returns mrca", {
  datelife_phylo <- datelife_search(input =
    c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"),
    summary_format="mrca")
  expect_equal(class(datelife_phylo),"numeric")
  expect_gte(length(datelife_phylo), 2)
})

test_that("get_datelife_result works", {
  # skip_on_cran()
  # skip_on_os("linux") #b/c no pathd8 on travis linux
  datelife_result.in <- get_datelife_result(input = c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"), partial=TRUE, use_tnrs=FALSE, approximate_match=TRUE, cache=get("opentree_chronograms"))
  expect_equal(typeof(datelife_result.in), "list")
  expect_s3_class(datelife_result.in, "datelifeResult")
  # expect_gte(length(datelife_result.in), 4) #as of Nov 4, 2016, had length 8
  expect_equal(class(datelife_result.in[[1]]), "matrix")
})

test_that("Making OToL and BOLD tree works", {
  skip_on_cran()
  skip_on_os("linux") #b/c no pathd8 on travis linux
  phy.pars <- make_bold_otol_tree(input=c("Rhea americana",  "Struthio camelus","Gallus gallus", "Pterocnemia pennata"), marker="COI", otol_version="v2", doML=FALSE)
  phy.ml <- make_bold_otol_tree(input=c("Rhea americana",  "Struthio camelus","Gallus gallus", "Pterocnemia pennata"), marker="COI", otol_version="v2", doML=TRUE)
  expect_equal(class(phy.pars), "phylo")
  expect_equal(class(phy.ml), "phylo")
  expect_gte(max(phy.pars$edge.length), 1)
  expect_lte(min(phy.ml$edge.length), 1)
})


# test_that("Congruification works with pathd8 and treepl", {
# 	skip_on_cran()
#   skip_on_os("linux") #b/c no pathd8 on travis linux
# it does not work anymmore from get_datelife_result, only from summarize with congruify = TRUE
# change this test acordingly ^
#   datelife_result <- get_datelife_result(ape::write.tree(ape::rcoal(3, tip.label=c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"))), partial=TRUE, use_tnrs=FALSE, approximate_match=TRUE)
#   # expect_gte(length(datelife_result), 2)
#   # length(datelife_result) is equal 6
#   #TODO add a different test here, testing length is not efficient
# })


# to test https://github.com/phylotastic/datelife/issues/11
test_that("That we don't get negative brlen from pathd8", {
  skip_on_cran()
  skip_on_os("linux") #b/c no pathd8 on travis linux
  tree <- datelife_search(input = "(((((((Homo sapiens,(Ara ararauna,Alligator mississippiensis)Archosauria)Amniota,Salamandra atra)Tetrapoda,Katsuwonus pelamis)Euteleostomi,Carcharodon carcharias)Gnathostomata,Asymmetron lucayanum)Chordata,(Echinus esculentus,Linckia columbiae)Eleutherozoa)Deuterostomia,(((((Procambarus alleni,Homarus americanus)Astacidea,Callinectes sapidus),(Bombus balteatus,Periplaneta americana)Neoptera)Pancrustacea,Latrodectus mactans)Arthropoda,((Lineus longissimus,(Octopus vulgaris,Helix aspersa)),Lumbricus terrestris))Protostomia);", summary_format = "phylo_median", partial = TRUE, use_tnrs = FALSE, approximate_match = TRUE, dating_method = "PATHd8")
  expect_true(min(tree$edge.length)>=0)
})

test_that("We can get trees of two taxa back", {
    skip_on_cran()
    res <- datelife_search(c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"), summary_format="data_frame")
    two.taxon <- which(res$Ntax==2)[1]
    tree <- ape::read.tree(text=as.character(res$Newick[two.taxon]))
    expect_equal(ape::Ntip(tree), 2)
    expect_gte(tree$edge.length[1], 10)
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
