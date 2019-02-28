test_that("tree_fix_brlen works", {
    utils::data(plant_bold_otol_tree, subset2_search)
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
    x4 <- summarize_datelife_result(datelife_result = subset2_search$datelife_result, summary_format = "phylo_sdm")
    skip("we need to debug tree_fix_brlen for this example")
    expect_true(ape::is.ultrametric(x4, option = 2))
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

test_that("tree_add_dates works", {
    skip_on_cran()
    skip_on_os("linux") #b/c no mrbayes on travis linux
    utils::data(felid_sdm)
    y <- tree_add_dates(dated_tree = felid_sdm$phy, missing_taxa = letters[1:5])
    missing_taxa <- felid_sdm$phy
    dated_tree <- ape::drop.tip(felid_sdm$phy, c(1,5,9,10,20))
    missing_taxa$edge.length <- NULL
    constraint_tree <- suppressWarnings(geiger::congruify.phylo(reference = phylo_tiplabel_space_to_underscore(dated_tree), target = phylo_tiplabel_space_to_underscore(missing_taxa), scale = NA))
    names(constraint_tree$calibrations)
})
