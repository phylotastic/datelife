test_that("get_otol_chronograms works", {
  xx5 <- get_otol_chronograms(verbose=TRUE, max_tree_count = 5)
  expect_true(all(c("trees", "authors", "curators", "studies", "dois") %in% names(xx5)))
  # xx <- get_otol_chronograms(verbose=TRUE)
  # table(unlist((sapply(xx$trees, "[", "mapped"))))
  # check the state of trees with ott_id problems:
  skip_on_cran()
  skip_on_travis()
  skip("check ott_id problems_500")  # read.csv always gives an error with check(), so just run this tests locally
  rr <- read.csv(file = "data-raw/ott_id_problems_500.csv", row.names = 1)
  tt <- opentree_chronograms$trees[[grep(rr$study.id[1], unlist(opentree_chronograms$studies))]] # get the first tree with ott_ids download problem
  # tt$tip.label
  # sapply(tt[c("tip.label", "mapped", "ott_ids", "original.tip.label")], length) == length(tt$tip.label)

})

test_that("is_good_chronogram works as expected", {
  t1 <- felid_gdr_phylo_all$phylo_all[[1]]
  expect_true(is_good_chronogram(t1))
  # test that all types of not.mapped are detected:
  t1$tip.label[1] <-  "*tip_#1_not_mapped_to_OTT._Original_label_-_Hagensia_havilandi"
  t1$tip.label[2] <-  "*tip_#1_not_mapped_to_OTT._Original_label_-_Felis_silvestris"
  expect_warning(xx <- is_good_chronogram(t1))
  expect_false(xx)
  t1$tip.label <- gsub("_", " ", t1$tip.label)
  expect_warning(xx <- is_good_chronogram(t1))
  expect_false(xx)
  # test presence of unmapped tip labels
  utils::data(problems)
  expect_warning(xx <- is_good_chronogram(problems[[5]]))
  expect_false(xx)
  # test tips with less than two characters as labels
  xx <- problems[[3]]
  xx$tip.label <- sub(" ", "", sub(".*-.", "", xx$tip.label))
  expect_warning(xx <- is_good_chronogram(xx))
  expect_false(xx)
  # enhance: test that there are no duplicated labels in chronogram:
  t1$tip.label[4] <-  "*tip_#1_not_mapped_to_OTT._Original_label_-_Hagensia_havilandi"
  # is_good_chronogram(t1)
  # enhance: test that ott_ids element is ok
  skip_on_cran()
  skip_on_travis()
  skip("check ott_id problems_500")  # read.csv always gives an error with check()
  rr <- read.csv(file = "data-raw/ott_id_problems_500.csv", row.names = 1)
  tt <- opentree_chronograms$trees[[grep(rr$study.id[1], unlist(opentree_chronograms$studies))]] # get the first tree with ott_ids download problem
  expect_true(is_good_chronogram(tt))

})

test_that("clean_ott_chronogram works as expected", {
    skip_on_cran()
    skip_on_travis()
    skip("this test takes too long")
  # test class of output:
  tree1 <- rotl::get_study_tree(study_id = "ot_1250", tree_id = "tree2", tip_label = "ott_taxon_name")
  # length(tree1$tip.label) # 31749
  xx <- clean_ott_chronogram(tree1)
  inherits(xx, "phylo")
  # test for duplicated
  tt <- ape::rcoal(10)
  tt$tip.label <- c("*tip_#1_not_mapped_to_OTT._Original_label_-_Elephas_maximus",
                    "Homo sapiens",
                    "Felis silvestris",
                    "*tip_#4_not_mapped_to_OTT._Original_label_-_Elephas_maximus",
                    "Unicorn",
                    "*tip #6 not mapped to OTT. Original label - Homo sapiens",
                    "*tip #7 not mapped to OTT. Original label - Homi sappiens",
                    "*tip #8 not mapped to OTT. Original label - Felix sylvestris",
                    "*tip #9 not mapped to OTT. Original label - Ave",
                    "*tip #10 not mapped to OTT. Original label - Eukarya")
  tt1 <- clean_ott_chronogram(tt)
  correct <- c("Elephas maximus", "Homo sapiens", "Felis silvestris", "Elephas maximus",
  "Unicorn", "Homo sapiens", "Homo sapiens", "Felis silvestris", "Are", "Eukaryota")
  expect_true(all(sapply(1:10, function(x) grepl(correct[x], tt1$tip.label[x]))))  # for one on one comparisons)
  # test that tip.labels after clean_ott_chronogram make sense with a real tree?
  # new.tree <- rotl::get_study_tree(study_id='ot_1000',tree_id='tree1', tip_label="ott_taxon_name")
  # try.tree <- clean_ott_chronogram(new.tree)
  # data.frame(new.tree$tip.label, try.tree$tip.label)
  skip("clean_ott_chronogram works as expected local test")
  # enhance: there's something wrong when trying to clean the following tree:
  # new.tree2 <- rotl::get_study_tree(study_id = "ot_1041", tree_id = "tree1", tip_label = "ott_id")
  # it is a problem of rotl function, make an issue!
  t1 <- rotl::get_study_tree(study_id = "ot_1041", tree_id = "tree1", tip_label = "original_label")
  # following gives an error:
  t2 <- rotl::get_study_tree(study_id = "ot_1041", tree_id = "tree1", tip_label = "ott_id", deduplicate = TRUE)
  # Error in rncl(file = file, ...) :
  # Taxon number 1662 (coded by the token 1662) has already been encountered in this tree. Duplication of taxa in a tree is prohibited.
  t3 <- rotl::get_study_tree(study_id = "ot_1041", tree_id = "tree1", tip_label = "ott_taxon_name")
  rotl::taxonomy_taxon_info(1662)
  grep("Rhinechis", t1$tip.label)
  # the following line takes too long for some reason:
  # tree1 <- rotl::get_study_tree(study_id = "ot_311", tree_id = "tree1", tip_label = "ott_taxon_name")
})

test_that("opentree_chronograms object is ok", {
  # write(paste(names(opentree_chronograms), collapse = '", "'), file = "data-raw/names.txt")
  # test that all expected elements are in opentree_chronograms:
  expect_true(all(c("trees", "authors", "curators", "studies", "dois") %in% names(opentree_chronograms)))
  # test that all opentree_chronograms elements have the same length:
  # length(opentree_chronograms[[1]])
  expect_true(all(sapply(opentree_chronograms, length) == length(opentree_chronograms$trees)))
  # opentree_chronograms$studies[1]
  # check mapping of labels
  skip_on_cran()
  skip_on_travis()
  skip("opentree_chronograms object is ok, local tests")
  # xx object is a tmp opentree_chronograms object
  names(xx$trees[[1]])
  class(xx$trees[[1]])
  yy <- sapply(xx$trees, function(x) x$tip.label)
  length(yy)
  yy[[1]]
  expect_false(any(grepl("not.mapped", unlist(yy))))
  yy <- sapply(xx$trees, function(x) x$mapped)
  xx$trees[[1]]$tip.label[which(!grepl("ott", yy[[1]]))]
})
