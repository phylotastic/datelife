test_that("get_otol_chronograms works", {
	# skip_on_cran()
  #	skip_on_travis()
  # utils::data(opentree_chronograms)
  xx <- get_otol_chronograms(verbose=TRUE, max_tree_count = 10)
  expect_true(all(c("trees", "authors", "curators", "studies", "dois") %in% names(xx)))
  xx <- get_otol_chronograms(verbose=TRUE)
  table(unlist((sapply(xx$trees, "[", "mapped"))))
  # check the state of trees with ott_id problems:
  rr <- read.csv(file = "data-raw/ott_id_problems_500.csv", row.names = 1)
  tt <- xx$trees[[grep(rr$study.id[1], unlist(xx$studies))]] # get the first tree with ott_ids download problem
  tt$tip.label
  sapply(tt[c("tip.label", "mapped", "ott_ids", "original.tip.label")], length) == length(tt$tip.label)

})

test_that("is_good_chronogram works as expected", {
  utils::data(felid_gdr_phylo_all)
  t1 <- felid_gdr_phylo_all$phylo_all[[1]]
  expect_true(is_good_chronogram(t1))
  # test that all types of not.mapped are detected:
  t1$tip.label[1] <-  "*tip_#1_not_mapped_to_OTT._Original_label_-_Hagensia_havilandi"
  t1$tip.label[2] <-  "*tip_#1_not_mapped_to_OTT._Original_label_-_Felis_silvestris"
  expect_false(is_good_chronogram(t1))
  t1$tip.label <- gsub("_", " ", t1$tip.label)
  expect_false(is_good_chronogram(t1))
  # test presence of unmapped tip labels
  utils::data(problems)
  expect_false(is_good_chronogram(problems[[5]]))
  # test tips with less that two characters as labels
  xx <- problems[[3]]
  xx$tip.label <- sub(" ", "", sub(".*-.", "", xx$tip.label))
  expect_false(is_good_chronogram(xx))
  # enhance: test that there are no duplicated labels in chronogram:
  t1$tip.label[4] <-  "*tip_#1_not_mapped_to_OTT._Original_label_-_Hagensia_havilandi"
  # is_good_chronogram(t1)
  # enhance: test that ott_ids element is ok
  rr <- read.csv(file = "data-raw/ott_id_problems_500.csv", row.names = 1)
  tt <- xx$trees[[grep(rr$study.id[1], unlist(xx$studies))]] # get the first tree with ott_ids download problem
  # is_good_chronogram(tt)

})

test_that("clean_ott_chronogram works as expected", {
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

  # enhance: there's something wrong when trying to clean the following tree:
  # new.tree2 <- rotl::get_study_tree(study_id = "ot_1041", tree_id = "tree1", tip_label = "ott_id")
  # it is a problem of rotl function, make an issue!
  t1 <- rotl::get_study_tree(study_id = "ot_1041", tree_id = "tree1", tip_label = "original_label")
  t2 <- rotl::get_study_tree(study_id = "ot_1041", tree_id = "tree1", tip_label = "ott_id", deduplicate = TRUE)
  t3 <- rotl::get_study_tree(study_id = "ot_1041", tree_id = "tree1", tip_label = "ott_taxon_name")
  rotl::taxonomy_taxon_info(1662)
  grep("Rhinechis", t1$tip.label)
  # the following line takes too long for some reason:
  # tree1 <- rotl::get_study_tree(study_id = "ot_311", tree_id = "tree1", tip_label = "ott_taxon_name")
})

test_that("opentree_chronograms object is ok", {
  utils::data(opentree_chronograms)
  # write(paste(names(opentree_chronograms), collapse = '", "'), file = "data-raw/names.txt")
  # test that all expected elements are in opentree_chronograms:
  expect_true(all(c("trees", "authors", "curators", "studies", "dois") %in% names(opentree_chronograms)))
  # test that all opentree_chronograms elements have the same length:
  # length(opentree_chronograms[[1]])
  expect_true(all(sapply(opentree_chronograms, length) == length(opentree_chronograms$trees)))
  opentree_chronograms$studies[1]
  # check mapping of labels
  names(xx$trees[[1]])
  class(xx$trees[[1]])
  yy <- sapply(xx$trees, function(x) x$tip.label)
  length(yy)
  yy[[1]]
  expect_false(any(grepl("not.mapped", unlist(yy))))
  yy <- sapply(xx$trees, function(x) x$mapped)
  xx$trees[[1]]$tip.label[which(!grepl("ott", yy[[1]]))]
})
