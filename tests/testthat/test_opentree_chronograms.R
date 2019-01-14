test_that("get_otol_chronograms works", {
	# skip_on_cran()
  #	skip_on_travis()
  # utils::data(opentree_chronograms)
  xx <- get_otol_chronograms(verbose=TRUE, max_tree_count = 10)
  attr(xx)

  # xx <- get_otol_chronograms(verbose=TRUE)
  # sapply(xx$trees, "[", "ott_ids")
  # test that the following makes sense:
  # new.tree <- rotl::get_study_tree(study_id='ot_1000',tree_id='tree1', tip_label="ott_taxon_name")
  # try.tree <- clean_ott_chronogram(new.tree)
  # data.frame(new.tree$tip.label, try.tree$tip.label)
  # is_good_chronogram(try.tree)
  # class(try.tree)
  # inherits(try.tree, "phylo")
  # test that rotl function is generating appropriate trees with ott_taxon_names and ott_ids:
  # t1 <- rotl::get_study_tree(study_id = "ot_1000", tree_id = "tree1", tip_label = "ott_taxon_name")
  # t2 <- rotl::get_study_tree(study_id = "ot_1000", tree_id = "tree1", tip_label = "ott_id")
  # data.frame(t1$tip.label[260:291], t2$tip.label[260:291])
  # match(t1$tip.label, t2$tip.label) # exactly the same as match(t2$tip.label, t1$tip.label)

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
})

test_that("clean_ott_chronogram works as expected", {
  # enhance: there's something wrong when trying to clean the following tree:
  # new.tree2 <- rotl::get_study_tree(study_id = "ot_1041", tree_id = "tree1", tip_label = "ott_id")
  # it is a problem of rotl function, make an issue
  # t1 <- rotl::get_study_tree(study_id = "ot_1041", tree_id = "tree1", tip_label = "original_label")
  # t2 <- rotl::get_study_tree(study_id = "ot_1041", tree_id = "tree1", tip_label = "ott_id", deduplicate = TRUE)
  # t3 <- rotl::get_study_tree(study_id = "ot_1041", tree_id = "tree1", tip_label = "ott_taxon_name")

  # tree1 <- rotl::get_study_tree(study_id = "ot_311", tree_id = "tree1", tip_label = "ott_taxon_name")
  tree1 <- rotl::get_study_tree(study_id = "ot_1250", tree_id = "tree2", tip_label = "ott_taxon_name")
  length(tree1$tip.label) # 31749
  xx <- clean_ott_chronogram(tree1)
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
})
