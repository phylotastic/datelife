
test_that("a known datelife_search run works as it should", {
  expect_no_error(
    datelife_query <- make_datelife_query(
      input = c("Delphinus_delphis", 
                "Gallus gallus", 
                "elephas Maximus", 
                "felis_catus", 
                "homo-sapiens"))
  )

  expect_no_error(
    datelife_result <- get_datelife_result_datelifequery(
      datelife_query = datelife_query,
      partial = TRUE,
      cache = "opentree_chronograms")
    )
  
  expect_true(length(datelife_result) > 0)
  
  expect_no_error(
    res <- summarize_datelife_result(
              datelife_result = datelife_result,
              datelife_query = datelife_query,
              summary_format = "phylo_median",
              na_rm = FALSE,
              summary_print = c("citations", "taxa"),
              taxon_summary = c("none", "summary", "matrix"),
              criterion = "taxa")
    )
  
  expect_no_error(
    taxon_summ <- get_taxon_summary(
      datelife_result = datelife_result,
      datelife_query = datelife_query)
  )
  
  expect_no_error(
    datelifeSearch <- datelife_search(
      datelife_query, 
      summary_format = "phylo_median")
  )
}) 

##########################################
test_that("datelife_use with bladj works", {
  du <- datelife_use(
    input = "Rhea americana, Struthio camelus, Gallus gallus",
    each = FALSE,
    dating_method = "bladj"
  )
  expect_true("phylo" %in% class(du))
  expect_error(datelife_use_datelifequery())
  # testing that phylo workflows work
  # ds <- datelife_search(input = du, summary_format = "citations")
  expect_warning(match_all_calibrations(phy = du, calibrations = NULL))
  du00 <- du0 <- du
  du0$edge.length <- NULL
  expect_warning(extract_calibrations_phylo(input = du0))
  expect_error(extract_calibrations_phylo(input = attributes(du)$datelife_query))
  # test a datelifeQuery input with NA as phy:
  dq <- attributes(du)$datelife_query
  dq$phy <- NA
  expect_warning(datelife_use_datelifequery(datelife_query = dq))
  # when input is a datelifeQuery it returns it:
  # make_datelife_query(input = attributes(du)$datelife_query)
  # test that other bladj workflows work:
  sapply(c("mean", "min", "max"), function(x) {
    use_calibrations_bladj(
      phy = du0,
      calibrations = attributes(du)$datelife_calibrations,
      type = x
    )
  })
  # test that absence of phylogeny or no branch lengths returns warning
  expect_warning(use_calibrations_pathd8(phy = NA))
  expect_error(use_calibrations_bladj(phy = NA))
  expect_warning(use_calibrations_pathd8(phy = du0, calibrations = NULL))
  expect_error(use_calibrations_bladj(phy = du0, calibrations = NULL))
  du00$edge.length <- c(NaN, NaN, NaN, NaN)
  expect_warning(use_calibrations_pathd8(phy = du00))
  # test get calibrations from a vector (calls datelife_search and extract_calibrations_phylo)
  gc <- get_calibrations_vector(input = c("Rhea americana", "Struthio camelus", "Gallus gallus"))
  # testing that you can get a datelife result from a tree:
  dr <- get_datelife_result(input = du)
  is_datelife_result_empty(dr)
  extract_calibrations_dateliferesult(input = dr)
  # test an empty datelifeResult object
  input_process(input = du0)
})

##########################################
test_that("datelife_use with pathd8 works", {
  # du <- datelife_use(
  #   input = "Rhea americana, Struthio camelus, Gallus gallus",
  #   each = FALSE,
  #   dating_method = "pathd8"
  # )  
  # Error in `if (nrow(calibs$present_calibrations) < 1) {
  #   warning("\nDating analysis is not possible with this set of calibrations.")
  #   return(NA)
  # }`: argument is of length zero
  # Backtrace:
  #     ▆
  #  1. └─datelife::datelife_use(...) at test_main.R:94:3
  #  2.   └─datelife::datelife_use_datelifequery(...) at datelife/R/datelife_use.R:56:3
  #  3.     └─datelife::use_all_calibrations(...) at datelife/R/datelife_use.R:98:3
  #  4.       └─datelife::use_calibrations(...) at datelife/R/calibrations_use.R:51:5
  #  5.         └─datelife::use_calibrations_pathd8(phy, calibrations, ...) at datelife/R/calibrations_use.R:167:7
})
#############################################################################
test_that("input processing a newick string and multiPhylo workflows work", {
  newick <- "(Gallus_gallus,(Rhea_americana,Struthio_camelus)Palaeognathae)Aves;"
  phylo <- input_process(newick)
  newickBL <- "(Gallus_gallus:165.8333333,(Rhea_americana:128,Struthio_camelus:128)Palaeognathae:37.83333333)Aves;"
  phyloBL <- input_process(newickBL)
  multiphy <- structure(list(phylo, phylo), class = "multiPhylo")
  # there are no ages in multiphy:
  expect_warning(extract_calibrations_phylo(input = multiphy))
  # multiphyBL has one tree with no branch lengths nad on etree with branch lengths:
  multiphyBL <- structure(list(phylo, phyloBL), class = "multiPhylo")
  expect_warning(calibs <- extract_calibrations_phylo(input = multiphyBL))
  # expect error when calibrations are not provided
  expect_error(use_calibrations(phy = multiphyBL, calibrations = NULL))
  # expect error when calibrations are not congruified:
  expect_error(use_calibrations(phy = multiphyBL, calibrations = calibs))
  # congruifying calibrations:
  # congruify_and_mrca_multiPhylo(phy = phylo, source_chronograms = multiphyBL)
  matched_calibs <- datelife:::match_all_calibrations(phy = phylo, calibrations = calibs)
  inherits(matched_calibs, "list")
  # FIX following, giving error:
  # expect_no_error(
  #   use_all_calibrations(phy = matched_calibs$phy, 
  #                        calibrations = matched_calibs$matched_calibrations, 
  #                        each = TRUE)
  # )
  # ... Using calibrations to date a tree topology.
  # ... Using median ages as secondary calibrations with BLADJ.
  # Error in if (nrow(calibrations) < 1) { : argument is of length zero
  # do we need the following test?
  # phylo$tip.label <- c("rooster", "nandu", "ostrich")
  # expect_warning(match_all_calibrations(phy = phylo, calibrations = calibs))
})

#################################
test_that("object checks work", {
  expect_warning(match_all_calibrations(phy = NULL))
  expect_error(use_calibrations(phy = NULL))
  expect_error(get_calibrations_datelifequery(datelife_query = NULL))
})

###############################################
test_that("you can load opentree_chronograms",{
  data(opentree_chronograms, package = "datelife")
  expect_true(all(c("authors","curators","dois","studies","trees") %in% ls(opentree_chronograms)))
})


# test match_all_calibrations when all(all_nodes < ape::Ntip(phy)) is not TRUE

# test make_bold_otol_tree that does not get a phylo object

# test all function check_conflicting_calibrations in use_calibrations_bladj.R
