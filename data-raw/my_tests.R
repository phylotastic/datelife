i <- "Rhea americana, Struthio camelus, Gallus gallus"
dq <- make_datelife_query(i)
mbot <- make_bold_otol_tree(dq)
dq$phy <- mbot
gcdq <- get_calibrations_datelifequery(
  input = dq,
  each = FALSE
)
phyloall <- datelife_search(
  input = dq,
  summary_format = "phylo_all"
)

dr <- get_datelife_result_datelifequery(
  datelife_query = dq,
  partial = TRUE,
  use_tnrs = FALSE,
  approximate_match = FALSE,
  cache = "opentree_chronograms",
  verbose = TRUE
)

res <- summarize_datelife_result(
  datelife_result = dr,
  datelife_query = dq,
  summary_format = "phyloall",
  partial = TRUE,
  cache = "opentree_chronograms",
  summary_print = c("citations", "taxa"),
  taxon_summary = "summary",
  verbose = TRUE,
  criterion = "taxa"
)
