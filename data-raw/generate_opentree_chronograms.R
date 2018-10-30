# modified get_otol_chronograms code to retain taxa ott_ids too
opentree_chronograms <- get_otol_chronograms()
length(opentree_chronograms[[1]])
devtools::use_data_raw()
devtools::use_data(opentree_chronograms)
