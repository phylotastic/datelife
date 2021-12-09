ls(opentree_chronograms)
# [1] "authors"  "curators" "dois"     "studies"  "trees"


tools::showNonASCIIfile("data/opentree_chronograms.rda")

tools::showNonASCII(opentree_chronograms$authors[[2]])

stringi::stri_trans_general(opentree_chronograms$authors[[2]], "latin-ascii")

stringi::stri_enc_toascii(opentree_chronograms$authors[[2]])

stringi::stri_enc_toutf8(opentree_chronograms$authors[[2]])


tools::showNonASCII(unlist(opentree_chronograms$curators))

tools::showNonASCII(unlist(opentree_chronograms$dois))

tools::showNonASCII(unlist(opentree_chronograms$studies))

tools::showNonASCII(unlist(opentree_chronograms$trees))
