test_that("tnrs_match.phylo works", {
  # phy <- tt
    phy <- ape::rcoal(10, tip.label = c("*tip_#1_not_mapped_to_OTT._Original_label_-_Elephas_maximus",
                                       "Homo sapiens",
                                       "Felis silvestris",
                                       "*tip_#4_not_mapped_to_OTT._Original_label_-_Elephas_maximus",
                                       "Unicorn",
                                       "*tip #6 not mapped to OTT. Original label - Homo sapiens",
                                       "*tip #7 not mapped to OTT. Original label - Homi sappiens",
                                       "*tip #8 not mapped to OTT. Original label - Felix sylvestris",
                                       "*tip #9 not mapped to OTT. Original label - Ave",
                                       "*tip #10 not mapped to OTT. Original label - Eukarya"))
    # start of clean_ott_chronogram functions
    phy.ori <- phy
  	phy <- phylo_tiplabel_underscore_to_space(phy)
   	unmapped.taxa <- unique(c(which(nchar(phy$tip.label) <= 2), which(grepl("not.mapped", phy$tip.label)))) # numeric of indices
   	phy$tip.label[unmapped.taxa] <- sub(".*-.", "", phy$tip.label[unmapped.taxa])  # this gets the original label and gets rid of the not.mapped tag
   	phy$tip.label[unmapped.taxa] <- gsub("aff ", "", phy$tip.label[unmapped.taxa])  # removes aff tag
       phy$tip.label[unmapped.taxa][stringr::str_count(phy$tip.label[unmapped.taxa], " ")>=2] <- gsub("^([^ ]* [^ ]*) .*$", "\\1", phy$tip.label[unmapped.taxa][stringr::str_count(phy$tip.label[unmapped.taxa], " ")>=2])
    tipstodrop <- c()
  	cond <- match(unique(phy$tip.label[unmapped.taxa]), phy$tip.label[-unmapped.taxa])
   	if(any(!is.na(cond))){
   		mm <- match(phy$tip.label[unmapped.taxa], unique(phy$tip.label[unmapped.taxa])[!is.na(cond)])
   		tipstodrop <- c(tipstodrop, unmapped.taxa[!is.na(mm)])
   		unmapped.taxa <- unmapped.taxa[is.na(mm)] # update unmapped.taxa object (removing taxa in mapped tips that are duplicated in unmapped.taxa, preventing unnecesary calls for tnrs_match_names)
   	}
   	dd <- duplicated(phy$tip.label[unmapped.taxa])
   	tipstodrop <- c(tipstodrop, unmapped.taxa[dd])
   	unmapped.taxa <- unmapped.taxa[!dd] # update unmapped.taxa object (by removing duplicated labels within unmapped.taxa)
    # end of clean_ott_chronogram function
    phy2 <- tnrs_match.phylo(phy, tip = unmapped.taxa)
    expect_true(all(c("edge", "edge.length", "tip.label", "Nnode", "mapped", "original.tip.label", "ott_ids") %in% names(phy2)))
    expect_true(all(phy2$original.tip.label == phy$tip.label))
    expect_true(all(grep("original", phy2$mapped) == which(is.na(phy2$ott_ids))))
    # test that ott_ids element generation is ok:
    rr <- read.csv(file = "data-raw/ott_id_problems_500.csv", row.names = 1)
    tt <- xx$trees[[grep(rr$study.id[1], unlist(xx$studies))]] # get the first tree with ott_ids download problem
    length(tt$ott_ids)
    is.null(tt$ott_ids)
})
