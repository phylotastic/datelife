#testing functions that process or check datelife inputs
# test_that("We get species names back from subspecies and var names",{
#   # rotl::tnrs_match_names("Ceratopogon slossonae")
#   NA
# })

test_that("make_datelife_query works from comma separated names", {
  input_processed <- make_datelife_query(input = c("Rhea americana, Pterocnemia pennata, Struthio camelus"),
                                         use_tnrs=FALSE,
                                         approximate_match=TRUE)
  expect_equal(length(input_processed$cleaned_names),
               3)
})

test_that("is_datelife query works", {
  expect_false(is_datelife_query(NA))
  expect_false(is_datelife_query(NULL))
  expect_true(is_datelife_query(threebirds_query))
  # remove class datelifeQuery, but format/structure of the object is ok, it should still work:
  class(threebirds_query) <- "random"
  expect_true(is_datelife_query(threebirds_query))

})

test_that("Mus higher-taxon search is giving species back", {
  skip_on_cran()
  skip_on_travis()
  make_datelife_query("Echinus", get_spp_from_taxon = TRUE)
  make_datelife_query("Mus", get_spp_from_taxon = TRUE)
  # expect_true(length(rphylotastic::taxon_get_species("Mus")) > 0)
})

test_that("make_datelife_query works from phylo object as input", {
    input.processed <- make_datelife_query(ape::write.tree(ape::rcoal(3, tip.label=c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"))), use_tnrs=FALSE, approximate_match=TRUE)
    expect_equal(class(input.processed$phy),"phylo")
})

# test_that("Processing complex newick works", {
# 	skip_on_cran()
# 	skip_on_travis()
#   utils::data(opentree_chronograms)
# 	expect_error(datelife_search("((((((Typha latifolia,(Phragmites australis,(Sporobolus alterniflorus,Sporobolus pumilus)Sporobolus)PACMAD clade)Poales,(((Hydrilla verticillata,Vallisneria americana)Hydrocharitaceae,Potamogeton perfoliatus),Zostera marina,Ruppia maritima)Alismatales),(Lythrum salicaria,Myriophyllum spicatum)),(Ulva,Caulerpa taxifolia))Chloroplastida,((Skeletonema,(Gomphonema,Didymosphenia geminata)Bacillariophyceae)Bacillariophytina,Prorocentrum)SAR),Microcystis)Eukaryota;", summary_format="phylo_all"), NA)
# })

test_that("make_datelife_query works from vector of taxa", {
  input.processed <- make_datelife_query(c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"), use_tnrs=FALSE, approximate_match=TRUE)
  expect_equal(length(input.processed$cleaned_names),3)
})

test_that("make_datelife_query works from newick input", {
  processed <- make_datelife_query("((Zea mays,Oryza sativa),((Arabidopsis thaliana,(Glycine max,Medicago sativa)),Solanum lycopersicum));")
  expect_equal(class(processed$phy), "phylo")
  expect_equal(ape::Ntip(processed$phy), 6)
  expect_equal(ape::Nnode(processed$phy), 5)
  expect_equal(length(processed$cleaned_names), 6)
})

test_that("input_process works", {
	new <- "(((((Pterois miles,Pterois volitans)Pteroinae)Teleostei)Chordata,Lymnaea))Metazoa;"
	phy <- ape::read.tree(text="((Zea mays,Oryza sativa),((Arabidopsis thaliana,(Glycine max,Medicago sativa)),Solanum lycopersicum)Pentapetalae);")
	notnew <- "a,b;"
	expect_s3_class(input_process(c(new, new)), "phylo") #trying to process two phylogenies will give a message
	input <- c(new, new)
	expect_s3_class(input_process(c(phy, phy)), "phylo") #trying to process two phylogenies will give a message
	expect_output(x <- input_process(new), NA)
	expect_output(x <- input_process(notnew), NA)
	expect_output(x <- input_process("purrr"), NA)
	expect_s3_class(x <- input_process(new), "phylo") # output is phylo
	expect_s3_class(x <- input_process(phy), "phylo") # output is phylo
})

test_that("datelife_query_check works with phylo as input", {
    datelife_query_check(datelife_query = felid_sdm$phy)
    datelife_query_check(datelife_query = threebirds_median)
})
