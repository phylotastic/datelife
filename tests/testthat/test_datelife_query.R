#testing functions that process or check datelife inputs
test_that("datelife_query works", {
	cleaned.input_tnrs <- tnrs_match(names = c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"))
	input <- c("Rhea americana", "Pterocnemia pennata", "Struthio camelus")
	df <- get_ott_children(ott_ids = cleaned.input_tnrs$ott_id, ott_rank = "species")
	# head(rownames(df[[1]])[grepl("species", df[[1]]$rank)])
	cleaned_names <- lapply(df, function (x) rownames(x)[grepl("species", x$rank)])
})
test_that("Processing input newick works", {
	skip_on_cran()
#	skip_on_travis() #b/c no pathd8
  skip_on_os("linux") #b/c no pathd8 on travis linux

  utils::data(opentree_chronograms)
  input.processed <- make_datelife_query(ape::write.tree(ape::rcoal(3, tip.label=c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"))), use_tnrs=FALSE, approximate_match=TRUE)
  expect_equal(class(input.processed$phy),"phylo")
})

# test_that("Processing complex newick works", {
# 	skip_on_cran()
# 	skip_on_travis()
#   utils::data(opentree_chronograms)
# 	expect_error(datelife_search("((((((Typha latifolia,(Phragmites australis,(Sporobolus alterniflorus,Sporobolus pumilus)Sporobolus)PACMAD clade)Poales,(((Hydrilla verticillata,Vallisneria americana)Hydrocharitaceae,Potamogeton perfoliatus),Zostera marina,Ruppia maritima)Alismatales),(Lythrum salicaria,Myriophyllum spicatum)),(Ulva,Caulerpa taxifolia))Chloroplastida,((Skeletonema,(Gomphonema,Didymosphenia geminata)Bacillariophyceae)Bacillariophytina,Prorocentrum)SAR),Microcystis)Eukaryota;", summary_format="phylo_all"), NA)
# })

test_that("Processing input string", {
	skip_on_cran()
  utils::data(opentree_chronograms)
  input.processed <- make_datelife_query(c("Rhea americana", "Pterocnemia pennata", "Struthio camelus"), use_tnrs=FALSE, approximate_match=TRUE)
  expect_equal(length(input.processed$cleaned_names),3)
})

test_that("Processing newick input works", {
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
	expect_s3_class(input_process(c(new, new), verbose=TRUE), "phylo") #trying to process two phylogenies will give a message
	input <- c(new, new)
	expect_s3_class(input_process(c(phy, phy), verbose=TRUE), "phylo") #trying to process two phylogenies will give a message
	expect_message(x <- input_process(new, verbose=TRUE)) # when verbose=TRUE it will give a printed message
	expect_message(x <- input_process(phy, verbose=TRUE)) # idem
	expect_message(x <- input_process(notnew, verbose=TRUE)) # idem
	expect_output(x <- input_process(new, verbose=FALSE), NA) # when verbose=FALSE there is no printed message, but it will work with expect_message too
	expect_output(x <- input_process(notnew, verbose=FALSE), NA) # idem
	expect_output(x <- input_process("purrr", verbose=FALSE), NA) # idem
	expect_s3_class(x <- input_process(new, verbose=FALSE), "phylo") # output is phylo
	expect_s3_class(x <- input_process(phy, verbose=FALSE), "phylo") # output is phylo
})

test_that("datelife_query_check works", {
    utils::data(felid_sdm)
    datelife_query_check(datelife_query  = felid_sdm$phy)
})
