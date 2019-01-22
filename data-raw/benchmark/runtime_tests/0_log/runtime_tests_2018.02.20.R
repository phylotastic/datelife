# Martes 20 de Febrero 2018

# fix building errors: packages not in cran will always give an error if noted as dependencies √√
# we should change them to suggests: √
devtools::use_package("phylocomr", "Suggests")
# I put paleotree and rphylotastic as suggests too √

# fix building errors: non-ascii encoding in author.pretty and author.results objects from contributorcache.RData file
# changed the code but also needed to update contributorcache.RData object √
make_contributor_cache()
[1] "1 of 808 done"
[1] "2 of 808 done"
[1] "3 of 808 done"
[1] "4 of 808 done"
[1] "5 of 808 done"
[1] "6 of 808 done"
[1] "7 of 808 done"
[1] "8 of 808 done"
[1] "9 of 808 done"
[1] "10 of 808 done"
[1] "11 of 808 done"
[1] "12 of 808 done"
[1] "13 of 808 done"
[1] "14 of 808 done"
[1] "15 of 808 done"
[1] "16 of 808 done"
[1] "17 of 808 done"
[1] "18 of 808 done"
[1] "19 of 808 done"
[1] "20 of 808 done"
[1] "21 of 808 done"
[1] "22 of 808 done"
[1] "23 of 808 done"
No encoding supplied: defaulting to UTF-8.
[1] "24 of 808 done"
[1] "25 of 808 done"
No encoding supplied: defaulting to UTF-8.
[1] "26 of 808 done"
[1] "27 of 808 done"
No encoding supplied: defaulting to UTF-8.
[1] "28 of 808 done"
[1] "29 of 808 done"
[1] "30 of 808 done"
[1] "31 of 808 done"
[1] "32 of 808 done"
[1] "33 of 808 done"
[1] "34 of 808 done"
[1] "35 of 808 done"
[1] "36 of 808 done"
[1] "37 of 808 done"
[1] "38 of 808 done"
[1] "39 of 808 done"
[1] "40 of 808 done"
[1] "41 of 808 done"
[1] "42 of 808 done"
[1] "43 of 808 done"
server error for doi “http://dx.doi.org/10.1016/j.ympev.2015.09.014”, you may want to try again.
no results with relavency score greater than ‘min.relevance’ successfully retrieved
Error in (function (bibtype, textVersion = NULL, header = NULL, footer = NULL,  : 
  argument "bibtype" is missing, with no default
[1] "44 of 808 done"
################################ ORIGINAL FUNCTION:
make_contributor_cache <- function(outputfile = "contributorcache.RData", verbose = TRUE) {
  all.sources <- rotl::tol_about(include_source_list=TRUE)
  all.studies <- as.list(rep(NA,length(all.sources$source_list)))
  author.results <- data.frame()
  curator.results <- data.frame()

  for (i in sequence(length(all.studies))) {
    study.id <- strsplit(all.sources$source_list[[i]], '@')[[1]][[1]]
    all.studies[[i]] <- rotl::get_study_meta(study.id)
    doi <- authors <- clade.name <- curator <- study.year <- NULL
    try(doi <- gsub('http://dx.doi.org/', '', attr(rotl::get_publication(all.studies[[i]]), "DOI")))
    try(authors <- as.character(knitcitations::bib_metadata(doi)$author))

    try(clade.name <- all.studies[i][[1]]$nexml$`^ot:focalCladeOTTTaxonName`)
    try(curators <- unique(rotl::get_study_meta(study.id)[["nexml"]][["^ot:curatorName"]][[1]]))
    if(!is.null(authors)) {
     # author.results <- rbind(author.results, expand.grid(author=authors, study=study.id, clade=clade.name, publication=get_publication(all.studies[[i]])[1], doi=attr(get_publication(all.studies[[i]]),"DOI")))
      author.results <- rbind(author.results, expand.grid(person=authors, study=study.id, clade=clade.name, stringsAsFactors=FALSE))
    }
    if(!is.null(curators)) {
      # author.results <- rbind(author.results, expand.grid(author=authors, study=study.id, clade=clade.name, publication=get_publication(all.studies[[i]])[1], doi=attr(get_publication(all.studies[[i]]),"DOI")))
      curator.results <- rbind(curator.results, expand.grid(person=curators, study=study.id, clade=clade.name, stringsAsFactors=FALSE))
    }
    if(verbose) {
      print(paste(i, "of", length(all.studies), "done"))
    }
  }
  author.pretty <- make_overlap_table(author.results)
  try(Encoding(author.pretty) <- "latin1")
  try(author.pretty <- iconv(author.pretty, "latin1", "UTF-8"))
  curator.pretty <- make_overlap_table(curator.results)
  try(Encoding(curator.pretty) <- "latin1")
  try(curator.pretty <- iconv(curator.pretty, "latin1", "UTF-8"))
  save(author.results, curator.results, author.pretty, curator.pretty, file=outputfile)
  return(list(author.results=author.results, curator.results=curator.results, author.pretty=author.pretty, curator.pretty=curator.pretty))
}
################################ TESTS:
make_contributor_cache <- function(
outputfile = "contributorcache.RData" 
verbose = TRUE
) {
 all.sources <- rotl::tol_about(include_source_list=TRUE)
  all.studies <- as.list(rep(NA,length(all.sources$source_list)))
  author.results <- data.frame()
  curator.results <- data.frame()
i <- 22
i <- 24
i <- 90
i <- 348
  for (i in sequence(length(all.studies))) {
    study.id <- strsplit(all.sources$source_list[[i]], '@')[[1]][[1]]
    all.studies[[i]] <- rotl::get_study_meta(study.id)
    doi <- authors <- clade.name <- curator <- study.year <- NULL
    try(doi <- gsub('http://dx.doi.org/', '', attr(rotl::get_publication(all.studies[[i]]), "DOI")))
    try(authors <- as.character(knitcitations::bib_metadata(doi)$author)) # here's one encoding warning!!
# No encoding supplied: defaulting to UTF-8.
?knitcitations::bib_metadata
length(authors)
	Encoding(authors) <- "latin1"
	authors <- iconv(authors, "latin1", "UTF-8")
	print(authors)
    try(clade.name <- all.studies[i][[1]]$nexml$`^ot:focalCladeOTTTaxonName`)
    try(curators <- unique(rotl::get_study_meta(study.id)[["nexml"]][["^ot:curatorName"]][[1]]))
    if(!is.null(authors)) {
     # author.results <- rbind(author.results, expand.grid(author=authors, study=study.id, clade=clade.name, publication=get_publication(all.studies[[i]])[1], doi=attr(get_publication(all.studies[[i]]),"DOI")))
      author.results <- rbind(author.results, expand.grid(person=authors, study=study.id, clade=clade.name, stringsAsFactors=FALSE))
    }
    if(!is.null(curators)) {
      # author.results <- rbind(author.results, expand.grid(author=authors, study=study.id, clade=clade.name, publication=get_publication(all.studies[[i]])[1], doi=attr(get_publication(all.studies[[i]]),"DOI")))
      curator.results <- rbind(curator.results, expand.grid(person=curators, study=study.id, clade=clade.name, stringsAsFactors=FALSE))
    }
    if(verbose) {
      print(paste(i, "of", length(all.studies), "done"))
    }
  }
  author.results2 <- author.results
  author.pretty <- make_overlap_table(author.results2)
  head(author.results2, 30)
  length(author.pretty)
  names(author.pretty)
  Encoding(author.pretty)
  head(author.pretty)
  levels(author.pretty$person)
  try(Encoding(author.pretty$person) <- "latin1")
  try(author.pretty <- iconv(author.pretty, "latin1", "UTF-8"))
  ?iconv
  curator.pretty <- make_overlap_table(curator.results)
  head(curator.pretty)
  try(Encoding(curator.pretty) <- "latin1")
  try(curator.pretty <- iconv(curator.pretty, "latin1", "UTF-8"))
devtools::load_all()
getwd()
make_contributor_cache(outputfile = '~/Desktop/datelife/data/contributorcache.RData')
1] "345 of 808 done"
[1] "346 of 808 done"
[1] "347 of 808 done"
Error in if (file.exists(x)) { : argument is of length zero
[1] "348 of 808 done"
[1] "349 of 808 done"

# there are non-ASCII encoding errors in object opentree_chronograms, too:
update_datelife_cache <- function(save = TRUE, file = "opentree_chronograms.RData", verbose = TRUE){
	if (save) {
		cache.update <- save_otol_chronograms(file = file, verbose = verbose)
	} else {
		cache.update <- get_otol_chronograms(verbose = verbose)
	}
	return(cache.update)
}

str(opentree_chronograms)
names(opentree_chronograms)
opentree_chronograms$authors
length(opentree_chronograms$authors)
lapply(opentree_chronograms$authors, Encoding)
Encoding(opentree_chronograms$authors[[1]]) <- "latin1"

for (i in sequence(length(opentree_chronograms$authors))){
	print(i)
	print(Encoding(opentree_chronograms$authors[[i]]))
}
length(opentree_chronograms$authors[[172]])
length(opentree_chronograms$curators)
devtools::load_all()
getwd()
update_datelife_cache(save = TRUE, file = "~/Desktop/datelife/data/opentree_chronograms.RData", verbose = TRUE)
[1] "Problematic combos"
[1] "tree_id='tree1', study_id='ot_1091'"    "tree_id='tree1', study_id='ot_1100'"   
[3] "tree_id='Tr95407', study_id='ot_1144'"  "tree_id='tree1', study_id='ot_1150'"   
[5] "tree_id='tree1', study_id='ot_1155'"    "tree_id='Tr48913', study_id='ot_1207'" 
[7] "tree_id='Tr99854', study_id='ot_1208'"  "tree_id='tree1', study_id='ot_209'"    
[9] "tree_id='tree1', study_id='ot_278'"     "tree_id='tree2', study_id='ot_525'"    
[11] "tree_id='tree2', study_id='ot_722'"     "tree_id='tree1', study_id='ot_722'"    
[13] "tree_id='tree1', study_id='ot_755'"     "tree_id='tree2', study_id='ot_806'"    
[15] "tree_id='tree3', study_id='ot_83'"      "tree_id='tree1', study_id='ot_896'"    
[17] "tree_id='tree1', study_id='ot_934'"     "tree_id='tree1', study_id='ot_939'"    
[19] "tree_id='tree1', study_id='ot_945'"     "tree_id='tree1', study_id='ot_952'"    
[21] "tree_id='tree1', study_id='ot_955'"     "tree_id='tree1', study_id='ot_965'"    
[23] "tree_id='tree1', study_id='ot_984'"     "tree_id='tree2', study_id='ot_986'"    
[25] "tree_id='tree2668', study_id='pg_1339'" "tree_id='tree2740', study_id='pg_1372'"
[27] "tree_id='tree3216', study_id='pg_1597'" "tree_id='tree3393', study_id='pg_1684'"
[29] "tree_id='tree3522', study_id='pg_1748'" "tree_id='tree3563', study_id='pg_1766'"
[31] "tree_id='tree3653', study_id='pg_1804'" "tree_id='tree3652', study_id='pg_1804'"
[33] "tree_id='tree3658', study_id='pg_1807'" "tree_id='tree3657', study_id='pg_1807'"
[35] "tree_id='tree3659', study_id='pg_1807'" "tree_id='tree4990', study_id='pg_2374'"
[37] "tree_id='tree4989', study_id='pg_2374'" "tree_id='tree6240', study_id='pg_2688'"
[39] "tree_id='tree6321', study_id='pg_2735'" "tree_id='tree6336', study_id='pg_2740'"
[41] "tree_id='tree5995', study_id='pg_435'" 
Error in `Encoding<-`(`*tmp*`, value = "latin1") : 
  a character vector argument expected

# Continue Runtime tests for:
# make_bold_otol_tree(input = c("Rhea americana",  "Struthio camelus", "Gallus gallus"), usetnrs = FALSE, approximatematch = TRUE, marker = "COI", otol_version = "v2", chronogram = TRUE, doML = FALSE, sppfromtaxon=FALSE, verbose=FALSE)
# # first with doML=FALSE

# finished with 300 spp

# gave following error with 400:
# Error: Request-URI Too Long (HTTP 414)
# In addition: There were 50 or more warnings (use warnings() to see the first 50)
devtools::install_github("phylotastic/datelife")
library(datelife)
setwd("~/Desktop/datelife")
devtools::load_all()
library(ape)
library(geiger)
library(testthat)
library(microbenchmark)
# # # already launched a first run that is not recorded, so cache is already loaded when evaluating the first set.

i <- 400
xname <- paste0("random_sample_",i, "_aves_spp")
setwd("~/Google Drive/datelife/runtime_tests/1_name_samples")
load(file=paste0(xname,".RData"))
make_bold_otol_tree(input=random_sample_400_aves_spp[[1]])
# Error: Request-URI Too Long (HTTP 414)

make_bold_otol_tree <- function(input = c("Rhea americana",  "Struthio camelus", "Gallus gallus"), 
input=random_sample_400_aves_spp[[1]]
use_tnrs = FALSE 
approximate_match = TRUE 
marker = "COI" 
otol_version = "v2" 
chronogram = TRUE 
doML = FALSE 
get_spp_from_taxon = FALSE 
verbose = FALSE
) {
	#otol returns error with missing taxa in v3 of rotl
	input <- input_check(input = input, use_tnrs = use_tnrs, approximate_match = approximate_match, get_spp_from_taxon = get_spp_from_taxon, verbose = verbose)
	input <- input$cleaned_names
	if (verbose) cat("Searching", marker, "sequences for these taxa in BOLD...", "\n")
	sequences <- bold::bold_seqspec(taxon = input, marker = marker)
for(i in 301:400){
	cat("query length:", i, "\t")
	bold::bold_seqspec(taxon = input[1:i], marker = marker)	
	cat("OK", "\n")
}
query length: 301	 OK 
query length: 302	 OK 
query length: 303	 OK 
query length: 304	 OK 
query length: 305	 OK 
query length: 306	 OK 
query length: 307	 OK 
query length: 308	 OK 
query length: 309	 OK 
query length: 310	 OK 
query length: 311	 OK 
query length: 312	 OK 
query length: 313	 OK 
query length: 314	 OK 
query length: 315	 OK 
query length: 316	 OK 
query length: 317	 OK 
query length: 318	 OK 
query length: 319	 OK 
query length: 320	 OK 
query length: 321	 OK 
query length: 322	 OK 
query length: 323	 OK 
query length: 324	 OK 
query length: 325	 OK 
query length: 326	 OK 
query length: 327	 OK 
query length: 328	 OK 
query length: 329	 OK 
query length: 330	 OK 
query length: 331	 OK 
query length: 332	 OK 
query length: 333	 OK 
query length: 334	 OK 
query length: 335	 OK 
query length: 336	 Error: Request-URI Too Long (HTTP 414)

sequences <- bold::bold_seqspec(taxon = input[1:300], marker = marker)
str(sequences)
length(sequences)
lapply(sequences, length)
sequences2 <- sequences
x <- cbind(sequences2, sequences) # cbind just concatenates the two lists, ending up with a 160 items length 
x <- rbind(sequences2, sequences) # this is what I want; concatenating elements within the list
x <- rbind(c(), sequences) 
length(x)
lapply(x, length)
732*2
make_bold_otol_tree(input=random_sample_400_aves_spp[[1]]) # it works now :)
		rr <- suppressWarnings(rotl::tnrs_match_names(names = input[1:250]))
		rr <- suppressWarnings(rotl::tnrs_match_names(names = input[250:300]))
rr

# caught an error while using tnrs_mathc_names, it will onlystick to teh last batch run, so it would analyse 50 names only on 300 input, 150 names in 400, and so on...
For example, went from this:

Phylogenetic tree with 44 tips and 43 internal nodes.

Tip labels:
	Gymnopithys rufigula, Montifringilla taczanowskii, Macronectes giganteus, Eupherusa nigriventris, Stercorarius skua, Calyptorhynchus latirostris, ...
Node labels:
	Aves ott81461, Neognathae ott241846, , , , , ...

Rooted; includes branch lengths.

to this:
Phylogenetic tree with 121 tips and 120 internal nodes.

Tip labels:
	Xenops minutus, Mergellus albellus, Branta canadensis, Egretta thula, Dryocopus lineatus, Calidris tenuirostris, ...
Node labels:
	Aves ott81461, Neognathae ott241846, , , , , ...

Rooted; includes branch lengths.
seq_len(c(1,35))
# restart the test from 300 again...: screen -r 16:07 hrs
ninput <- c(300, 400,500,700,1000,1500,2000,3000,4000, 5000, 6000,7000,8000,9000,10000)
for(i in ninput){
	xname <- paste0("random_sample_",i, "_aves_spp")
	setwd("~/Google Drive/datelife/runtime_tests/1_name_samples")
	load(file=paste0(xname,".RData"))
	y <- microbenchmark(make_bold_otol_tree(input=get(xname)[[1]]), times=1L)
	levels(y$expr)[1] <- as.character(i)
	for(j in 2:100){
		yy <- microbenchmark(make_bold_otol_tree(input=get(xname)[[j]]), times=1L)
		levels(yy$expr)[1] <- as.character(i)
		y <- rbind(y, yy)
	}
	rm(list=xname)
	xnameobj <- paste0("gbot_runtime_2018.02.20_", i,"_aves_spp")
	assign(xnameobj, y)
	setwd("~/Google Drive/datelife/runtime_tests/2_tests/2_random_spp_names/3_gbot")
	save(list=xnameobj, file=paste0(xnameobj,".RData"))
	rm(list=xnameobj)
	print(i)
}
