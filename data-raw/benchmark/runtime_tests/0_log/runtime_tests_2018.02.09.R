# Viernes 9 de Febrero 2018
# Continue Runtime tests for:
# make_bold_otol_tree(input = c("Rhea americana",  "Struthio camelus", "Gallus gallus"), usetnrs = FALSE, approximatematch = TRUE, marker = "COI", otol_version = "v2", chronogram = TRUE, doML = FALSE, sppfromtaxon=FALSE, verbose=FALSE)
# # first with doML=FALSE

devtools::install_github("phylotastic/datelife")
library(datelife)
devtools::load_all()
library(ape)
library(geiger)
library(testthat)
library(microbenchmark)

# Fix error in make_bold_otol_tree tests, when analyzing more than 250 names:
Error: HTTP failure: 400
Queries containing more than 250 strings are not supported. You may submit multiplesmaller queries to avoid this limit.BadInputExceptionorg.neo4j.server.rest.repr.BadInputExceptionlist("org.opentree.taxonomy.plugins.tnrs_v3.match_names(tnrs_v3.java:241)", "java.lang.reflect.Method.invoke(Method.java:498)", "org.neo4j.server.plugins.PluginMethod.invoke(PluginMethod.java:57)", "org.neo4j.server.plugins.PluginManager.invoke(PluginManager.java:168)", "org.neo4j.server.rest.web.ExtensionService.invokeGraphDatabaseExtension(ExtensionService.java:300)", "org.neo4j.server.rest.web.ExtensionService.invokeGraphDatabaseExtension(ExtensionService.java:122)", "java.lang.reflect.Method.invoke(Method.java:498)",
    "org.neo4j.server.rest.security.SecurityFilter.doFilter(SecurityFilter.java:112)")
In addition: There were 50 or more warnings (use warnings() to see the first 50)
> xnameobj
[1] "gbot_runtime_2018.01.12_200_aves_spp"
> xbname
Error: object 'xbname' not found
> xname
[1] "random_sample_300_aves_spp"


# # # Launch a first run that is not recorded, so cache is already loaded when evaluating the first set.
asd <- "((Zea mays,Oryza sativa),((Arabidopsis thaliana,(Glycine max,Medicago sativa)),Solanum lycopersicum)Pentapetalae);"
asd.pi <- make_datelife_query(asd)
asd.names <- asd.pi$cleaned.names
make_bold_otol_tree(input=asd.names)
get_datelife_result(asd)
i <- 300
xname <- paste0("random_sample_",i, "_aves_spp")
setwd("~/Google Drive/datelife/runtime_tests/1_name_samples")
load(file=paste0(xname,".RData"))
make_bold_otol_tree(input=random_sample_300_aves_spp[[1]])
# Warning messages:
# 1: In check_tnrs(res) : Chrysococcyx rufomerus are not matched
# 2: In rotl::tnrs_match_names(names = input[xx[i]:yy[i]]) :
#   NAs introduced by coercion
# 3: In make_bold_otol_tree(input = random_sample_300_aves_spp[[1]]) :
#   	Negative branch lengths in BOLD chronogram.
rotl::tnrs_match_names(names = "Chrysococcyx rufomerus")
# now start the test runs from 300 names:
ninput <- c(10,100,200,300,400,500,700,1000,1500,2000,3000,4000, 5000, 6000,7000,8000,9000,10000)
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
	xnameobj <- paste0("gbot_runtime_2018.02.09_", i,"_aves_spp")
	assign(xnameobj, y)
	setwd("~/Google Drive/datelife/runtime_tests/2_tests/2_random_spp_names/3_gbot")
	save(list=xnameobj, file=paste0(xnameobj,".RData"))
	rm(list=xnameobj)
    print(i)
}
