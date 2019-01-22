# Viernes 12 de Enero 2018
# Runtime tests for:
# make_bold_otol_tree(input = c("Rhea americana",  "Struthio camelus", "Gallus gallus"), usetnrs = FALSE, approximatematch = TRUE, marker = "COI", otol_version = "v2", chronogram = TRUE, doML = FALSE, sppfromtaxon=FALSE, verbose=FALSE)
# # first with doML=FALSE
# summarizeResults
# make_datelife_query √
# evaluate cache loading? with different amount of trees?
# evaluate usetnrs functions

devtools::install_github("phylotastic/datelife")
library(datelife)
source('~/Desktop/datelife/R/datelife.R')
library(ape)
library(geiger)
library(testthat)
library(microbenchmark)
# Warning message:
# In strptime(x, fmt, tz = "GMT") :
  # unknown timezone 'default/America/New_York'
# # # Use Randomized spp names for each run (100 times)
# # # Launch a first run that is not recorded, so cache is already loaded when evaluating the first set.
asd <- "((Zea mays,Oryza sativa),((Arabidopsis thaliana,(Glycine max,Medicago sativa)),Solanum lycopersicum)Pentapetalae);"
asd.pi <- make_datelife_query(asd)
asd.names <- asd.pi$cleaned.names
make_bold_otol_tree(input=asd.names)
get_datelife_result(asd, process_input=TRUE)
# now start the test runs:
source('~/Desktop/datelife/R/datelife.R', chdir = TRUE)
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
	xnameobj <- paste0("gbot_runtime_2018.01.12_", i,"_aves_spp")
	assign(xnameobj, y)
	setwd("~/Google Drive/datelife/runtime_tests/2_tests/2_random_spp_names/3_gbot")
	save(list=xnameobj, file=paste0(xnameobj,".RData"))
	rm(list=xnameobj)
}
# screen -r 6417
# Error in i <- 10 and j <- 3
make_bold_otol_tree(input=get(xname)[[j]])
# Error in data[, tree$tip.label] : subscript out of bounds
	input=get(xname)[[j]]
	input_check(input=get(xname)[[j]])  # this is not the problem...OK
	row.index <- 0
	taxa.to.drop <- c()
	for (i in input){
		row.index <- row.index + 1
		taxon.index <- which(grepl(i, sequences$species_name)) # what happens here if there are no sequences from the taxon????
		if (length(taxon.index)>0){
			seq.index <- which.max(sequences$nucleotide_ATGC_length[taxon.index])
			# sequences[taxon.index,][seq.index,]
			seq <- strsplit(sequences$nucleotides[taxon.index][seq.index], split = "")[[1]]
			final.sequences[row.index, sequence(length(seq))] <- seq
		} else {
			taxa.to.drop <- c(taxa.to.drop, i)
		}
		final.sequences.names[row.index] <- i
	}
	rownames(final.sequences) <- gsub(" ", "_", final.sequences.names)

phy$tip.label
names(alignment)
phy$tip.label %in% names(alignment)
final.sequences.names <- gsub(" ", "_", final.sequences.names)
phy$tip.label %in% final.sequences.names
phy$tip.label[4]  # this lineage is not in alignment matrix, why??????????
# [1] "Charmosyna_aureicincta"
# fixed this in file datelife_practice_2018.01.14, yay!

# run again:
# screen -r
# another error:
# Error in make_bold_otol_tree(input = get(xname)[[j]]) :
  # Not enough sequence coverage in BOLD to perform analysis of this set of taxa
# fixed it...

# run 3rd time:
# screen -r
# Error in check_ott_ids(ott_ids) : NAs are not allowed
# In addition: Warning messages:
# 1: In make_bold_otol_tree(input = get(xname)[[j]]) :
  # Not enough sequences available in BOLD. No tree was constructed.
# 2: In check_tnrs(res) :
  # Mesopicos elliotii, Primobucco olsoni, Motacilla feldegg are not matched
# 3: In rotl::tnrs_match_names(names = input) : NAs introduced by coercion

# run 4th time:
# Error: $ operator is invalid for atomic vectors

# run 5th time:
# Error in data[, tree$tip.label] : subscript out of bounds

# run 6th time:
# Error: HTTP failure: 400
# Queries containing more than 250 strings are not supported. You may submit multiplesmaller queries to avoid this limit.BadInputExceptionorg.neo4j.server.rest.repr.BadInputExceptionlist("org.opentree.taxonomy.plugins.tnrs_v3.match_names(tnrs_v3.java:241)", "java.lang.reflect.Method.invoke(Method.java:498)", "org.neo4j.server.plugins.PluginMethod.invoke(PluginMethod.java:57)", "org.neo4j.server.plugins.PluginManager.invoke(PluginManager.java:168)", "org.neo4j.server.rest.web.ExtensionService.invokeGraphDatabaseExtension(ExtensionService.java:300)", "org.neo4j.server.rest.web.ExtensionService.invokeGraphDatabaseExtension(ExtensionService.java:122)", "java.lang.reflect.Method.invoke(Method.java:498)",
#     "org.neo4j.server.rest.security.SecurityFilter.doFilter(SecurityFilter.java:112)")
# In addition: There were 50 or more warnings (use warnings() to see the first 50)


load('~/Google Drive/datelife/runtime_tests/aves.spp.RData')
length(aves.spp$cleaned.names)
yy <- microbenchmark(make_bold_otol_tree(input=aves.spp$cleaned.names), times=100L)
levels(yy$expr)[1] <- "12750"
assign("gbot_runtime_2018.01.12_12750_aves_spp", yy)
setwd("~/Google Drive/datelife/runtime_tests/2_tests/0_all_names")
save(gbot_runtime_2018.01.12_12750_aves_spp, file="gbot_runtime_2018.01.12_12750_aves_spp.Rdata")

# # # make a plot:
setwd("~/Google Drive/datelife/runtime_tests/2_tests/2_random_spp_names/3_gbot")
ninput <- c(10,100,200,300,400,500,700,1000,1500,2000,3000,4000, 5000, 6000,7000,8000,9000,10000)
length(ninput)
res <- c()
for(i in ninput){
	xname <- paste0("gbot_runtime_2018.01.12_",i,"_aves_spp")
	x <- paste0(xname, ".RData")
	load(x)
	res <- rbind(res, get(xname))
}
length(res)
str(res)
load('~/Google Drive/datelife/runtime_tests/2_tests/0_all_names/gbot_runtime_2018.01.12_12750_aves_spp.Rdata')
gbot_runtime_2018.01.12_12750_aves_spp
# levels(gbot_runtime_2018.01.12_12750_aves_spp$expr) <- "all Aves (12750)"
res <- rbind(res, gbot_runtime_2018.01.12_12750_aves_spp)
ggplot2::autoplot(res)
setwd("~/Google Drive/datelife/runtime_tests/2_tests/2_random_spp_names/3_gbot")
pdf(file="make_bold_otol_tree_microbenchmark_2018.01.12_xtime_good_labels.pdf")
y_min <- 200
y_max <- 1e+5
    res$Time <- microbenchmark:::convert_to_unit(res$time, "t") #changing the name of the element itself is the easiest way to make it appear as axis label
    # object$'Query Length' <- object$expr #changing for a name with spaces won't work...
    plt <- ggplot2::ggplot(res, ggplot2::aes_string(x = "expr", y = "Time"))
    plt <- plt + ggplot2::coord_cartesian(ylim = c(y_min, y_max))
    plt <- plt + ggplot2::stat_ydensity()
    plt <- plt + ggplot2::scale_x_discrete(name = "Query Length")
    plt <- plt + ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1))
    plt <- plt + ggplot2::theme(axis.text.y = ggplot2::element_text(angle=0))
	plt <- 	plt + ggplot2::scale_y_log10(name="Time (seconds)", breaks=c(1e+03, 1e+035, 1e+04, 1e+045, 1e+05), labels=c("1e+03"="1 s", "1e+035"="", "1e+04"="10 s", "1e+045"="", "1e+05"="100 s"), position="left")
    plt
dev.off()
