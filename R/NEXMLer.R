## write NEXML version of a tree to file - quick ass-ugly version
## takes in object of either phylo or multiPhylo
## adapted from APE's write.nexus
## embarassingly ugly: like write.nexus does not create a nexus object, NEXMLer does not
##   create an xml object, but just prints out elements of phylo object to file
##   - an advantage (for large trees) is that new variables need not be created
##      - i.e. save memory
##   - ideally an xml object would be created. for later.
## assumption is that all trees have the same number and specific tip names
## tree(s) need not be binary (but see below. for now.)

# i'm working from an example nexml file which doesn't seem to use referencing...
# distinct possibility i'm fucking up id vs. label...

## TODO: have not yet considered case where:
##   1) root edge exists
##   2) no edge lengths exist
## also have not dealt with meta data, but seems straightforward

DEBUG <- FALSE;

NEXMLer <- function(phy, file = "") {
	if (class(phy) != "phylo" && class(phy) != "multiPhylo") {
		stop("Wrong input, fucktard.\n");
	}
	
	nTrees <- 0L;
	nTaxa <- 0L;
	
# the rooting state of all trees should(?) be identical. may be efficient...
# also assuming constant number of taxa
	rootage <- NULL;
	obj <- list(phy); # this is the ape strategy; force into list and determine nTrees later
	
	if (class(phy) == "phylo") {
		nTrees <- 1;
		nTaxa <- length(phy$tip.label);
		if (is.rooted(phy)) {
			rootage <- "rooted";
		} else {
			rootage <- "unrooted";
		}
	} else { # dealing with multiple trees
		nTrees <- length(phy);
		obj <- obj[[1]];
		if (all(unlist(lapply(phy, is.rooted)))) {
			rootage <- "rooted";
		} else if (all(!unlist(lapply(phy, is.rooted)))) {
			rootage <- "unrooted";
		} else {
			rootage <- "variable";
		}
		nTaxa <- length(phy[[1]]$tip.label);
	}
	
	if (DEBUG) cat("\nRootage=", rootage, "\n\n", sep="");
	
	
## store spacing as a variable; or create true xml object
	
	cat("<?xml version=\"1.0\"?>\n", file = file);
	
# header: move this bit <outside> for nicer updating
	cat(paste("<nex:nexml\n",
	"\tgenerator=\"NEXMLer\"\n",
	"\tversion=\"0.9\"\n",
	"\txml:base=\"http://example.org/\"\n",
	"\txmlns:nex=\"http://www.nexml.org/2009\"\n",
	"\txmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns\"\n",
	"\txmlns:xml=\"http://www.w3.org/XML/1998/namespace\"\n",
	"\txmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n",
	"\txsi:schemaLocation=\"http://www.nexml.org/2009 ../xsd/nexml.xsd\"\n",
	"\txmlns=\"http://www.nexml.org/2009\">\n"), file = file, append = TRUE);
# </outside>
	
	cat("\t<otus id=\"OTUs\" label=\"fooTaxaBlock\">\n", file = file, append = TRUE);
	for (i in 1:nTaxa) {
		cat(paste("\t\t<otu id=\"", obj[[1]]$tip.label[i], "\"/>\n", sep=""), file = file, append = TRUE);
	}
	cat("\t</otus>\n", file = file, append = TRUE);
	
	cat("\t<trees otus=\"OTUs\" id=\"Trees\" label=\"TreesBlockFromXML\">\n", file = file, append = TRUE);
	
	for (i in 1:nTrees) {
		treeLabel <- paste("tree", i, sep="");
	# label is optional, although could be passed in via multiPhylo; make flexible (later)
		cat(paste("\t\t<tree id=\"", treeLabel, "\" xsi:type=\"nex:FloatTree\"",
			" label=\"", treeLabel, "\">\n", sep=""), file = file, append = TRUE);
		
		tree <- obj[[i]] # ugh. ugly, but safe(?)
		tips <- tree$tip.label;
		edges <- tree$edge;
		lengths <- tree$edge.length;
		numNode <- 0L;
		
	# NOTE: this here implicitly assumes binary trees; make more general
		if (rootage == "rooted") {
			numNode <- (2 * nTaxa) - 1;
		} else if (rootage == "unrooted") {
			numNode <- (2 * nTaxa) - 2;
		} else { # rooting is apparently variable. check.
			numNode <- ifelse(is.rooted, (2 * nTaxa) - 1, (2 * nTaxa) - 2);
		}
		
		for (j in 1:numNode) {
			nodeLabel <- paste("n", j, sep="");
			cat(paste("\t\t\t<node id=\"", nodeLabel, "\" label=\"", nodeLabel, "\"", sep=""));
			if (j <= nTaxa) {
				cat(paste(" otu=\"", tips[j], "\"", sep=""));
			} else if (j == (nTaxa + 1) && rootage == "rooted") {
				cat(" root=\"true\"");
			}
			cat("/>\n", sep="");
		}
		
		for (j in 1:length(edges[,1])) {
			cat(paste("\t\t\t<edge source=\"n", edges[j,1], "\" target=\"n", edges[j,2],
				"\" id=\"e", j, "\" length=\"", lengths[j], "\"/>\n", sep=""));
		}
		cat("\t\t</tree>\n", file = file, append = TRUE);
	}
	
	cat("\t</trees>\n", file = file, append = TRUE);
	
	cat("</nex:nexml>", file = file, append = TRUE);
}