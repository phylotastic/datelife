# functions to fix negative or null branch lengths and to add dates of missing taxa at random


#' Takes a tree and fixes negative or zero length branches in several ways
#' @param tree A tree either as a newick character string or as a phylo object
#' @param fixing_criterion A character vector specifying the type of branch length to be fixed: "negative" or "zero"
#' @param fixing_method A character vector specifying the method to fix branch lengths: "bladj", "mrbayes" or a number to be assigned to all branches meeting fixing_criterion
#' @return A phylo object with fixed branch lengths
#' @export
tree_fix_brlen <- function(tree = NULL, fixing_criterion = "negative", fixing_method = 0){
	phy <- tree_check(tree = tree)
	fixing_criterion <- match.arg(arg = fixing_criterion, choices = c("negative", "zero"), several.ok = FALSE)
	if(fixing_criterion == "negative"){
		index <- which(phy$edge.length < 0)  # identifies edge numbers with negative edge lengths value
	} else {
		index <- which(phy$edge.length == 0)  # identifies edge numbers with null/zero edge lengths value
	}
	if(!is.numeric(fixing_method)){
		fixing_method <- match.arg(fixing_method, c("bladj", "mrbayes"))
	} else { # chunk for neg or zero br len to zero or any number determined by user
		# stop()
		for (i in index){
			# snode <- pos.phy$edge[i,1]
			# pool  <- pos.phy$edge[seq(nrow(pos.phy$edge))[-i], 1]
			# sisedge <- which(pool==snode) # determines position of sister edge
			# pos.phy$edge.length[sisedge] <- pos.phy$edge.length[sisedge] - pos.phy$edge.length[i]
			# adds neg branch length to sister branch, should add error to both sides???? or only to the daughter branches??
			cnode <- phy$edge[i,2]
			dauedge <- which(phy$edge[,1] == cnode)
			phy$edge.length[dauedge] <- phy$edge.length[dauedge] + phy$edge.length[i] + fixing_method[1]
			phy$edge.length[i] <- fixing_method[1]
			fixed.phy <- phy
		}
	}

	if(any(fixing_method == c("bladj", "mrbayes"))) { #chunk for bladj and mrbayes
		phy <- tree_add_nodelabels(tree = phy)  # all nodes need to be named
		cnode <- phy$edge[index,2]  # we assume that the negative edge length is the one that needs to be changed (but it could be the sister edge that should be shorter)
		tofix <- cnode-length(phy$tip.label)  # so, we use the node number of the crown node of the negative branch lengths
		if(fixing_method == "bladj")
			fixed.phy <- make_bladj_tree(tree = phy,
				nodenames = phy$node.label[-tofix],
				nodeages = tree_get_node_data(tree = phy,
					node_data = "node_age")$node_age[-tofix])
		if(fixing_method == "mrbayes") {
			mrbayes.file <- paste0(deparse(substitute(tree)), "_brlen_fixed.nexus")  # make "phylo" an argument?
			ncalibration <- tree_get_node_data(tree = phy,
				node_data = c("descendant_tips_label", "node_age"))
			ncalibration <- lapply(ncalibration, "[", seq(phy$Nnode)[-tofix])
			phy <- tree_add_outgroup(tree = phy, outgroup = "fake_outgroup")
			fixed.phy <- make_mrbayes_tree(constraint = phy, taxa = phy$tip.label,
				ncalibration = ncalibration, age_distribution = "fixed",
				root_calibration = FALSE,
				missing_taxa = NULL, mrbayes_output_file = mrbayes.file)
			fixed.phy <- ape::drop.tip(fixed.phy, "fake_outgroup")
		}
	}
	return(fixed.phy)
}

#' Takes a tree and uses bladj to estimate node ages and branch lengths given a set of fixed node ages and respective node names
#' @param nodenames A character vector with node names from tree with fixed ages
#' @param nodeages A numeric vector with known or fixed node ages from tree
#' @inheritParams tree_fix_brlen
#' @return A phylo tree
#' @export
make_bladj_tree <- function(tree = NULL, nodenames = NULL, nodeages = NULL){
	phy <- tree_check(tree = tree, dated = FALSE)
	if(is.null(phy$node.label)) {
		stop("phy must have node labels")
	}
	if(!is.null(phy$edge.length)) {
		phy$edge.length <- NULL
	}
	m <- match(nodenames, phy$node.label)
	if(any(is.na(m))) {
		stop("all nodenames must be in phy$node.label") # add a printed line saying which nodenames are not in phy$node.label
	}
	if(length(nodenames) != length(nodeages)) {
		stop("nodenames and nodeages must have the same length")
	}
	if(!is.character(nodenames)) {
		stop("nodenames must be a character vector")
	}
	if(!is.numeric(nodeages)) {
		stop("nodeages must be a numeric vector")
	}
	ages_df <- data.frame(
		a = nodenames,
		b = nodeages
	)
	new.phy <- phylocomr::ph_bladj(ages = ages_df, phylo = phy)
	attributes(new.phy) <- NULL
	new.phy <- phytools::read.newick(text = new.phy)
	# plot(new.phy)
	return(new.phy)
}

#' Takes a constraint tree and uses mrBayes to get node ages and branch lengths given a set of node calibrations without any data.
# we can add the option to use data and no constraint tree.
#' @inheritParams make_mrbayes_runfile
#' @return A phylo tree with branch lengths proportional to time. It will save all mrBayes outputs in the working directory.
#' @export
make_mrbayes_tree <- function(constraint = NULL, taxa = NULL, ncalibration = NULL, missing_taxa = NULL, age_distribution = "fixed", root_calibration = FALSE, mrbayes_output_file = "mrbayes_run.nexus"){
	make_mrbayes_runfile(constraint = constraint, taxa = taxa, ncalibration = ncalibration, age_distribution = age_distribution, root_calibration = root_calibration, missing_taxa = missing_taxa, mrbayes_output_file = mrbayes_output_file)
	message("Running MrBayes. This might take a while.")
	run_mrbayes(mrbayes_output_file = mrbayes_output_file)
	mrbayes_contre <-	tryCatch(ape::read.nexus(paste0(mrbayes_output_file, ".con.tre")), error = function(e) NA)
	message("Done.")
	if(length(mrbayes_contre) == 1){
		if(is.na(mrbayes_contre)){
			stop("MrBayes ran but output files cannot be found.", "\n", "  Please check the log file for errors.")
		}
	}
	return(mrbayes_contre)
}

#' Makes a mrBayes run block file with a constraint topology and a set of node calibrations and missing taxa.
#' @inheritParams get_mrbayes_node_constraints
# add #' @param outgroup = NULL argument
#' @param mrbayes_output_file A character vector specifying the name of mrBayes run file and outputs (can specify directory too).
#' @return A MrBayes block run file in nexus format.
#' @export
make_mrbayes_runfile <- function(constraint = NULL, taxa = NULL, ncalibration = NULL, missing_taxa = NULL, age_distribution = "fixed", root_calibration = FALSE, mrbayes_output_file = "mrbayes_run.nexus"){
  if(!is.null(constraint)) {
  	# constraint <- phylo_tiplabel_space_to_underscore(constraint)
		node_constraints <- get_mrbayes_node_constraints(constraint = constraint, taxa = taxa, ncalibration = ncalibration, age_distribution = age_distribution, root_calibration = root_calibration, missing_taxa = missing_taxa)
		og <- tree_get_singleton_outgroup(tree = constraint) #if(is.null(outgroup))
	} else {
		stop("constraint is NULL")
	}
	ogroup <- c()
	if(!is.na(og)) {
		ogroup <- paste0("outgroup ", og, ";")
	}
	# start mrbayes block:
	bayes_data <- c(paste("   Begin DATA; \nDimensions ntax=", length(node_constraints$taxa), "nchar = 1;"),
	"Format datatype = DNA gap=- missing=?;",
	"Matrix\n",
	paste(node_constraints$taxa, "?"),
	";")
  bayes_set <- c("   Begin MRBAYES;",
  	"unlink shape=(all) tratio=(all) statefreq=(all) revmat=(all) pinvar=(all);\n",
  	node_constraints$node_constraints[-length(node_constraints$node_constraints)],
  	ogroup, "",
  	node_constraints$node_constraints[length(node_constraints$node_constraints)], "",
  	if(!is.null(ncalibration)) "prset nodeagepr = calibrated;", "",
  	node_constraints$node_calibrations, "\n",
  	"   set usebeagle = no Beaglesse = no;", "",
  	paste("prset ", c("brlenspr = clock:birthdeath", "Extinctionpr = Fixed(0)",
  	"Speciationpr = exponential(1)", "clockvarpr = ibr", "ibrvarpr = exponential(10)"), ";", sep=""),
  	"mcmcp nruns = 1 nchains = 1 ngen = 50000000 samplefreq = 1000;",
  	"mcmc;", "",
  	paste0("sumt filename=", mrbayes_output_file, " burnin = 5000000 contype = halfcompat;\n"),
  	"end;"
  )
	all <- c(bayes_data, "\n", bayes_set)
	write(all, mrbayes_output_file)
	return(all)
}

#' Fabricates dates of missing taxa (with no data) on an already dated tree.
#' @param dated_tree a tree (newick or phylo) with branch lengths proportional to absolute time
#' @inheritParams  missing_taxa_check
#' @inheritParams datelife_search
#' @param adding_criterion Only valid when dating_method = "mrbayes". A character vector to specify how missing_taxa should be added to dated_tree.
#' 	 Choose one of:
#' \describe{
#'	\item{adding_method = "random"}{missing_taxa will be added at random to dated_tree.
#'	}
#'	\item{adding_method = "taxonomy"}{taxa will be added to dated_tree following a dataframe with taxonomic assignations given in missing_taxa argument. If no dataframe is given, OToLs reference taxonomy will be used.
#'	}
#'	\item{adding_method = "tree"}{taxa will be added to dated_tree following a tree given in missing_taxa argument. If no tree is given, OToL synthetic tree will be used.
#'	}
#' }
#' @inheritParams make_mrbayes_runfile
#' @return A phylo object
#' @export
tree_add_dates <- function(dated_tree = NULL, missing_taxa = NULL, dating_method = "mrbayes", adding_criterion = "random", mrbayes_output_file = "mrbayes_tree_add_dates.nexus"){
	dated_tree <- tree_check(tree = dated_tree, dated = TRUE)
	missing_taxa <- missing_taxa_check(missing_taxa = missing_taxa, dated_tree = dated_tree)
	dating_method <- match.arg(dating_method, c("bladj", "mrbayes"))
	adding_criterion <- tryCatch(match.arg(adding_criterion, c("random", "taxonomy", "tree")), error = function(e) "random")  # if it does not match any it is assigned to NULL
	if(dating_method == "bladj"){
		# we need to add a missing_taxa check here. We can only use bladj if missing_taxa is a tree
		if(inherits(missing_taxa, "phylo")){
			missing_taxa_phy <- missing_taxa
		} else {
			if(is.data.frame(missing_taxa)){
				all_taxa <- unique(c(dated_tree$tip.label, levels(missing_taxa$taxon)))
			}
			if (is.vector(missing_taxa)){
				all_taxa <- unique(c(dated_tree$tip.label, missing_taxa))
			}
			missing_taxa_phy <- get_otol_synthetic_tree(input = all_taxa)  # tip labes have underscores already
			# this does not always recovers all taxa missing
			# add a warning and a suggestion to rerun tree_add_dates with remaining absent taxa
			# add absent_taxa element here too...
		}
		constraint_tree <- geiger::congruify.phylo(reference = phylo_tiplabel_space_to_underscore(dated_tree), target = missing_taxa_phy, scale = NA)

		dated_tree <- tree_add_nodelabels(tree = dated_tree)  # all nodes need to be named
		new.phy <- make_bladj_tree(tree = missing_taxa, nodenames = dated_tree$node.label, nodeages = tree_get_node_data(tree = dated_tree, node_data = "node_age")$node_age)
	}
	if(dating_method == "mrbayes"){
		dated_tree <- tree_add_outgroup(tree = dated_tree, outgroup = "an_outgroup")  # we need to add a fake outgroup, otherwise mrbayes won't respect the root age
		ncalibration <- tree_get_node_data(tree = dated_tree, node_data = c("node_age", "descendant_tips_label"))
		# we need to be more specific in the way it uses missing taxa next. If it is a tree, then missing_taxa goes as the constraint. If it is a vector, then it just goes as missing taxa. If it is a data frame, we should call pastis.
		new.phy <- make_mrbayes_tree(constraint = dated_tree, ncalibration = ncalibration, missing_taxa = missing_taxa, mrbayes_output_file = mrbayes_output_file)
		new.phy <- ape::drop.tip(new.phy, "an_outgroup")
	}
	return(new.phy)
}
#' Checks that missing_taxa argument is ok to be used by make_mrbayes_runfile inside tree_add_dates functions.
#' @param  missing_taxa A tree, a data frame or a vector enlisting all missing taxa you want to include.
#' \describe{
#'
#' 	\item{A tree.}{Either as a phylo object or as a newick character string.
#'		It contains all taxa that you want at the end, both missing and non missing.
#'		This tree will be used as a hard constraint.
#'	}
#' 	\item{A data frame.}{It contains two columns.
#'		taxon. A character vector of missing taxon names.
#'		clade. A character or numeric vector of nodes from a constraint tree to which each taxon will be assigned.
#'	}
#' 	\item{A character vector.}{It contains the names of the missing taxa.
#'		They will be added at random to the constraint tree.
#'	}
#' }
#'
#' @inheritParams tree_add_dates
#' @return A phylo object, a newick character string or a dataframe with taxonomic assignations
#' @export
missing_taxa_check <- function(missing_taxa = NULL, dated_tree = NULL){
	badformat <- TRUE
	if(is.data.frame(missing_taxa)){ # or is.matrix??
		if("taxon" %in% names(missing_taxa)){
			badformat <- FALSE
			if(length(missing_taxa) > 1 & !"clade" %in% names(missing_taxa)){
				badformat <- TRUE
			}
		}
		# check PastisData format to further check data frame
	} else {
		missing_taxa_phy <- input_process(missing_taxa)  # process input if it is newick or phylo
		if(inherits(missing_taxa_phy, "phylo")){
			phylo_check(phy = dated_tree, dated = TRUE)
			dtINmt <- dated_tree$tip.labels %in% missing_taxa$tip.labels
			mtINdt <- missing_taxa$tip.labels %in% dated_tree$tip.labels
			if (!all(dtINmt)) {
				warning("not all taxa from dated_tree are in missing_taxa tree")
			}
			missing_taxa_pruned <- ape::drop.tip(missing_taxa, missing_taxa$tip.labels[mtINdt])
			# phylo_prune_missing_taxa(phy = , taxa = ) # use this one??
			# dated_tree == missing_taxa_pruned # check that both trees are equal?
			# we don't need to
			# we can congruify if tree are not equal
			# we just need to make a tree with all lineages on dated_tree and all lineages in missing taxa vector
			missing_taxa <- missing_taxa_phy
			badformat <- FALSE
		} else {
			# if missing_taxa is provided as a vector, it must be a character vector:
			# it does not matter if original vector has no real names, if you want to add numbers or NAs or booleans at random, you can
			missing_taxa <- as.character(missing_taxa)
			missing_taxa[which(is.na(missing_taxa))] <- "NA"
			badformat <- FALSE
		}
	}
	if(length(missing_taxa) == 0){ # this is only valid if we don't wanna accept NULL as missing_taxa, but tere are cases in which it is, e.g. tree_fix_brlen
		missing_taxa <- NULL
		badformat <- FALSE
	}
	if(badformat){
		stop("missing_taxa must be a character vector with species names,
		a data frame with taxonomic assignations, a newick character string, a phylo object, or NULL")
	}
	# IMPORTANT: Add a check that taxa in dated.trees is in reference_tree and viceversa
	return(missing_taxa)
	# if (is.null(reference_tree)){
		# if (is.null(add_taxon_distribution)) {
			# cat("specify a reference_tree or add_taxon_distribution to be added to the dated.trees")
			# stop("")
		# } else {
			# # construct a tree with phylotastic or take that from otol?
			# cat("Constructing a reference_tree with taxa from dated.tree and add_taxon_distribution")
		# }
	# } else {
		# if(!is.null(add_taxon_distribution)){
			# cat("A reference_tree was given, add_taxon_distribution argument is ignored")
		# }
}

#' Runs MrBayes from R.
#' @inheritParams make_mrbayes_tree
#' @return A phylo object with the consensus tree. MrBayes output files are stored in the working directory.
#' @export
run_mrbayes <- function(mrbayes_output_file = NULL){
	file <- mrbayes_output_file
	# chunck of code borrowed from phyloch::mrbayes()
	if(is.null(file)) {
		stop("You must provide a block file for MrBayes run")
	}
	if (.Platform$OS.type == "unix"){
		system(paste("mb > execute", file))
	} else {
		system(paste0("mrbayes ", file, ".bayes"))
	}
}

# #' Makes a node calibrations block for a MrBayes run file from a list of taxa and ages or from a dated tree. It can follow a constraint tree.
# #' @inheritParams make_mrbayes_tree
# #' @inheritParams tree_fix_brlen
# #' @param mrbayes_calibration_file NULL or a character vector indicating the name of mrbayes calibration block file.
# #' @param age_distribution A character string specifying the type of calibration. Only "fixed" is implemented for now.
# #' @return A set of MrBayes calibration commands printed in console as character strings or as a text file with name specified in file.
# #' @export
# # This function is set to match node names with constraints obtained from paleotree::GetMrBayesConstraints
# get_mrbayes_node_calibrations <- function(constraint = NULL, ncalibration = NULL, age_distribution = "fixed", 	mrbayes_calibration_file = NULL){
# 	if(!is.null(constraint)){
# 		phy <- tree_check(tree = constraint)
# 	}
# 	if(length(ncalibration) == 2){  # if it is a list of descendant tips labels and node ages, from tree_get_node_data function
# 		if(!is.list(ncalibration)) {
# 			stop("ncalibration must be a newick character string, a phylo object or a list with taxon names and dates")
# 		}
# 		includes.ncalibration <- lapply(ncalibration$descendant_tips_label, function(x) gsub(" ", "_", x))
# 		nages <- ncalibration$node_age
#
# 	} else {  # if it is a tree
# 			ncalibration <- tree_check(tree = ncalibration, dated = TRUE)
# 	    phy <- phylo_tiplabel_space_to_underscore(phy)
# 	    ncalibration <- phylo_tiplabel_space_to_underscore(ncalibration)
# 	    splits.ncalibration <- ape::prop.part(ncalibration)
# 	    includes.ncalibration <- lapply(splits.ncalibration, function(x) ncalibration$tip.label[x])
# 			nages <- ape::branching.times(ncalibration)
# 	}
# 	nodes <- sapply(includes.ncalibration, function(tax)
# 				phytools::findMRCA(phy, tax, type="node")) - length(phy$tip.label)
# 	calibrations <- paste0("calibrate node", nodes, " = ", age_distribution, "(", nages, ");")
# 	root <- which(nodes == 1)  # tests for the presence of a root calibration, which should be implemented with treeagepr and not with calibrate
# 	if(length(root) != 0){
# 		nodes <- nodes[-root]
# 		nages <- nages[-root]
# 		calibrations <- c(calibrations, paste0("prset treeagepr = ", age_distribution, "(",
# 		nages[root], ");"))
# 	}
#     if (!is.null(mrbayes_calibration_file)) {
#         write(calibrations, mrbayes_calibration_file)
#     }
#     else {
#         return(calibrations)
#     }
# }

#' Makes a node constraints and node calibrations block for a MrBayes run file from a list of taxa and ages or from a dated tree
#' @param constraint The constraint tree: a phylo object or a newick character string, with or without branch lengths.
#' @param taxa A character vector with taxon names to be mantained in tree
#' @inheritParams missing_taxa_check
#' @param ncalibration The node calibrations: a phylo object with branch lengths proportional to time; in this case all nodes from ncalibration will be used as calibration points. Alternatively, a list with two elements: the first is a character vector with node names from phy to calibrate; the second is a numeric vector with the corresponding ages to use as calibrations.
#' @param age_distribution A character string specifying the type of calibration. Only "fixed" and "uniform" are implemented for now.
#' \describe{
#'
#' 	\item{fixed}{ The age given in ncalibration will be used as fixed age.
#'	}
#' 	\item{lognormal}{The age given in ncalibration will be used as mean age.
#' 		The standard deviation can be provided. # still need to add this option.
#' 		By default, a 95 CI sd is used.
#'	}
#' 	\item{uniform}{The age given in ncalibration will be used as mean age.
#' 		Where min_age = 0.9 * mean age, and max_age = 1.1 * mean age.
#'	}
#' }
#'
#' @param root_calibration Used to set a calibration at the root or not. Default to FALSE. Only relevant if ncalibration is specified.
#' @param mrbayes_constraints_file NULL or a character vector indicating the name of mrbayes costraint and/or calibration block file.
#' @param clockratepr A character vector indicating the clockrateprior to be used.
#' @return A set of MrBayes constraints and/or calibration commands printed in console as character strings or as a text file specified in mrbayes_constraints_file.
#' @export
get_mrbayes_node_constraints <- function(constraint = NULL, taxa = NULL, missing_taxa = NULL, ncalibration = NULL, age_distribution = "fixed", root_calibration = FALSE, mrbayes_constraints_file = NULL, clockratepr = "prset clockratepr = fixed(1);"){
	stop_flag <- TRUE
  if(is.list(constraint) & "descendant_tips_label" %in% names(constraint)){
    stop_flag <- FALSE
  } else {
    constraint <- try(tree_check(tree = constraint), silent = TRUE)
    if(inherits(constraint, "phylo")) {  # if it is a tree and we want all its nodes to be used as calibrations over a constraint tree
			stop_flag <- FALSE
      constraint <- tree_get_node_data(tree = constraint, node_data = c("descendant_tips_label"))
      # constraint <- # if there are problems with setting a constraint at the root node, we should eliminate it. But check it, since it might not affect results
    }
  }
  if(stop_flag){
    stop("constraint must be a newick character string, a phylo object or a list of taxon names. See details.")
  }
	# all names with underscores:
  constraint$descendant_tips_label <- sapply(constraint$descendant_tips_label, gsub, pattern = " ", replacement = "_")
	if(!is.null(taxa)){
		taxa <- gsub(" ", "_", taxa)
	}
	# make sure that all taxa are included:
	taxa <- unique(c(taxa, unique(unlist(constraint$descendant_tips_label))))  # if taxa is NULL it will keep names from constraint OK
	missing_taxa <- missing_taxa_check(missing_taxa)
	if(!is.null(missing_taxa)){
		# for now missing taxa is only accepted as a vector:
		if(is.vector(missing_taxa)){
			missing_taxa <- gsub(" ", "_", missing_taxa)
			taxa <- unique(c(taxa, missing_taxa))
		}
	}
	# remove a constraint including all taxa, which is not allowed by mrbayes
	# only remove the lines that have all taxa in them or names not included in taxa:
	taxaINconst <- lapply(constraint$descendant_tips_label, match, taxa) # match names in taxa over constraints
	# we should think about the following line; maybe we should keep taxon names that are in constraints but not in taxa
	# I think we already make sure above that we have all taxon names in constraints in taxa vector; it is impossible to obtain the following scenario then:
	remove <- which(sapply(taxaINconst, anyNA))  # removes constraints with taxon names in constraints but not in taxa
	remove <- c(remove, which(length(taxa) == sapply(taxaINconst, length)))  # keep constraints that contain less than all taxa
	if(length(remove) > 0) {
	  constraint$descendant_tips_label <- constraint$descendant_tips_label[-remove]
	}
	# text block:
  node_constraints <- paste0("constraint node",
                             names(constraint$descendant_tips_label),
                             " = ",
                             sapply(constraint$descendant_tips_label, function(x) paste0(x, collapse = " ")),
                             ";")

  node_constraints <- c(node_constraints, " ",
	          				paste0("prset topologypr = constraints(",
										paste0("node", names(constraint$descendant_tips_label), collapse = ", "),
										");"))
  node_constraints <- list(node_constraints = node_constraints)

  if(!is.null(ncalibration)){
    stop_flag <- TRUE
    if(is.list(ncalibration) & "node_age" %in% names(ncalibration) &
       "descendant_tips_label" %in% names(ncalibration)){  # if it is a list of descendant tips labels and node ages, from tree_get_node_data function
      stop_flag <- FALSE
    } else {
      ncalibration <- try(tree_check(tree = ncalibration, dated = TRUE), silent = TRUE)
      if(inherits(ncalibration, "phylo")) {  # if it is a tree and we want all its nodes to be used as calibrations over a constraint tree
        ncalibration <- tree_get_node_data(tree = ncalibration, node_data = c("node_age", "descendant_tips_label"))
        stop_flag <- FALSE
      }
    }
    if(stop_flag){
      stop("ncalibration must be a newick character string or a phylo object with branch lengths, or a list of two elements containing taxon names and ages. See details.")
    }
		# all names with underscores too:
    ncalibration$descendant_tips_label <- sapply(ncalibration$descendant_tips_label, gsub, pattern = " ", replacement = "_")
    index <- match(constraint$descendant_tips_label, ncalibration$descendant_tips_label)  # this is useful to match ncalibration nodes to constraint nodes, it gets the order of ncalibration nodes
		index <- index[!is.na(index)]  # removes NA that appear when a constraint has no calibration
		if(length(index) != length(ncalibration$descendant_tips_label)){
      		warning("not all calibrations are present as constraints")
    	} else {
			message("all calibrations are present as constraints")
		}
		age_distribution <- match.arg(arg = age_distribution, choices = c("fixed", "lognormal", "gamma", "uniform"), several.ok = FALSE)
		if(age_distribution == "fixed"){
			age_distribution_set <- paste0(age_distribution,
															"(",
															ncalibration$node_age[index],  # we make ncalibration nodes match constraint nodes
															");")
		}
		if(age_distribution == "lognormal"){
			# tried changing to mean_log, but that is not the problem with mrbayes
			# meanlog <- function(age_mu, age_var){
			# 	log(age_mu / sqrt(1 + (age_var / age_mu^2)))
			# }
			# sdlog <- function(age_mu, age_var){
			# 	sqrt(log(1 + age_var / age_mu^2))
			# }
			# mean_log <- meanlog(ncalibration$node_age[index], age_var = 1)
			# sd_log <- sdlog(ncalibration$node_age[index], age_var = 1)
			# exp(mean_log + sd_log^2/2)
			# sqrt((exp(sd_log^2) - 1) * (exp(2 * mean_log +sd_log^2)))
			age_distribution_set <- paste0(age_distribution,
															"(",
															ncalibration$node_age[index],  # ncalibration nodes match constraint nodes
															",",
															ncalibration$node_age[index] * 0.1,  # sd is set to 10% of mean age
															");")
		}
		if(age_distribution == "uniform"){
			age_distribution_set <- paste0(age_distribution,
															"(",
															ncalibration$node_age[index] * 0.9,  # index makes ncalibration nodes match constraint nodes
															",",
															ncalibration$node_age[index] * 1.1,
															");")
		}

		if(age_distribution == "lognormal_uniform"){
			my_rlnorm <- function(n, age_mu, age_var){
        res <- stats::rlnorm(n = n, meanlog = log(age_mu / sqrt(1 + (age_var / age_mu^2))), sdlog = sqrt(log(1 + age_var / age_mu^2)))
        return(res)
      }
			min_max_age <- lapply(ncalibration$node_age[index],
				function(x) {
					y <- my_rlnorm(n = 100, age_mu = x, age_var = x * 0.1)
					return(c(min(y), max(y)))
				}
			)
			age_distribution_set <- paste0(age_distribution,
															"(",
															sapply(min_max_age, "[", 1),  # index makes ncalibration nodes match constraint nodes
															",",
															sapply(min_max_age, "[", 2),
															");")
		}

    node_calibrations <- paste0("calibrate node",
                                names(constraint$descendant_tips_label[index]),  # we maintain the nodes from constraint that have calibrations in index
                                " = ",
                                age_distribution_set)

		if(is.character(clockratepr)) {
			node_calibrations <- c(node_calibrations, clockratepr)
		}

    if(root_calibration){
      node_calibrations <- c(node_calibrations, paste0("prset treeagepr = ", age_distribution, "(",
                              max(ncalibration$node_age), ");"))
    } else if (is.numeric(root_calibration)){
      node_calibrations <- c(node_calibrations, paste0("prset treeagepr = ", age_distribution, "(",
                            root_calibration, ");"))
    }
    node_constraints <- c(node_constraints, list(node_calibrations = node_calibrations))
  }
  # nodes <- sapply(includes.ncalibration, function(tax)
    # phytools::findMRCA(phy, tax, type="node")) - length(phy$tip.label)

  if (!is.null(mrbayes_constraints_file)) {
    for (i in 1:length(node_constraints)) {
      write(node_constraints[[i]], sep = "\n", file = mrbayes_constraints_file, append = TRUE)
    }
  }
	node_constraints$taxa <- taxa
  return(node_constraints)
}

#' Identifies the presence of a single lineage outgroup in a phylogeny.
#' @inheritParams make_mrbayes_tree
#' @inheritParams tree_fix_brlen
#' @return A character vector with the name of the single lineage outgroup. Returns NA if there is none.
#' @export
tree_get_singleton_outgroup <- function(tree = NULL){
  phy <- tree_check(tree = tree)
	phy <- phylo_tiplabel_space_to_underscore(phy)
    outgroup <- NA
    splits <- ape::prop.part(phy)
    if(length(splits)>1){
    	index <- which.max(sapply(splits, length))
    	s1 <- splits[[index]]
    	index2 <- which.max(sapply(splits[-index], length))
    	s2 <- splits[-index][[index2]]
    	if(length(s1)-length(s2) == 1){
    		outgroup <- phy$tip.label[s1[!(s1 %in% s2)]]
    	}
    }
    return(outgroup)
}
