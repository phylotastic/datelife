#' Chronogram database
#'
#' Now storing >200 chronograms from Open Tree of Life
#'
#' @name opentree_chronograms
#' @docType data
#' @format A list of four elements, containing data from OpenTree of Life chronograms
#' \describe{
#'   \item{authors}{A list of lists of author names of the original studies that
#'     published chronograms in the Open Tree of Life database.}
#'   \item{curators}{A list of lists of curator names that uplodaded chronograms
#'     to the Open Tree of Life database.}
#'   \item{studies}{A list of study identifiers.}
#'   \item{trees}{A `multiPhylo` object storing the chronograms from Open Tree of
#'     Life database.}
#' }
#' @source \url{http://opentreeoflife.org}
#' @keywords opentree dates myrs million years time phylogeny chronogram
#' @details
#' Generated with
#' opentree_chronograms <- get_otol_chronograms()
#' usethis::use_data(opentree_chronograms, overwrite = T)
#' and updated with update_datelife_cache()
"opentree_chronograms"


#' Create an updated OpenTree chronograms database object

#' @description The function calls [get_otol_chronograms()] to update the OpenTree
#'  chronograms database cached in datelife. It has the option to write the updated
#' object as an .Rdata file, that will be independent of the `opentree_chronograms`
#' data object that you can load with `data("opentree_chronograms", package = "datelife")`.
#' @param write Defaults to `TRUE`, it saves an .Rdata file named indicated by argument `name`,
#'   containing available chronograms from Open Tree of Life. Saves to path indicated by argument `path`.
#' @param updated_name Used if `write = TRUE`. Defaults to `"opentree_chronograms_updated"`. A character
#'   vector of length one indicating the name to assign to both the updated OpenTree chronogram
#'   database object and the ".Rdata" file. For example, if `name = "my_database"`, the
#'   function will assign the updated chronogram database to an object named `my_database`
#'   and will write it to a file named "my_database.Rdata" in the path indicated
#'   by argument `file_path`.
#' @param file_path Used if `write = TRUE`. A character vector of length 1 indicating
#'   the path to write the updated database ".Rdata" file to, excluding file name.
#'   Defaults to temporary directory obtained with [base::tempdir()] and formatted with
#'   [base::file.path()].
#' @inheritDotParams get_otol_chronograms
#' @inherit get_otol_chronograms return
#' @export
update_datelife_cache <- function(write = TRUE,
                                  updated_name = "opentree_chronograms_updated",
                                  file_path = file.path(tempdir()),
                                  ...) {
  # , new_studies_only = TRUE
  # enhance: I think we can change the name to update_opentree_chronograms
  updated <- get_otol_chronograms(...)
  if (write) {
    assign(updated_name, updated) # assign the name indicated in updated_name to updated object
    file_name <- paste0(updated_name, ".RData")
    save(list = updated_name, file = file.path(file_path, file_name), compress = "xz")
    message("\nUpdated opentree_chronograms object was saved to file:\n\t", file_name, "\n\nIn directory:\n\t", file_path)
  }
  return(updated)
}

#' Update all data files as data objects for the package
#'
#' This includes opentree chronograms, contributors, treebase and curators
#' For speed, datelife caches chronograms and other information. Running this (within the checked out version of datelife) will refresh these. Then git commit and git push them back
#' @return None
#' @export
update_all_cached <- function() {
  opentree_chronograms <- get_otol_chronograms()
  usethis::use_data(opentree_chronograms, overwrite = TRUE)
  contributor_cache <- make_contributor_cache(outputfile = paste0(tempdir(), "/contributor.rda"))
  usethis::use_data(contributor_cache, overwrite = TRUE)
  treebase_cache <- make_treebase_cache(outputfile = paste0(tempdir(), "/treebase.rda"))
  usethis::use_data(treebase_cache, overwrite = TRUE)
  depositor_cache <- make_all_associations(outputfile = NULL)
  usethis::use_data(depositor_cache, overwrite = TRUE)
}

#' Check if a tree has branch lengths
#' @inheritParams phylo_check
#' @return A TRUE or FALSE
#' @export
phylo_has_brlen <- function(phy) {
  brlen <- TRUE
  if (is.null(phy$edge.length)) {
    brlen <- FALSE
  }
  return(brlen)
}


# #' Allow trees with duplicate taxa to be read in, courtesy David Winter
# #' @param study_id Open Tree study id
# #' @param tree_id Open Tree tree id
# #' @param tip_label Which Open Tree tree label format you want
# #' @return A phylo object
# get_study_tree(study_id="pg_710", tree_id="tree1277",
#                tip_label='ott_taxon_name', file = "/tmp/tree.tre",
#                file_format = "newick")
# tr <- ape::read.tree(file = "/tmp/tree.tre")
# study_id <- "ot_1207"
# tree_id <- "Tr48913"
# get_study_tree_with_dups <- function(study_id, tree_id, tip_label="ott_taxon_name") {
# 	tr <- rotl::get_study_tree(study_id=study_id, tree_id=tree_id, tip_label=tip_label,
# 		file_format="newick", file = "/data-raw/tree.tre")
# 	phy <- ape::read.tree(text=gsub(" ", "_", tr))
# 	phy$tip.label <- gsub("'", "", gsub("_", " ", phy$tip.label))
# 	return(	phy)
# }

#' Get all chronograms from Open Tree of Life database
#' @param max_tree_count Default to "all", it gets all available chronograms. For testing purposes, a numeric value indicating the max number of trees to be cached.
#' @return A list of 4 elements:
#' \describe{
#'   \item{authors}{A list of lists of author names of the original studies that
#'     published chronograms currently stored in the Open Tree of Life database.}
#'   \item{curators}{A list of lists of curator names that uplodaded chronograms
#'     to the Open Tree of Life database.}
#'   \item{studies}{A list of study identifiers from original studies that
#'     published chronograms currently stored in the Open Tree of Life database.}
#'   \item{trees}{A `multiPhylo` object storing the chronograms from Open Tree of
#'     Life database.}
#' }
#' @export
get_otol_chronograms <- function(max_tree_count = "all") {
  options(warn = 1)
  start_time <- Sys.time() # to register run time
  chronogram_matches <- rotl::studies_find_trees(property = "ot:branchLengthMode",
                                                 value = "ot:time",
                                                 detailed = TRUE,
                                                 verbose = FALSE)
  # rotl::studies_find_trees throws this Warning in rbind(c(n_trees = "1", tree_ids = "tree6160", candidate = "tree6160",  :
  # number of columns of result is not a multiple of vector length (arg 1)  colnames(chronogram_matches)
  # rotl::list_trees(chronogram_matches)
  trees <- list()
  authors <- list()
  curators <- list()
  studies <- list()
  dois <- list()
  tree_count <- 0
  bad_ones <- c()
  ott_id_problems <- data.frame(study_id = character(), tree_id = character()) # nrow(fix_negative_brlen) is 0
  # utils::opentree_chronograms
  for (study_index in sequence(dim(chronogram_matches)[1])) {
    # the only purpose for the following conditional is testing:
    if (is.numeric(max_tree_count)) {
      if (tree_count > max_tree_count) {
        break
      }
    }
    message("Downloading tree(s) from study ", study_index, " of ", dim(chronogram_matches)[1])

    for (chrono_index in sequence((chronogram_matches$n_matched_trees[study_index]))) {
      study_id <- chronogram_matches$study_ids[study_index]
      # study_id <- "ot_409" # largest tree in Hedges et al. 2015
      # chronogram_matches[chronogram_matches$study_ids == study_id,]
      # tree_id <- "tree2"
      # if(study_id %in% opentree_chronograms$studies)  # unlist(opentree_chronograms$studies) if new_studies_only
      # 	new_tree <- get_study_tree(study_id=study_id, tree_id=tree_id, tip_label='ott_taxon_name')
      try_tree <- new_tree2 <- new_tree <- NULL
      tree_id <- strsplit(chronogram_matches$match_tree_ids[study_index], ", ")[[1]][chrono_index]
      if (is.null(tree_id)) {
        ott_id_problems <- rbind(ott_id_problems, data.frame(study_id, tree_id))
      }
      potential_bad <- paste("tree_id='", tree_id, "', study_id='", study_id, "'", sep = "")

      if (!grepl("\\.\\.\\.", tree_id) & !is.na(tree_id)) { # to deal with ellipsis bug
        message("tree_id = '", tree_id, "', study_id = '", study_id, "'")
        new_tree <- tryCatch(rotl::get_study_tree(study_id = study_id,
                                                  tree_id = tree_id,
                                                  tip_label = "ott_taxon_name",
                                                  deduplicate = TRUE),
                             error = function(e) NULL)
        # would like to dedup; don't use get_study_subtree, as right now it doesn't take tip_label args
        # try(new_tree <- datelife:::get_study_tree_with_dups(study_id=study_id,tree_id=tree_id ))
        # try(new_tree <- rotl::get_study_subtree(study_id=study_id,tree_id=tree_id, tip_label="ott_taxon_name", subtree_id="ingroup")) #only want ingroup, as it's the one that's been lovingly curated.
        if (!is.null(new_tree) & phylo_has_brlen(phy = new_tree)) {
          # add ott_ids
          # right now the function is having trouble to retrieve trees with ott ids as tip labels from certain studies
          new_tree2 <- tryCatch(rotl::get_study_tree(study_id = study_id,
                                                     tree_id = tree_id,
                                                     tip_label = "ott_id"),
                                error = function(e) NULL)
          utils::head(new_tree2$tip.label)
          if (!inherits(new_tree2, "phylo")) {
            problem <- "otol database does not have ott ids for this tree (names have not been curated)"
            ott_id_problems <- rbind(ott_id_problems, data.frame(study_id, tree_id))
            # do tnrs here?
            # new.tree_tnrs <- tnrs_match(input = new_tree$tip.label)
            # new.tree_tnrs <- tnrs_match.phylo(input = new_tree)
          }
          new_tree$ott_ids <- gsub("_.*", "", new_tree2$tip.label) # if new_tree2 is null it will generate an empty vector
          try_tree <- clean_ott_chronogram(new_tree)
          # previous line will give NA if just one or no tip labels are mapped to ott???
          if (phylo_has_brlen(phy = try_tree)) {
            new_tree <- try_tree
            if (is_good_chronogram(new_tree)) {
              new_tree$tip.label <- gsub("_", " ", new_tree$tip.label)
              message("\t", "has tree with branch lengths")
              doi <- NULL
              try(doi <- gsub("https?://(dx\\.)?doi.org/", "", attr(rotl::get_publication(rotl::get_study_meta(study_id)), "DOI")))
              authors <- append(authors, NA)
              if (length(doi) == 0) {
                warning(paste(study_id, "has no DOI attribute, author names will not be retrieved."))
              } else {
                try(authors[length(authors)] <- list(paste(as.character(knitcitations::bib_metadata(doi)$author))))
              }
              curators <- append(curators, NA)
              try(curators[length(curators)] <- list(rotl::get_study_meta(study_id)[["nexml"]][["^ot:curatorName"]]))
              try(studies <- append(studies, study_id))
              tree_count <- tree_count + 1
              # print(tree_count)
              try(dois <- append(dois, chronogram_matches$study_doi[study_index]))
              trees[[tree_count]] <- new_tree
              names(trees)[tree_count] <- rotl::get_publication(rotl::get_study_meta(study_id))[1]
              message("\t", "was good tree")
              potential_bad <- NULL
            }
          }
          # add taxonomic nodelabels to trees object here
          new_tree <- map_nodes_ott(tree = new_tree)
        }
      } else {
        warning("Not all trees could be loaded from this study due to ellipsis bug, https://github.com/ropensci/rotl/issues/85")
      }
      if (!is.null(potential_bad)) {
        bad_ones <- append(bad_ones, potential_bad)
      }
    }
  }
  message("Problematic combos:")
  message(paste0(utils::capture.output(bad_ones), collapse = "\n"))
  for (i in sequence(length(authors))) {
    if (!any(is.na(authors[[i]])) & length(authors[[i]]) > 0) {
      Encoding(authors[[i]]) <- "latin1"
      authors[[i]] <- iconv(authors[[i]], "latin1", "UTF-8")
    }
    for (j in sequence(length(curators[[i]]))) {
      if (!any(is.na(curators[[i]][[j]])) & length(curators[[i]][[j]]) > 0) {
        Encoding(curators[[i]][[j]]) <- "latin1"
        curators[[i]][[j]] <- iconv(curators[[i]][[j]], "latin1", "UTF-8")
      }
    }
  }
  if (nrow(ott_id_problems) > 0) {
    problems_file <- file.path(tempdir(), paste0("ott_id_problems_", max_tree_count, ".csv"))
    utils::write.csv(ott_id_problems,
      file = problems_file,
      quote = FALSE, row.names = FALSE
    )
    message("Problematic chronograms were saved to ", problems_file)
  } else {
    message("There were no problematic chronograms.")
  }
  tot_time <- Sys.time() - start_time # end of registering function running time
  class(trees) <- "multiPhylo"
  result <- list(trees = trees, authors = authors, curators = curators, studies = studies, dois = dois)
  attr(result, "running_time") <- tot_time
  message(tot_time)
  return(result)
}


#' Check if a tree is a valid chronogram.
#' @inheritParams phylo_check
#' @return `TRUE` if it is a valid tree.
#' @export
is_good_chronogram <- function(phy) {
  passing <- TRUE
  if (!inherits(phy, "phylo")) {
    passing <- FALSE
    warning("tree failed over not being class phylo")
  }
  if (ape::Ntip(phy) <= ape::Nnode(phy)) {
    passing <- FALSE
    warning("tree failed over not having more internal nodes than tips")
  }
  if (length(which(grepl("not.mapped", phy$tip.label))) > 0) {
    warning("tree failed over having not mapped taxa that should have been checked")
    passing <- FALSE # not cleaned properly
  }
  # enhance: test that there are no duplicated labels in chronogram:
  if (any(is.na(phy$tip.label))) {
    passing <- FALSE
    warning("tree failed over having NA for tips")
  } else {
    # remove this if we want to preserve all original labels:
    if (min(nchar(phy$tip.label)) <= 2) {
      passing <- FALSE
      warning("tree failed for having names of two or fewer characters") # why is this important?
    }
  }
  if (!ape::is.rooted(phy)) {
    passing <- FALSE
    warning("tree failed over not being rooted")
  }
  if (!ape::is.ultrametric(phy, tol = 0.01, option = 2)) {
    passing <- FALSE
    warning("tree failed over not being ultrametric (NOTE: this condition should be removed for paleo trees)")
  }
  return(passing)
}


#' Problematic chronograms from Open Tree of Life.
#'
#' @name problems
#' @docType data
#' @format A list of trees with unmapped taxa
#' @source \url{http://opentreeoflife.org}
#' @keywords otol tree chronogram unmapped tnrs
#' @details
#' Before we developed tools to clean and map tip labels for our cached trees
#' we found some trees that were stored with unmapped tip labels
#' we extracted them and saved them to be used for testing functions.
#' Generated with
#' problems <- opentree_chronograms$trees[sapply(sapply(opentree_chronograms$trees, "[", "tip.label"), function(x) any(grepl("not.mapped", x)))]
#' usethis::use_data(problems)
#' opentree_chronograms object from commit https://github.com/phylotastic/datelife/tree/be894448f6fc437241cd0916fab45e84ac3e09c6
"problems"

#' Clean up some issues with Open Tree of Life chronograms
#' For now it 1) checks unmapped taxa and maps them with tnrs_match.phylo, 2) roots the chronogram if unrooted
#'
#' @inheritParams phylo_check
# # ' @return A cleaned up phylo object
#' @inherit tnrs_match.phylo return details
#' @export
clean_ott_chronogram <- function(phy) {
  phylo_check(phy, dated = FALSE)
  phy.ori <- phy
  # homogenize tip labels to blanks (no underscores)
  # it will speed up tnrs_match_names considerably bc no approximate matching needed
  phy <- phylo_tiplabel_underscore_to_space(phy)
  unmapped.taxa <- unique(c(which(nchar(phy$tip.label) <= 2), which(grepl("not.mapped", phy$tip.label)))) # numeric of indices
  phy$tip.label[unmapped.taxa] <- sub(".*-.", "", phy$tip.label[unmapped.taxa]) # this gets the original label and gets rid of the not.mapped tag
  phy$tip.label[unmapped.taxa] <- gsub("aff ", "", phy$tip.label[unmapped.taxa]) # removes aff tag
  phy$tip.label[unmapped.taxa][stringr::str_count(phy$tip.label[unmapped.taxa], " ") >= 2] <- gsub("^([^ ]* [^ ]*) .*$", "\\1", phy$tip.label[unmapped.taxa][stringr::str_count(phy$tip.label[unmapped.taxa], " ") >= 2])
  # grepl("^([^_]*_[^_]*)_.*$", tip.labels[stringr::str_count(tip.labels, "_")>=2])  # this works with underscores
  # drop duplicated tips now:
  tipstodrop <- c()
  # identify not mapped taxa that have mapped duplicates first
  # in this case, drop the unmapped tips and keep the mapped ones.
  # identify unmapped tips with duplicates among mapped tips (meaning that it does not have to be remapped):
  cond <- match(unique(phy$tip.label[unmapped.taxa]), phy$tip.label[-unmapped.taxa])
  if (any(!is.na(cond))) {
    # the following only gets the first match:
    # mm <- match(unique(phy$tip.label[unmapped.taxa])[!is.na(cond)], phy$tip.label[unmapped.taxa])
    # this gets all unmapped.taxa that match:
    mm <- match(phy$tip.label[unmapped.taxa], unique(phy$tip.label[unmapped.taxa])[!is.na(cond)])
    # tipstodrop <- c(tipstodrop, unmapped.taxa[!is.na(mm)]) # no need for this now
    unmapped.taxa <- unmapped.taxa[is.na(mm)] # update unmapped.taxa object (removing taxa in mapped tips that are duplicated in unmapped.taxa, preventing unnecesary calls for tnrs_match_names)
  }
  # identify duplicated taxa within unmapped tips (necessary to drop or modify tips afterwards)
  # dd <- duplicated(phy$tip.label[unmapped.taxa])
  # tipstodrop <- c(tipstodrop, unmapped.taxa[dd])
  # this ones should be mapped anyways (tnrs_match.phylo won't spend more time if we leave duplicates to be mapped)
  # and we actually need to check the duplicates after running tnrs_match.phylo
  phy$mapped <- rep("ott", length(phy$tip.label))
  if (length(unmapped.taxa) > 0) {
    phy <- tnrs_match.phylo(phy, tip = unmapped.taxa) # this also creates an original tip labels element
  }
  dd <- duplicated(phy$tip.label)
  while (any(dd)) {
    # phy <- ape::drop.tip(phy, tip = unique(tipstodrop)) # do not drop duplicated tips, just modify them
    phy$tip.label[dd] <- paste0(phy$tip.label[dd], "_dup") # enhance: do we actually need to modify dups???
    dd <- duplicated(phy$tip.label) # checks if there are still any dups
  }
  if (!ape::is.rooted(phy) & ape::is.ultrametric(phy, option = 2)) {
    phy$root.edge <- 0
  }
  return(phy)
}

# #' Sample of raw open tree chronograms for testing
# #'
# #' @name raw_opentree_chronograms
# #' @docType data
# #' @format A multiPhylo object
# #' @source \url{http://opentreeoflife.org}
# #' @keywords otol chronogram test
# #' @details
# #'
# #' Generated with
# #' opentree_chronograms_ntip <- unname(sapply(opentree_chronograms$trees, ape::Ntip))
# #' study.ids <- unique(unlist(opentree_chronograms$studies[which(opentree_chronograms_ntip < 100 & opentree_chronograms_ntip > 50)]))
# #'
# #' study.ids <- c("ot_122")
# #' tree.ids <- c("tree1")
# #' raw_opentree_chronograms <- lapply(seq(length(study.ids)), function(x){
# #' 	rotl::get_study_tree(study_id=study.ids[x],tree_id=tree.ids[x], tip_label="ott_taxon_name")
# #' })
# "raw_opentree_chronograms"


#' Add Open Tree of Life Taxonomy to tree nodes.
#' @inheritParams tree_fix_brlen
#' @return A `phylo` object with "nodelabels".
#' @examples
#' \dontrun{ # This is a flag for package development. You are welcome to run the example.
#'
#' # Load the Open Tree chronograms database cached in datelife:
#' utils::data(opentree_chronograms)
#'
#' # Get the small chronograms (i.e., chronograms with less that ten tips) to generate a pretty plot:
#' small <- opentree_chronograms$trees[unlist(sapply(opentree_chronograms$trees, ape::Ntip)) < 10]
#'
#' # Now, map the Open Tree taxonomy to the nodes of the first tree
#' phy <- map_nodes_ott(tree = small[[1]])
#' # and plot it:
#' # plot_phylo_all(phy)
#' library(ape)
#' plot(phy)
#' nodelabels(phy$node.label)
#'
#' } #end dontrun
#' @export
map_nodes_ott <- function(tree) {
  # utils::data(opentree_chronograms)
  # ss <- opentree_chronograms$trees[unlist(sapply(opentree_chronograms$trees, ape::Ntip)) < 10]
  # tree <- ss[[1]]
  # phy <- opentree_chronograms$trees[[i]]
  phy <- tree_check(tree = tree, dated = FALSE)
  cc <- tryCatch(classification_paths_from_taxonomy(phy$tip.label, sources = "Open Tree of Life Reference Taxonomy"),
    error = function(e) list(resolved = c())
  ) # this traps the error when phy = opentree_chronograms$trees[[50]] from load(data-raw/opentree_chronograms_oct2018.rda)
  if (length(cc$resolved) == 0) { # when no tip label could be mapped to Opentree Taxonomy
    message("Tip labels could not be mapped to Open Tree of Life Taxonomy")
    return(phy)
  }
  paths <- cc$resolved$classification_path
  ranks <- cc$resolved$classification_path_ranks
  # paths <- gsub('Not assigned', "NA", paths) # CoL gives Not assigned for some taxa. We need these in here. Then we split
  paths <- strsplit(paths, "\\|")
  ranks <- strsplit(ranks, "\\|")
  # sapply(strsplit(aa$classification_path, "\\|"), matrix, nrow = 1, byrow = T)
  bb <- lapply(paths, matrix, nrow = 1, byrow = T)
  bb <- lapply(seq(length(paths)), function(x) {
    colnames(bb[[x]]) <- ranks[[x]]
    bb[[x]]
    # as.data.frame(bb[[x]])
    # bb[[x]]
  })
  # dt <- as.matrix(data.table::rbindlist(as.data.frame(bb), fill = TRUE))
  dt <- plyr::rbind.fill.matrix(bb) # works better than previous line
  rownames(dt) <- cc$resolved$user_supplied_name
  # dt <- dt[,-grep("species", colnames(dt))]  # not needed
  # mm <- matrix(unlist(strsplit(aa$classification_path, "\\|")), nrow = length(aa$classification_path), byrow=T)
  # enhance: nodelabel.phylo only keeps one taxon name per node, and more than often there are several names per nodes
  # find a way to keep them all
  gg <- suppressWarnings(geiger::nodelabel.phylo(phy = phy, taxonomy = dt, ncores = 1))
  # got <- get_ott_lineage(phy)
  return(gg)
}
