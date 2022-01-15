# functions (from phunding package) to deal with open tree Taxonomy rotl
# functions, namely, rotl::taxonomy_taxon_info()

#' Identifies, extracts and cleans children names from a [taxonomy_taxon_info()]
#' output
#'
#' @param taxon_info An output of [rotl::taxonomy_taxon_info()].
#' @param invalid A character vector of flags used by Open Tree of Life Taxonomy
#'   to detect invalid taxon names.
#' @return A list with valid children unique OTT names, ott_ids and ranks.
#' @details Eliminates all taxa that will give problems when trying to retrieve
#'   an induced subtree from Open Tree of Life.
#' @examples
#' # Some suppressed taxa within "Gleicheniales" family of ferns:
#' tt <- rotl::taxonomy_taxon_info(3866048)
#' rotl::taxonomy_taxon_info(5262766)
#' @export
clean_taxon_info_children <- function(taxon_info, invalid = c("barren", "extinct", "uncultured", "major_rank_conflict", "incertae_sedis", "unplaced", "conflict", "environmental", "not_otu", "hidden", "hybrid")) {
  # invalid <- c("BARREN", "EXTINCT", "UNCULTURED", "MAJOR_RANK_CONFLICT", "INCERTAE", "UNPLACED", "CONFLICT")
  # names(taxon_info[[2]])
  for (i in seq(taxon_info)) {
    # length(tax_info[[i]][[1]]$children)
    # sapply(tax_info[[i]][[1]]$children, length)
    # length(sapply(sapply(tax_info[[i]][[1]]$children, "[", "flags"), function(y) unlist(tolower(y))))
    # sapply(sapply(tax_info[[i]][[1]]$children, "[", "flags"), function(y) unlist(tolower(y)))[1]
    ii <- lapply(sapply(sapply(taxon_info[[i]]$children, "[", "flags"), unlist), function(y) any(invalid %in% y))
    # ii <- lapply(sapply(sapply(taxon_info[[i]][[1]]$children, "[", "flags"), unlist), function(y) any(toupper(invalid) %in% y))
    if (length(ii) > 0) {
      taxon_info[[i]]$children <- taxon_info[[i]]$children[!unlist(ii)]
    }
    # now clean names with sp. or environmental
    ii <- unlist(sapply(c("sp\\.", "environmental", "unclassified", "incertae"), function(x) grep(x, sapply(sapply(taxon_info[[i]]$children, "[", "unique_name"), unlist))))

    if (length(ii) > 0) {
      taxon_info[[i]]$children <- taxon_info[[i]]$children[-unique(ii)]
    }
    # taxon_info[[i]][[1]]$children <- taxon_info[[i]][[1]]$children[!unlist(ii)]
  }
  return(taxon_info)
}

#' Checks input for [get_ott_clade()], [get_ott_children()], and
#' [get_otol_synthetic_tree()] functions
#'
#' @param input Optional. A character vector of names or a `datelifeQuery` object
#' @param ott_ids If not NULL, it takes this argument and ignores input. A
#'   numeric vector of ott ids obtained with [rotl::taxonomy_taxon_info()] or
#'   [rotl::tnrs_match_names()] or [tnrs_match()]
#' @inheritDotParams make_datelife_query
#' @return A named numeric vector of valid Open Tree Taxonomy (OTT) ids.
#' @export
#' @details By default, it uses the `ott_id` argument if it is not NULL.
check_ott_input <- function(input = NULL, ott_ids = NULL, ...) {
  # for debugging:
  # utils::data(threebirds.rda)
  # input <- threebirds_query$cleaned_names
  # ott_ids <- threebirds_query$ott_ids

  if (is.null(input) & is.null(ott_ids)) {
    message("Both input arguments are empty; please specify one.")
    return(NA)
  }
  if (is.null(ott_ids) | all(is.na(ott_ids))) {
    # checks that input is a datelifeQuery object, otherwise it uses make_datelife_query on input
    if (!is_datelife_query(input)) {
      input <- make_datelife_query(input, ...)
    }
    if (is.numeric(input$ott_id)) {
      ott_ids <- input$ott_ids
      names(ott_ids) <- input$cleaned_names
    } else {
      input_tnrs <- datelife::tnrs_match(input = input$cleaned_names)
      # should we clean tnrs from invalid? what about NA's? Maybe only NA's at the end, so we can tell users what was unsuccessful
      # if we decide to clean, the two following lines should be uncommented:
      # input_tnrs <- clean_tnrs(tnrs = input_tnrs)
      # input_tnrs <- input_tnrs[!is.na(input_tnrs$unique_name),]  # gets rid of names not matched with rotl::tnrs_match_names; otherwise rotl::tol_induced_subtree won't run
      ott_ids <- input_tnrs$ott_id
      names(ott_ids) <- input_tnrs$unique_name
    }
  }
  # ott_ids <- suppressWarnings(as.numeric(ott_ids))
  if (!is.numeric(ott_ids)) {
    ott_ids2 <- ott_ids
    message("ott_ids is not a numeric vector; coercing to numeric.")
    ott_ids <- as.numeric(ott_ids)
    if (any(is.na(is.numeric(ott_ids2)))) {
      message(paste0("\nInput '", paste(names(input)[which(is.na(ott_ids2))], collapse = "', '"), "' not numeric."))
    }
  }
  if (is.null(names(ott_ids))) {
    # ott_ids <- c(1, ott_ids) # if ott_ids does not exist it will give an error
    # so we're catching it in an sapply, and giving "error" when it errors.
    names(ott_ids) <- unlist(sapply(seq(length(ott_ids)), function(i) {
      tryCatch(rotl::tax_name(rotl::taxonomy_taxon_info(ott_ids = ott_ids[i])),
        error = function(e) "error"
      )
    }))
    # then we subset it:
    if (any(grepl("error", names(ott_ids)))) {
      message(paste0("\nott_id '", paste(ott_ids[grepl("error", names(ott_ids))], collapse = "', '"), "' not found in Open Tree of Life Taxonomy."))
      ott_ids <- ott_ids[!grepl("error", names(ott_ids))]
    }
  }
  # now check for NAs:
  if (any(is.na(ott_ids))) {
    message(paste0("\nInput '", paste(names(input)[which(is.na(ott_ids))], collapse = "', '"), "' not found in Open Tree of Life Taxonomy."))
    ott_ids <- ott_ids[which(!is.na(ott_ids))]
  }
  if (length(ott_ids) < 1) {
    message("At least one valid input name or numeric ott IDs are needed to get any information")
    return(NA)
  }
  return(ott_ids)
}
#' Get the ott id and name of all lineages from one or more input taxa
#' @inheritParams check_ott_input
#' @return A list of named numeric vectors of ott ids from input and all the clades it belongs to.
#' @examples
#' \dontrun{ # This is a flag for package development. You are welcome to run the example.
#' taxa <- c("Homo", "Bacillus anthracis", "Apis", "Salvia")
#' lin <- get_ott_lineage(taxa)
#' lin
#' # Look up an unknown OTT id:
#' get_ott_lineage(ott_id = 454749)
#' } # end dontrun
#' @export
get_ott_lineage <- function(input = NULL, ott_ids = NULL) {
  # ott_ids <- c(335590, 555178, 748370, 1070795, 3942422, 907458, 472526, 820645, 31926, 756728, 605194, 490099)
  # ott_ids <- c(417116, 120960)
  # ott_ids <- 417116
  input_ott_match <- check_ott_input(input, ott_ids)
  tax_info <- .get_ott_lineage(input_ott_match)
  lin <- lapply(tax_info, "[", "lineage")
  ott_names <- lapply(lin, function(x) unlist(sapply(x[[1]], "[", "unique_name")))
  # unlist(sapply(lin[[1]][[1]], "[", "unique_name"))
  # length(ott_names) == length(tax_info)
  # ott_ids <- sapply(seq(length(lin)), function(x) stats::setNames(unlist(sapply(lin[[x]][[1]], "[", "ott_id")), ott_names[[x]]))
  ott_ids <- lapply(seq(length(lin)), function(x) unlist(sapply(lin[[x]][[1]], "[", "ott_id")))
  ott_ranks <- lapply(seq(length(lin)), function(x) unlist(sapply(lin[[x]][[1]], "[", "rank")))
  # if(length(ott_ids) == 1){
  #     res <- list(matrix(c(ott_ids, ott_ranks), ncol =2, dimnames = list(ott_names, c("ott_ids", "ott_ranks"))))
  # } else {
  mat <- function(x) {
    matrix(c(ott_ids[[x]], ott_ranks[[x]]), ncol = 2, dimnames = list(ott_names[[x]], c("ott_ids", "ott_ranks")))
  }
  res <- lapply(seq(length(lin)), mat)
  # }
  # enhance: ott_id names should be the name of the rank, look at the example to see why
  names(res) <- names(input_ott_match)
  # stats::setNames(ott_ids, names(input_ott_match))
  return(res)
}
#' Get the lineage of a set of taxa using rotl:taxonomy_taxon_info(include_lineage = TRUE)
#' @param input_ott_match An Output of check_ott_input function.
#' @return A taxonomy_taxon_info object
.get_ott_lineage <- function(input_ott_match) {
  tax_info <- vector(mode = "list", length = length(input_ott_match))
  progression <- utils::txtProgressBar(min = 0, max = length(tax_info), style = 3)
  for (i in seq(length(input_ott_match))) {
    tax_info[i] <- tryCatch(rotl::taxonomy_taxon_info(input_ott_match[i], include_lineage = TRUE),
      error = function(e) NA
    )
    utils::setTxtProgressBar(progression, i)
  }
  cat("\n") # just to make the progress bar look better
  tax_info
}


#' Get the ott id and name of one or several given taxonomic ranks from one or
#' more input taxa
#' @inheritParams check_ott_input
#' @inheritParams get_ott_children
#' @return A list of named numeric vectors with ott ids from input and all requested ranks.
#' @export
get_ott_clade <- function(input = NULL, ott_ids = NULL, ott_rank = "family") {
  # ott_ids= c('649007', '782239', '782231', '1053057', '372826', '766272', '36015', '914245', '873016', '684051')
  # ott_ids = c('431493', '431493', '431493', '431493', '431493', '431493', '431493', '429482', '429482', '429482')
  rank <- ott_rank
  input_ott_match <- check_ott_input(input, ott_ids)
  tax_info <- .get_ott_lineage(input_ott_match)
  # names(tax_info[[10]])
  # sapply(tax_info[[10]]$lineage, "[", "rank")
  # length(tax_info[[10]]$lineage)
  # names(tax_info[[10]]$lineage[[1]])
  # tax_info[[10]]$lineage[[1]]$flags
  input_ott_names <- unlist(sapply(tax_info, "[", "unique_name"))
  rank_ott_ids <- rank_names <- lapply(seq(length(rank)), function(x) rep(NA, length(input_ott_names)))
  # I still need to drop all invalid lineages first here!!!
  for (i in seq(length(rank))) {
    ready <- grepl(rank[i], unlist(sapply(tax_info, "[", "rank")))
    if (any(ready)) {
      rank_names[[i]][ready] <- unlist(sapply(tax_info[ready], "[", "unique_name"))
      rank_ott_ids[[i]][ready] <- unlist(sapply(tax_info[ready], "[", "ott_id"))
    }
    if (!all(ready)) {
      rank_names[[i]][!ready] <- sapply(tax_info[!ready], function(x) {
        lin <- tryCatch(x$lineage[grep(paste0("^", rank[i], "$"), sapply(x$lineage, "[", "rank"))][[1]]$unique_name,
          error = function(e) NA
        )
        return(lin)
      })
      rank_ott_ids[[i]][!ready] <- sapply(tax_info[!ready], function(x) {
        lin <- tryCatch(x$lineage[grep(paste0("^", rank[i], "$"), sapply(x$lineage, "[", "rank"))][[1]]$ott_id,
          error = function(e) NA
        )
        return(lin)
      })
    }
    names(rank_ott_ids[[i]]) <- rank_names[[i]]
    # length(rank_ott_ids[i])
    # stop()
  }
  names(input_ott_match) <- input_ott_names
  res <- c(list(input_ott_match), rank_ott_ids)
  names(res) <- c("input", rank)
  return(res)
}

# .get_ott_clade <- function(input, rank){
#     tax_info <- tryCatch(rotl::taxonomy_taxon_info(input, include_lineage = TRUE),
#     error = function(e) NA)
#     input_ott_names <- tax_info$unique_name
#     rank_ott_ids <- rank_names <- vector(mode = "list", length = length(rank))
#     # I still need to drop all invalid lineages first here!!!
#     for (i in seq(length(rank))){
#         if(rank[i] %in% tax_info[[1]]$rank){
#             rank_names <- tax_info[[1]]$unique_name
#             rank_ott_ids <- tax_info[[1]]$ott_id
#         } else {
#             rank_names[[i]] <- sapply(tax_info, function(x) {
#                 lin <- tryCatch(x$lineage[grep(paste0("^", rank[i], "$"), sapply(x$lineage, "[", "rank"))][[1]]$unique_name,
#                 error = function(e) NA)
#                 return(lin)
#             })
#             rank_ott_ids[[i]] <- sapply(tax_info, function(x) {
#                 lin <- tryCatch(x$lineage[grep(paste0("^", rank[i], "$"), sapply(x$lineage, "[", "rank"))][[1]]$ott_id,
#                 error = function(e) NA)
#                 return(lin)
#             })
#         }
#         names(rank_ott_ids[[i]]) <- rank_names[[i]]
#         # length(rank_ott_ids[i])
#         # stop()
#     }
#     return(rank_ott_ids)
# }

#'
#' Extract valid children from given taxonomic name(s) or OpenTree Taxonomic ids (not from a taxonomy
#' taxon info object).
#' @inheritParams check_ott_input
#' @param taxonomic_source A character vector with the desired taxonomic sources. Options are "ncbi", "gbif" or "irmng". Any other value will retrieve data from all taxonomic sources. The function defaults to "ncbi".
#' @return A named list containing valid taxonomic children of given taxonomic name(s).
#' @details
#' GBIF and other taxonomies contain deprecated taxa that are not marked as such in the Open Tree Taxonomy.
#' We are relying mainly in the NCBI taxonomy for now.
#' @examples
#' # genus Dictyophyllidites with ott id = 6003921 has only extinct children
#' # in cases like this the same name will be returned
#' tti <- rotl::taxonomy_taxon_info(6003921, include_children = TRUE)
#' gvc <- get_valid_children(ott_ids = 6003921)
#' # More examples:
#' get_valid_children(ott_ids = 769681) # Psilotopsida
#' get_valid_children(ott_ids = 56601) # Marchantiophyta
# input= "Erythrospiza" # obsolete genus that is still on ott
# input = c("Felis", "Homo", "Malvaceae")
# input = "Telespiza"
#' @export
get_valid_children <- function(input = NULL, ott_ids = NULL, taxonomic_source = "ncbi") {
  taxonomic_source_here <- tryCatch(match.arg(taxonomic_source, c("ncbi", "gbif", "irmng")),
    error = function(e) NULL
  )
  input_ott_match <- check_ott_input(input, ott_ids)
  all_children <- vector(mode = "list", length = length(input_ott_match))
  # monotypic <- vector(mode = "logical", length = length(input_ott_match))
  progression <- utils::txtProgressBar(min = 0, max = length(all_children), style = 3)
  for (i in seq(length(input_ott_match))) {
    tt <- tryCatch(rotl::taxonomy_taxon_info(input_ott_match[i], include_children = TRUE),
      error = function(e) NA
    )
    # length(tt[[1]])
    tt <- clean_taxon_info_children(tt) # removes all invalid children
    if (length(tt[[1]]$children) > 0) {
      # sapply(tt[[1]]$children, "[", "flags")
      # sapply(tt[[1]]$children, "[", "unique_name")
      # which(unlist(sapply(tt[[1]]$children, "[", "unique_name")) == "Mesangiospermae")
      # tt[[1]]$children[108]
      rr <- unname(unlist(sapply(tt[[1]]$children, "[", "rank")))
      child <- unname(unlist(sapply(tt[[1]]$children, "[", "ott_id")))
      tax_sources <- unname(lapply(seq(length(tt[[1]]$children)), function(k) {
        unlist(tt[[1]]$children[[k]]$tax_sources)
      }))
      # names(tt[[1]]$children[[1]])
      # unname(unlist(sapply(tt[[1]]$children, "[", "flags")))
      # unname(unlist(sapply(tt[[1]]$children, "[", "is_suppressed")))
      # unname(unlist(sapply(tt[[1]]$children, "[", "source")))
      names(tax_sources) <- names(child) <- names(rr) <- unname(unlist(sapply(tt[[1]]$children, "[", "unique_name")))
      monotypic <- FALSE
    } else {
      child <- tt[[1]]$ott_id
      rr <- tt[[1]]$rank
      # tax_sources <- list(unlist(tt[[1]]$children[[2]]$tax_sources))
      tax_sources <- list(unlist(tt[[1]]$tax_sources))
      names(tax_sources) <- names(child) <- names(rr) <- tt[[1]]$unique_name
      monotypic <- TRUE
    }
    if (inherits(taxonomic_source_here, "character")) {
      keep <- sapply(tax_sources, function(x) any(grepl(taxonomic_source_here, x)))
      child <- child[keep]
      rr <- rr[keep]
    }
    all_children[[i]] <- list(children = data.frame(ott_id = child, rank = rr), is_monotypic = monotypic)
    utils::setTxtProgressBar(progression, i)
  }
  cat("\n") # just to make the progress bar look better
  names(all_children) <- names(input_ott_match)
  return(all_children)
}

#' Use this instead of [rotl::tol_subtree()] when taxa are not in synthesis tree and
#' you still need to get all species or an induced OpenTree subtree
#' @inheritParams check_ott_input
#' @param ott_rank A character vector with the ranks you wanna get lineage children from.
#' @param ... Other arguments to pass to [get_valid_children()].
#' @return A `data.frame` object.
#' @examples
#' # An example with the dog genus:
#'
#' # It is currently not possible to get an OpenTree subtree of a taxon that is
#' #  missing from the OpenTree synthetic tree.
#' # The dog genus is not monophyletic in the OpenTree synthetic tree, so in
#' #  practice, it has no node to extract a subtree from.
#' tnrs <- tnrs_match("Canis")
#' \dontrun{ # This is a flag for package development. You are welcome to run the example.
#' rotl::tol_subtree(tnrs$ott_id[1])
#' #> Error: HTTP failure: 400
#' #> [/v3/tree_of_life/subtree] Error: node_id was not found (broken taxon).
#' } # end dontrun
#' ids <- tnrs$ott_id[1]
#' names(ids) <- tnrs$unique_name
#' children <- get_ott_children(ott_ids = ids) # or
#' children <- get_ott_children(input = "Canis")
#' str(children)
#' ids <- children$Canis$ott_id
#' names(ids) <- rownames(children$Canis)
#' tree_children <- datelife::get_otol_synthetic_tree(ott_ids = ids)
#' plot(tree_children, cex = 0.3)
#'
#' # An example with flowering plants:
#'
#' \dontrun{ # This is a flag for package development. You are welcome to run the example.
#' oo <- get_ott_children(input = "magnoliophyta", ott_rank = "order")
#' # Get the number of orders of flowering plants that we have
#' sum(oo$Magnoliophyta$rank == "order")
#' } # end dontrun
#' @export
# children <- get_ott_children(input = "Fringillidae")
# rotl::tol_subtree(ott_id = 839319) #broken taxa
# children <- get_ott_children(input = "Cetaceae")
# nrow(children[[1]])
# sum(children[[1]]$rank == "species")
get_ott_children <- function(input = NULL, ott_ids = NULL, ott_rank = "species", ...) {
  # iput <- c("Felis", "Homo", "Malvaceae")
  input_ott_match <- check_ott_input(input, ott_ids)
  # names(input_ott_match)
  all_children <- vector(mode = "list", length = length(input_ott_match))
  input_ott_match_ranks <- unlist(sapply(rotl::taxonomy_taxon_info(input_ott_match), "[", "rank"))
  # grepl(paste0("\\b", ott_rank, "\\b"), c("species", "subspecies"))
  ott_ranki <- grepl(paste0("\\b", ott_rank, "\\b"), input_ott_match_ranks)
  df <- data.frame(ott_id = input_ott_match[ott_ranki], rank = rep(ott_rank, sum(ott_ranki)))
  all_children[ott_ranki] <- lapply(seq(nrow(df)), function(i) df[i, ])
  # all_children[c(T,F,T)] <- input[c(T,F,T)]
  # sum(c(T,F,T))
  # progression <- utils::txtProgressBar(min = 0, max = length(all_children), style = 3)
  if (any(sapply(all_children, is.null))) {
    input_ott_match <- input_ott_match[sapply(all_children, is.null)]
    for (i in seq(length(input_ott_match))) {
      mm <- data.frame(ott_ids = vector(mode = "numeric", length = 0), rank = vector(mode = "logical", length = 0))
      vv <- get_valid_children(ott_ids = input_ott_match[i], ...)
      success <- vv[[1]]$children$rank == ott_rank | vv[[1]]$is_monotypic
      if (any(success)) {
        mm <- rbind(mm, vv[[1]]$children[success, ])
      }
      while (!all(success)) {
        vv <- get_valid_children(ott_ids = unlist(sapply(sapply(
          vv, "[",
          "children"
        ), "[", "ott_id"))[!success], ...)
        if (any(unlist(sapply(vv, "[", "is_monotypic")))) {
          dd <- do.call("rbind", sapply(vv[unlist(sapply(vv, "[", "is_monotypic"))], "[", "children"))
          rownames(dd) <- unname(unlist(sapply(sapply(vv[unlist(sapply(vv, "[", "is_monotypic"))], "[", "children"), rownames)))
          # rownames(dd) <- gsub("\\..*", "", rownames(dd))
          mm <- rbind(mm, dd)
          vv <- vv[!unlist(sapply(vv, "[", "is_monotypic"))]
        }
        success <- unlist(sapply(sapply(vv, "[", "children"), "[", "rank")) == ott_rank
        if (any(success)) {
          dd <- do.call("rbind", sapply(vv, "[", "children"))[success, ]
          rownames(dd) <- unname(unlist(sapply(sapply(vv, "[", "children"), rownames)))[success]
          # rownames(dd) <- gsub("\\..*", "", rownames(dd))
          mm <- rbind(mm, dd)
        }
      }
      # utils::setTxtProgressBar(progression, i)
      # its easier to do the row naming in the previous steps, bc the following is much time consuming:
      # rownames(mm) <- unname(unlist(sapply(rotl::taxonomy_taxon_info(mm$ott_ids), "[", "unique_name")))
      all_children[[i]] <- mm
    }
  }
  names(all_children) <- names(input_ott_match)
  return(all_children)
}

#' Get an mrcaott tag from an OpenTree induced synthetic tree and get its name and ott id
#' @param tag A character vector with the mrca tag
#' @return A numeric vector with ott id from original taxon named with the corresponding ott name
#' @export
recover_mrcaott <- function(tag) {
  # for debugging:
  # tag <- "mrcaott7813ott454749"
  # tag = "mrcaott4291ott4292"
  # tag = "mrcaott27233ott44289"
  # tag = "mrcaott80776ott602508"
  if (!grep("mrcaott", tag)) {
    message("incorrect tag")
    return(NA)
  }
  mrca_ottids1 <- gsub("mrcaott", "", tag)
  mrca_ottids1 <- gsub("ott.*", "", mrca_ottids1)
  ott_ids <- c(mrca_ottids1, gsub("mrcaott.*ott", "", tag))
  mrca_lin <- datelife::get_ott_lineage(ott_ids = as.numeric(ott_ids))
  mm <- match(rownames(mrca_lin[[1]]), rownames(mrca_lin[[2]]))
  mm <- mm[!is.na(mm)]
  return(stats::setNames(mrca_lin[[2]][mm[1], "ott_ids"], rownames(mrca_lin[[2]])[mm[1]]))
}

# TODO: all these functions could be wrapped up into two or three single ones
# one to get one or all the lineages above a taxon
# another one to get some or all lineages below a taxon
# one to get everything: all_upper, all_lower, all (both upper and lower), or th eones specified by the users
# for example:
# get_ott_lineage(input = c("Felis", "Homo"), ott_ids = NULL, ott_rank = c("all", "all_upper", "all_lower"))
# get a list of all ranks available in ott
