# functions (from phunding package) to deal with open tree Taxonomy rotl functions, namely:
# rotl::taxonomy_taxon_info

#' Identifies, extracts and cleans children names from a taxonomy_taxon_info() output
#' @param taxon_info An output of rotl::taxonomy_taxon_info()
#' @param invalid A character vector of flags used by Open Tree Taxonomy to detect invalid taxon names
#' @return a list with valid children unique OTT names, ott_ids and ranks
#' @details eliminates all taxa that will give problems when trying to retrieve an induced subtree from otol.
#' @examples
#' # some suppressed taxa within "Gleicheniales" family of ferns:
#' rotl::taxonomy_taxon_info(3866048)
#' rotl::taxonomy_taxon_info(5262766)
#' @export
clean_taxon_info_children <- function(taxon_info, invalid = c("barren", "extinct", "uncultured", "major_rank_conflict", "incertae_sedis", "unplaced", "conflict", "environmental", "not_otu", "hidden", "hybrid")){
    # invalid <- c("BARREN", "EXTINCT", "UNCULTURED", "MAJOR_RANK_CONFLICT", "INCERTAE", "UNPLACED", "CONFLICT")
    # names(taxon_info[[2]])
    for (i in seq(taxon_info)){
      # length(tax_info[[i]][[1]]$children)
      # sapply(tax_info[[i]][[1]]$children, length)
      # length(sapply(sapply(tax_info[[i]][[1]]$children, "[", "flags"), function(y) unlist(tolower(y))))
      # sapply(sapply(tax_info[[i]][[1]]$children, "[", "flags"), function(y) unlist(tolower(y)))[1]
      ii <- lapply(sapply(sapply(taxon_info[[i]]$children, "[", "flags"), unlist), function(y) any(toupper(invalid) %in% y))
      # ii <- lapply(sapply(sapply(taxon_info[[i]][[1]]$children, "[", "flags"), unlist), function(y) any(toupper(invalid) %in% y))
      if(length(ii)>0){
        taxon_info[[i]]$children <- taxon_info[[i]]$children[!unlist(ii)]
      }
      # now clean names with sp. or environmental
      ii <- unlist(sapply(c("sp\\.","environmental", "unclassified", "incertae"), function(x) grep(x, sapply(sapply(taxon_info[[i]]$children, "[", "unique_name"), unlist))))

      if(length(ii)>0){
        taxon_info[[i]]$children <- taxon_info[[i]]$children[-unique(ii)]
      }
      # taxon_info[[i]][[1]]$children <- taxon_info[[i]][[1]]$children[!unlist(ii)]
    }
    return(taxon_info)
}

#' checks input for get_ott_clade and get_ott_children functions
#' returns a numeric vector of ott ids
#' @param input A character vector of names
#' @param ott_id A numeric vector of ott ids obtained with rotl::taxonomy_taxon_info or rolt::tnrs_match_names or datelife::batch_tnrs_match_names
#' @export
check_ott_input <- function(input, ott_id){
    if(is.null(ott_id)){
      if(!any(c("ott_id", "ott_ids") %in% names(input))){
        # input <- datelife::datelife_query_check(input)$cleaned_names
        input_tnrs <- datelife::batch_tnrs_match_names(names = input)
        # should we clean tnrs from invalid? what about NA's?
        # if we decide to clean, the two following lines should be uncommented:
        # df <- clean_tnrs(tnrs = df)
        # df <- df[!is.na(df$unique_name),]  # gets rid of names not matched with rotl::tnrs_match_names; otherwise rotl::tol_induced_subtree won't run
        input_ott_match <- suppressWarnings(as.numeric(input_tnrs$ott_id))
        names(input_ott_match) <- input_tnrs$unique_name
      } else {
        input_ott_match <- suppressWarnings(as.numeric(input$ott_id))
        if(is.data.frame(input)){
          names(input_ott_match) <- rownames(input)
        }
        # else {
        #     names(input_ott_match) <- names(input$ott_id)
        # }
      }

      if(any(is.na(input_ott_match))){
          message(paste0("\nInput '", paste(input[which(is.na(input_ott_match))], collapse = "', '"), "' not found in Open Tree of Life Taxonomy."))
          input_ott_match <- input_ott_match[which(!is.na(input_ott_match))]
      }
    } else {
          input_ott_match <- suppressWarnings(as.numeric(ott_id))
          # add a check for existing/valid ott ids??
      if(any(is.na(input_ott_match))){
          message(paste0("\nOtt id '", paste(ott_id[which(is.na(input_ott_match))], collapse = "', '"), "' not numeric and will be excluded from the search."))
          input_ott_match <- input_ott_match[which(!is.na(input_ott_match))]
      }
      # is the following enough check for valid ott ids? may be we should move it above the previous chunk
      names(input_ott_match) <- rotl::tax_name(rotl::taxonomy_taxon_info(ott_ids = input_ott_match))
    }
    if(length(input_ott_match) < 1){
      message("At least one valid input name or numeric ott_id are needed to get any information")
      return(NA)
    }
    return(input_ott_match)
}
#' gets the ott id and name of all lineages from one or more input taxa
#' @inheritParams check_ott_input
#' @return a list of named numeric vectors of ott ids from input and all the clades it belongs to.
#' @examples
#' taxa <- c("Homo", "Bacillus anthracis", "Apis", "Salvia")
#' lin <- get_ott_lineage(taxa)
#' @export
get_ott_lineage <- function(input, ott_id = NULL){
  input_ott_match <- check_ott_input(input, ott_id)
  tax_info <- .get_ott_lineage(input_ott_match)
  lin <- lapply(tax_info, "[", "lineage")
  ott_names <- sapply(lin, function(x) unlist(sapply(x[[1]], "[", "unique_name")))
  # unlist(sapply(lin[[1]][[1]], "[", "unique_name"))
  # length(ott_names) == length(tax_info)
  ott_ids <- sapply(seq(length(lin)), function(x) setNames(unlist(sapply(lin[[x]][[1]], "[", "ott_id")), ott_names[[x]]))
  setNames(ott_ids, names(input_ott_match))
}
#' Gets the lineage of a set of taxa using rotl:taxonomy_taxon_info(include_lineage = TRUE)
#' @param input_ott_match An Output of check_ott_input function.
#' @return A taxonomy_taxon_info object
.get_ott_lineage <- function(input_ott_match){
  tax_info <- vector(mode = "list", length = length(input_ott_match))
  progression <- utils::txtProgressBar(min = 0, max = length(tax_info), style = 3)
  for (i in seq(length(input_ott_match))){
      tax_info[i] <- tryCatch(rotl::taxonomy_taxon_info(input_ott_match[i], include_lineage = TRUE),
        error = function(e) NA)
      utils::setTxtProgressBar(progression, i)
  }
  tax_info
}


#' gets the ott id and name of one or several given taxonomic rank from one or more input taxa
#' @inheritParams check_ott_input
#' @inheritParams get_ott_children
#' @return a list of named numeric vectors with ott ids from input and all requested ranks
#' @export
get_ott_clade <- function(input = c("Felis", "Homo"), ott_id = NULL, ott_rank = "family"){
    # ott_id= c('649007', '782239', '782231', '1053057', '372826', '766272', '36015', '914245', '873016', '684051')
    # ott_id = c('431493', '431493', '431493', '431493', '431493', '431493', '431493', '429482', '429482', '429482')
  rank <- ott_rank
  input_ott_match <- check_ott_input(input, ott_id)
  tax_info <- .get_ott_lineage(input_ott_match)
  # names(tax_info[[10]])
  # sapply(tax_info[[10]]$lineage, "[", "rank")
  # length(tax_info[[10]]$lineage)
  # names(tax_info[[10]]$lineage[[1]])
  # tax_info[[10]]$lineage[[1]]$flags
  input_ott_names <- unlist(sapply(tax_info, "[", "unique_name"))
  rank_ott_ids <- rank_names <- lapply(seq(length(rank)), function(x) rep(NA, length(input_ott_names)))
  # I still need to drop all invalid lineages first here!!!
  for (i in seq(length(rank))){
      ready <- grepl(rank[i], unlist(sapply(tax_info, "[", "rank")))
      if(any(ready)){
        rank_names[[i]][ready] <- unlist(sapply(tax_info[ready], "[", "unique_name"))
        rank_ott_ids[[i]][ready] <- unlist(sapply(tax_info[ready], "[", "ott_id"))
      }
      if(!all(ready)){
        rank_names[[i]][!ready] <- sapply(tax_info[!ready], function(x) {
            lin <- tryCatch(x$lineage[grep(paste0("^", rank[i], "$"), sapply(x$lineage, "[", "rank"))][[1]]$unique_name,
            error = function(e) NA)
            return(lin)
        })
        rank_ott_ids[[i]][!ready] <- sapply(tax_info[!ready], function(x) {
            lin <- tryCatch(x$lineage[grep(paste0("^", rank[i], "$"), sapply(x$lineage, "[", "rank"))][[1]]$ott_id,
            error = function(e) NA)
            return(lin)
        })
      }
      names(rank_ott_ids[[i]]) <- rank_names[[i]]
      # length(rank_ott_ids[i])
      # stop()
  }
  names(input_ott_match) <-  input_ott_names
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
#' extracts valid children from a set of input names or ott ids (not from a taxonomy taxon info object)
#' @inheritParams check_ott_input
#' @export
#' @details
#' # genus Dictyophyllidites with ott id = 6003921 has only extinct children
#' # in cases like this the same name will be returned
#' tti <- rotl::taxonomy_taxon_info(6003921, include_children = TRUE)
#' gvc <- get_valid_children(ott_id = 6003921)
#' # More examples:
#' get_valid_children(ott_id = 769681) # Psilotopsida
#' get_valid_children(ott_id = 56601)  # Marchantiophyta
get_valid_children <- function(input = c("Felis", "Homo", "Malvaceae"), ott_id = NULL){
    input_ott_match <- check_ott_input(input, ott_id)
    all_children <- vector(mode = "list", length = length(input_ott_match))
    # monotypic <- vector(mode = "logical", length = length(input_ott_match))
    progression <- utils::txtProgressBar(min = 0, max = length(all_children), style = 3)
    for (i in seq(length(input_ott_match))){
        tt <- tryCatch(rotl::taxonomy_taxon_info(input_ott_match[i], include_children = TRUE),
          error = function(e) NA)
          # length(tt[[1]])
        tt <- clean_taxon_info_children(tt) # removes all invalid children
        if(length(tt[[1]]$children) > 0){
          # sapply(tt[[1]]$children, "[", "flags")
          # sapply(tt[[1]]$children, "[", "unique_name")
          # which(unlist(sapply(tt[[1]]$children, "[", "unique_name")) == "Mesangiospermae")
          # tt[[1]]$children[108]
          rr <- unname(unlist(sapply(tt[[1]]$children, "[", "rank")))
          # ii <- grep(paste0("^", ott_rank, "$"), unname(unlist(rr)))  # need to unlist rr
          child <- unname(unlist(sapply(tt[[1]]$children, "[", "ott_id")))
          # if(length(child)>0){
              names(child) <- names(rr) <- unname(unlist(sapply(tt[[1]]$children, "[", "unique_name")))
          # }
          monotypic <- FALSE
        } else {
          child <- tt[[1]]$ott_id
          rr <- tt[[1]]$rank
          names(child) <- names(rr) <- tt[[1]]$unique_name
          monotypic <- TRUE
        }
        all_children[[i]] <- list(children = data.frame(ott_id = child, rank = rr), is_monotypic = monotypic)
        utils::setTxtProgressBar(progression, i)
    }
    names(all_children) <- names(input_ott_match)
    return(all_children)
}

#' use this instead of rotl::tol_subtree when taxa are not in synthesis tree and you still need to get all species or an induced otol subtree
#' @inheritParams check_ott_input
#' @param ott_rank A character vector with the ranks you wanna get lineage children from.
#' @examples
#' tnrs <- rotl::tnrs_match_names("Gleicheniales")) # Polypodiales are missing too
#' \dontrun{
#'   rotl::tol_subtree(tnrs$ott_id[1])
#'   Error: HTTP failure: 400
#'   [/v3/tree_of_life/subtree] Error: node_id was not found (broken taxon).
#' }
#' goc1 <- get_ott_children(ott_id = tnrs$ott_id[1])
#' rownames(goc1[[1]])
#' tree_goc1 <- datelife::get_otol_synthetic_tree(goc1[[1]])
#' plot(tree_goc1, cex = 0.3)
#' # Other examples:
#' oo <- get_ott_children(input= "magnoliophyta", ott_rank = "order")
#' sum(oo$Magnoliophyta$rank == "order") # to know how many orders we got
#' @export
get_ott_children <- function(input = c("Felis", "Homo", "Malvaceae"), ott_id = NULL, ott_rank = "species"){
    input_ott_match <- check_ott_input(input, ott_id)
    all_children <- vector(mode = "list", length = length(input_ott_match))
    # progression <- utils::txtProgressBar(min = 0, max = length(all_children), style = 3)
    for (i in seq(length(input_ott_match))){
        mm <- data.frame(ott_id = vector(mode = "numeric", length = 0), rank = vector(mode = "logical", length = 0))
        vv <- get_valid_children(ott_id = input_ott_match[i])
        success <- vv[[1]]$children$rank == ott_rank | vv[[1]]$is_monotypic
        if(any(success)){
          mm <- rbind(mm, vv[[1]]$children[success,])
        }
        while(!all(success)){
          vv <- get_valid_children(ott_id = unlist(sapply(sapply(vv, "[", "children"), "[", "ott_id"))[!success])
          if(any(unlist(sapply(vv, "[", "is_monotypic")))){
            dd <- do.call("rbind", sapply(vv[unlist(sapply(vv, "[", "is_monotypic"))], "[", "children"))
            rownames(dd) <- unname(unlist(sapply(sapply(vv[unlist(sapply(vv, "[", "is_monotypic"))], "[", "children"), rownames)))
            # rownames(dd) <- gsub("\\..*", "", rownames(dd))
            mm <- rbind(mm, dd)
            vv <- vv[!unlist(sapply(vv, "[", "is_monotypic"))]
          }
          success <- unlist(sapply(sapply(vv, "[", "children"), "[", "rank")) == ott_rank
          if(any(success)){
            dd <- do.call("rbind", sapply(vv, "[", "children"))[success,]
            rownames(dd) <- unname(unlist(sapply(sapply(vv, "[", "children"), rownames)))[success]
            # rownames(dd) <- gsub("\\..*", "", rownames(dd))
            mm <- rbind(mm, dd)
          }
        }
        # utils::setTxtProgressBar(progression, i)
        # its easier to do the row naming in the previous steps, bc the following is much time consuming:
        # rownames(mm) <- unname(unlist(sapply(rotl::taxonomy_taxon_info(mm$ott_id), "[", "unique_name")))
        all_children[[i]] <- mm
    }
    names(all_children) <- names(input_ott_match)
    return(all_children)
}
