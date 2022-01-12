#' Taxon name resolution service (tnrs) applied to a vector of names by batches
#' @param input A character vector of taxon names, or a phylo object with tip names, to be matched to taxonomy.
#' @param reference_taxonomy A character vector specifying the reference taxonomy to use for tnrs.
#' @inheritDotParams rotl::tnrs_match_names -names
#' @details There is no limit to the number of names that can be queried and matched.
#' @return An object of class  data frame or phylo, with the added class match_names.
#' @export
tnrs_match <- function(input, reference_taxonomy, tip, ...) {
  UseMethod("tnrs_match", input)
}

#' @return \code{NULL}
#'
#' @rdname tnrs_match
#' @method tnrs_match default
#' @examples
#' tnrs_match(input = c("Mus"))
#' tnrs_match(input = c("Mus", "Mus musculus"))
#' tnrs_match(input = c("Mus", "Echinus", "Hommo", "Mus"))
#' @export
tnrs_match.default <- function(input, reference_taxonomy = "otl", ...) { # enhance: add other reference taxonomies in the future
  # names <- unique(names) # it is best to allow processing everything, i.e., without modifying the original input vector
  # names <- names[!is.na(names)]  # tnrs_match_names does not allow NAs, but they're caught after with tryCatch
  # for debugging:
  # input <- c("Mus", "Mus musculus")
  # input = c("cetaceae", "felis", "luke skywalker")
  input <- stringr::str_trim(input, side = "both") # cleans the input of lingering unneeded white spaces
  allinput <- input
  input <- unique(input)
  # enhance: infer taxonomic contexts:
  # tnrs_infer_context(names = names)
  # rotl::tnrs_contexts()
  if (reference_taxonomy == "otl") {
    df <- suppressWarnings(rotl::tnrs_match_names(names = "Apis mellifera")) # this just generates the dummy for the table, it will be removed at the end
    df <- df[nrow(df) + 1, ]
    df[nrow(df) + length(input) - 1, ] <- NA
    progression <- utils::txtProgressBar(min = 0, max = length(input), style = 3)
    for (i in seq(length(input))) {
      df[i, ] <- tryCatch(suppressMessages(rotl::tnrs_match_names(names = input[i], ...)),
        error = function(e) {
          no_match <- rep(NA, length(df[1, ]))
          no_match[1] <- input[i]
          no_match
        }
      )
      utils::setTxtProgressBar(progression, i)
    }
    cat("\n") # just to make the progress bar look better
    # hardcoding Mus:
    if (sum("mus" == tolower(input)) > 0) {
      df["mus" == tolower(input), ] <- list("mus", "Mus (genus in Deuterostomia)", FALSE, 1068778, FALSE, "SIBLING_HIGHER", 3)
    }
    rownames(df)[1] <- "1"
    # df[is.na(df$unique_name),1] <- input[is.na(df$unique_name)]  # in case the unmatched input are dropped from final df
    # df <- rotl::tnrs_match_names(names = input)
    # df has the same order as input
    # when some input are not matched it gives a warning: NAs introduced by coercion
    # but if no input is matched, it gives an error that needs to be caught
    ii <- match(tolower(allinput), df$search_string)
    result <- df[ii, ]
  }
  return(result) # returns the whole data frame
}
#' @return \code{NULL}
#'
#' @rdname tnrs_match
#' @method tnrs_match phylo
# #
# # ' Taxon name resolution service (tnrs) applied to tips of a phylogeny
# # ' @inheritParams phylo_check
#' @param tip A vector of mode numeric or character specifying the tips to match. If left empty all tips will be matched.
# #' @param reference_taxonomy A character vector specifying the reference taxonomy to use for tnrs.
# #' @inheritDotParams rotl::tnrs_match_names -names
# #' @return An object of class phylo and match_names. See details.
#' @details
#' The output will preserve all elements from original input phylo object and will add
#' \describe{
#'     \item{phy$mapped}{A character vector indicating the state of mapping of phy$tip.labels:}
#' 		\describe{
#' 		    \item{original}{Tnrs matching was not attempted. Original labeling is preserved.}
#' 		    \item{ott}{Matching was manually made by a curator in Open Tree of Life.}
#' 		    \item{tnrs}{Tnrs matching was attempted and successful with no approximate matching. Original label is replaced by the matched name.}
#' 		    \item{approximated}{Tnrs matching was attempted and successful but with approximate matching. Original labeling is preserved.}
#' 		    \item{unmatched}{Tnrs matching was attempted and unsuccessful. Original labeling is preserved.}
#' 		}
#'     \item{phy$original.tip.label}{A character vector preserving all original labels.}
#'     \item{phy$ott_ids}{A numeric vector with ott id numbers of matched tips. Unmatched and original tips will be NaN.}
#' }
#' if tips are duplicated, tnrs will only be run once (avoiding increases in function running time) but the result will be applied to all duplicated tip labels
#' @export
tnrs_match.phylo <- function(input, reference_taxonomy = "otl", tip = NULL, ...) { # we can add other reference taxonomies in the future
  # enhance_aproximates: add an argument in case we want to give the choice to users of changing only direct matches or also approximated matches
  phylo_check(input, dated = FALSE)
  phy.ori <- phy <- input
  if (missing(tip) | is.null(tip)) {
    tomaptip <- 1:ape::Ntip(phy)
  }
  if (is.character(tip)) {
    tomaptip <- match(tip, phy$tip.labels)
    # enhance: add a more general grep instead of match(). Test which one is faster and more accurate.
  }
  if (is.numeric(tip)) {
    tomaptip <- tip
  }
  if (any(tomaptip > ape::Ntip(phy)) | any(is.na(tomaptip))) {
    fails <- unique(c(which(is.na(tomaptip)), which(tomaptip > ape::Ntip(phy))))
    message(paste0("Values in tip argument do not correspond to any tip.labels in input: '", paste(tip[fails], collapse = "', '"), "'"))
    tomaptip <- tomaptip[-fails]
  }
  if (is.null(phy$mapped)) {
    phy$mapped <- rep("original", length(phy$tip.label))
  }
  if (is.null(phy$ott_ids) | length(phy$ott_ids) == 0) {
    phy$ott_ids <- rep(NaN, length(phy$tip.label))
  }
  new.names <- tnrs_match(input = unique(phy$tip.label[tomaptip]), reference_taxonomy = reference_taxonomy, ...) # a data.frame
  # not needed to clean new.names with clean_tnrs, bc we are leaving taxa matched with flags
  # new.names <- tnrs_match(c(unique(phy$tip.label[tomaptip]), "NotAtaxon")) # see what happens when we have something that does not match
  # in case we have NAs in approximate match (otherwise following conditionals would not work):
  new.names$approximate_match[is.na(new.names$approximate_match)] <- TRUE
  new.names$approximate_match <- as.logical(new.names$approximate_match)
  # update the tnrs_match data.frame in case there were duplicated taxa in tips:
  new.names <- new.names[match(tolower(phy$tip.label[tomaptip]), new.names$search_string), ]
  # update the mapped element:
  phy$mapped[tomaptip[is.na(new.names$unique_name)]] <- "unmatched"
  phy$mapped[tomaptip[new.names$approximate_match]] <- "approximated"
  # get the rows that have good matches and that are not approximate:
  matched <- !is.na(new.names$unique_name) & !new.names$approximate_match
  phy$mapped[tomaptip[matched]] <- "tnrs"
  # change phy$tip.labels to tnrs matched names
  # enchance_aproximates: following two lines of code are useful in case we want to give the choice to users of changing only direct matches or also approximated matches
  # phy$tip.label[tomaptip[matched]] <- new.names$unique_name[matched]
  # phy$tip.label[tomaptip[new.names$approximate_match]] <- new.names$unique_name[new.names$approximate_match]
  phy$tip.label[tomaptip[!is.na(new.names$unique_name)]] <- new.names$unique_name[!is.na(new.names$unique_name)]
  # add element of original tip labels
  phy$original.tip.label <- phy.ori$tip.label
  # add ott_ids
  phy$ott_ids[tomaptip[!is.na(new.names$unique_name)]] <- new.names$ott_id[!is.na(new.names$unique_name)]
  # add a class to phylo object
  class(phy) <- append(class(phy), "match_names")
  return(phy)
}
#' Eliminates unmatched (NAs) and invalid taxa from a [rotl::tnrs_match_names()] or [tnrs_match()] output
#' Useful to get ott ids to retrieve an induced synthetic Open Tree of Life.
#' Needed because using `include_suppressed = FALSE` in [rotl::tnrs_match_names()] does not drop all invalid taxa.
#'
#' @param tnrs A data frame, usually an output from datelife::tnrs_match or rotl::tnrs_match_names functions, but see details.
#' @param invalid A character string with flags to be removed from final object.
#' @param remove_nonmatches Boolean, whether to remove unsuccessfully matched names or not.
#' @details Input can be any data frame or named list that relates taxa stored in an element named "unique" to a validity category stored in "flags".
#' @return A data frame or named list (depending on the input) with valid taxa only.
#' @export
clean_tnrs <- function(tnrs, invalid = c("barren", "extinct", "uncultured", "major_rank_conflict", "incertae", "unplaced", "conflict", "environmental", "not_otu"),
                       remove_nonmatches = FALSE) {
  if (!"flags" %in% names(tnrs)) {
    message("tnrs should be a data.frame from datelife::tnrs_match or rotl::tnrs_match_names functions")
    if (!is.data.frame(tnrs)) {
      message("Or at least contain a flags element.")
    }
    return(tnrs)
  }
  if (length(unique(sapply(tnrs, length))) != 1) {
    message("elements in tnrs are of different length, check that out")
    return(tnrs)
  }
  tt <- as.data.frame(tnrs)
  inv <- sapply(1:nrow(tt), function(i) any(tolower(invalid) %in% tolower(tt$flags[i])))
  tt <- tt[!inv, ]
  if (remove_nonmatches) {
    tt <- tnrs[!is.na(tt$unique_name), ]
  }
  return(tt)
}

#' Get OTT ids from a character vector containing species names and OTT ids.
#' @param x A character vector of taxon names, or a phylo object with tip names containing OTT ids.
#' @param na.rm A logical value indicating whether `NA` values should be stripped from the output.
#' @return An object of class numeric containing OTT ids only.
#' @export
extract_ott_ids <- function(x, na.rm = TRUE) {
  UseMethod("extract_ott_ids", x)
}

#' @return \code{NULL}
#'
#' @rdname extract_ott_ids
#' @method extract_ott_ids default
#' @examples
#' canis <- rotl::tnrs_match_names("canis")
#' canis_taxonomy <- rotl::taxonomy_subtree(canis$ott_id)
#' my_ott_ids <- extract_ott_ids(x = canis_taxonomy$tip_label)
#' # get the problematic elements from input
#' canis_taxonomy$tip_label[attr(my_ott_ids, "na.action")]
#' @export
extract_ott_ids.default <- function(x, na.rm = TRUE) {
  res <- stringr::str_extract(x, "_ott\\d+")
  res <- gsub("_ott", "", res)
  res <- as.numeric(res)
  if (anyNA(res)) {
    message("After extracting ott ids, there are some non numeric elements:\n")
    message(paste(paste("\t", x[which(is.na(res))]), collapse = "\n"))
  }
  if (na.rm & anyNA(res)) {
    message("\nNAs removed.")
    res <- stats::na.omit(res)
  }
  return(res)
}
