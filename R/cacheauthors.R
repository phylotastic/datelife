
#' Create a cache  from Open Tree of Life
#' @param outputfile Path including file name
#' @return List containing author and curator results
#' @export
make_contributor_cache <- function(outputfile = "contributorcache.RData") {
  all.sources <- rotl::tol_about(include_source_list = TRUE)
  all.studies <- as.list(rep(NA, length(all.sources$source_list)))
  author.results <- data.frame()
  curator.results <- data.frame()
  missed_doi <- c()
  for (i in sequence(length(all.studies))) {
    study.id <- strsplit(all.sources$source_list[[i]], "@")[[1]][[1]]
    try(all.studies[[i]] <- rotl::get_study_meta(study.id)) # some studies have no metadata
    doi <- authors <- clade.name <- curator <- study.year <- NULL
    try(doi <- gsub("https?://(dx\\.)?doi.org/", "", attr(rotl::get_publication(all.studies[[i]]), "DOI")))
    if (length(doi) == 0) {
      missed_doi <- c(missed_doi, study.id)
      warning(paste(study.id, "has no DOI attribute, author names will not be retrieved.")) # where else could we obtain a study's doi???
    } else {
      try(authors <- as.character(knitcitations::bib_metadata(doi)$author))
      if (is.null(authors)) {
        warning(paste(study.id, "DOI page was not found, author names will not be retrieved."))
        missed_doi <- c(missed_doi, study.id)
      }
    }
    # The following is necessary to establish good encoding in author.results object:
    if (!is.null(authors) & length(authors) > 0) {
      Encoding(authors) <- "latin1"
      authors <- iconv(authors, "latin1", "UTF-8")
    }
    try(clade.name <- all.studies[i][[1]]$nexml$`^ot:focalCladeOTTTaxonName`)
    try(curators <- unique(rotl::get_study_meta(study.id)[["nexml"]][["^ot:curatorName"]][[1]]))
    if (!is.null(authors)) {
      # author.results <- rbind(author.results, expand.grid(author=authors, study=study.id, clade=clade.name, publication=get_publication(all.studies[[i]])[1], doi=attr(get_publication(all.studies[[i]]),"DOI")))
      author.results <- rbind(author.results, expand.grid(person = authors, study = study.id, clade = clade.name, stringsAsFactors = FALSE))
    }
    if (!is.null(curators)) {
      # author.results <- rbind(author.results, expand.grid(author=authors, study=study.id, clade=clade.name, publication=get_publication(all.studies[[i]])[1], doi=attr(get_publication(all.studies[[i]]),"DOI")))
      curator.results <- rbind(curator.results, expand.grid(person = curators, study = study.id, clade = clade.name, stringsAsFactors = FALSE))
    }
    message(i, " of ", length(all.studies), " done ")
  }
  message("Authors from studies with the following id's could not be retrieved:")
  for (i in seq_len(missed_doi)) {
    message("\t", missed_doi[i])
  }
  author.pretty <- make_overlap_table(author.results)
  try(Encoding(author.pretty$person) <- "latin1")
  try(author.pretty <- iconv(author.pretty$person, "latin1", "UTF-8"))
  curator.pretty <- make_overlap_table(curator.results)
  try(Encoding(curator.pretty$person) <- "latin1")
  try(curator.pretty <- iconv(curator.pretty$person, "latin1", "UTF-8"))
  save(author.results, curator.results, author.pretty, curator.pretty, missed_doi, file = outputfile)
  return(list(author.results = author.results, curator.results = curator.results, author.pretty = author.pretty, curator.pretty = curator.pretty, missed_doi = missed_doi))
}


#' Create a cache from TreeBase
#' @param outputfile Path including file name
#' @return List containing author and curator results
make_treebase_cache <- function(outputfile = "treebasecache.RData") { # nocov start
  InvertNames <- function(author) {
    return(paste(sapply(strsplit(author, ", "), rev), collapse = " "))
  }
  all.studies <- treebase::download_metadata("", by = "all")
  unlist(all.studies[[1]])[which(names(unlist(all.studies[[1]])) == "creator")]
  author.results <- data.frame()
  for (i in sequence(length(all.studies))) {
    authors <- study.id <- NULL
    try(authors <- unlist(all.studies[[i]])[which(names(unlist(all.studies[[i]])) == "creator")])
    try(study.id <- unlist(all.studies[[i]])["identifier"])
    if (!is.null(authors)) {
      authors <- sapply(authors, InvertNames)
      author.results <- rbind(author.results, expand.grid(person = authors, study = study.id, stringsAsFactors = FALSE))
    }
    message(paste(i, "of", length(all.studies), "done."))
  }
  tb.author.pretty <- make_overlap_table(author.results)[, -3]
  tb.author.results <- author.results
  save(tb.author.results, tb.author.pretty, file = outputfile)
  return(list(tb.author.results = tb.author.results, tb.author.pretty = tb.author.pretty))
} # nocov end

#' Create an overlap table
#' @param results_table An "author.results" or "curator.results" `data.frame`
#' @return A `data.frame` with information on curators and what clades they've worked on
#' @export
make_overlap_table <- function(results_table) {
  unique.person <- unique(results_table$person)
  final.table <- data.frame()
  for (person.index in sequence(length(unique.person))) {
    local.df <- subset(results_table, results_table$person == unique.person[person.index])
    clades <- table(as.character(local.df$clade))
    names(clades)[which(nchar(names(clades)) == 0)] <- "Unknown"
    clade.string <- paste(sort(paste0(names(clades), " (", clades, ")")), collapse = ", ")
    person.split <- strsplit(unique(local.df$person), " ")[[1]]
    last.name.guess <- person.split[length(person.split)]
    final.table <- rbind(final.table, data.frame(lastname = last.name.guess, person = unique.person[person.index], study.count = nrow(local.df), clades = clade.string, stringsAsFactors = FALSE))
  }
  final.table <- final.table[order(final.table$lastname), ]
  final.table <- final.table[, -1]
  return(final.table)
}

#' Associate TreeBase authors with studies
#' @return `data.frame` with author last name, author first and other names, and comma delimited URLs for TreeBase studies
#'
make_treebase_associations <- function() {
  all.studies <- treebase::download_metadata("", by = "all")
  unlist(all.studies[[1]])[which(names(unlist(all.studies[[1]])) == "creator")]
  author.results <- data.frame()
  for (i in sequence(length(all.studies))) {
    authors <- study.id <- NULL
    try(authors <- unlist(all.studies[[i]])[which(names(unlist(all.studies[[i]])) == "creator")])
    try(study.id <- unlist(all.studies[[i]])["identifier"])
    try(study.url <- gsub("purl.org/phylo/treebase/phylows/study/TB2:S", "https://treebase.org/treebase-web/search/study/summary.html?id=", unlist(all.studies[[i]])["identifier"]))

    try(publisher <- unlist(all.studies[[i]])["publisher"])
    try(year <- unlist(all.studies[[i]])["date"])

    if (!is.null(authors)) {
      author.results <- rbind(author.results, expand.grid(person = authors, url = paste0('<a href="', study.url, '">', year, " ", publisher, " (TreeBase)</a>"), stringsAsFactors = FALSE))
    }
    message(paste(i, "of", length(all.studies), "done."))
  }

  author.lumped <- stats::aggregate(author.results$url, by = list(author = author.results$person), FUN = paste, collapse = ", ")
  colnames(author.lumped) <- c("person", "urls")
  return(author.lumped)
}

#' Associate Open Tree of Life authors with studies
#' @return `data.frame` with author last name, author first and other names, and comma delimited URLs for OToL studies
make_otol_associations <- function() {
  # .res <- system("curl -X POST https://api.opentreeoflife.org/v3/studies/find_studies -H 'content-type:application/json'", intern=TRUE)
  # .res <- httr::POST("https://api.opentreeoflife.org/v3/studies/find_studies")
  #  res <- vapply(.res[["matched_studies"]], function(x) x[["ot:studyId"]],character(1))

  all.trees <- rbind(rotl::studies_find_trees(property = "is_deprecated",
                                              value = "false",
                                              verbose = TRUE,
                                              detailed = TRUE),
                     rotl::studies_find_trees(property = "is_deprecated",
                                              value = "true",
                                              verbose = TRUE,
                                              detailed = TRUE)
  )
  # trees <- list()
  authors <- list()
  curators <- list()
  studies <- list()
  dois <- list()
  years <- list()
  tree.count <- 0
  bad.ones <- c()
  studyids <- list()

  for (study.index in sequence(dim(all.trees)[1])) {
    study.id <- all.trees$study_ids[study.index]
    doi <- NULL
    try(doi <- gsub("https?://(dx\\.)?doi.org/", "", attr(rotl::get_publication(rotl::get_study_meta(study.id)), "DOI")))
    authors <- append(authors, NA)
    dois <- append(dois, NA)
    years <- append(years, NA)
    studyids <- append(studyids, NA)
    try(studyids[length(studyids)] <- study.id)
    try(years[length(years)] <- rotl::get_study_meta(study.id)$nexml$`^ot:studyYear`)
    if (length(doi) == 0) {
      warning(paste(study.id, "has no DOI attribute, author names will not be retrieved."))
    } else {
      try(authors[length(authors)] <- list(paste(as.character(knitcitations::bib_metadata(doi)$author))))
      try(dois[length(dois)] <- doi)
    }
    try(studies <- append(studies, study.id))
  }
  for (i in sequence(length(authors))) {
    if (!any(is.na(authors[[i]])) & length(authors[[i]]) > 0) {
      Encoding(authors[[i]]) <- "latin1"
      authors[[i]] <- iconv(gsub("^\\s+|\\s+$", "", authors[[i]]), "latin1", "UTF-8")
    }
  }
  studies <- unlist(studies)
  combined.df <- data.frame()
  invert_names <- function(x) {
    x.split <- strsplit(x, " ")[[1]]
    final.names <- paste0(gsub("^\\s+|\\s+$", "", x.split[length(x.split)]), ", ", paste0(x.split[-length(x.split)], collapse = " "))
    return(final.names)
  }
  for (i in sequence(length(authors))) {
    urls <- paste0("<a href='https://tree.opentreeoflife.org/curator/study/view/", studyids[[i]], "'>", years[[i]], " ", studyids[[i]], " (OpenTree)</a>")
    for (author.index in sequence(length(authors[[i]]))) {
      combined.df <- rbind(combined.df, data.frame(person = invert_names(authors[[i]][author.index]), url = urls, stringsAsFactors = FALSE))
    }
  }
  author.lumped <- stats::aggregate(combined.df$url, by = list(author = combined.df$person), FUN = paste, collapse = ", ")
  colnames(author.lumped) <- c("person", "urls")
  return(author.lumped)
}

#' Find all authors and where they have deposited their trees
#' @param outputfile Path including file name. NULL to prevent saving.
#' @return a `data.frame` of "person" and "urls".
#' @export
make_all_associations <- function(outputfile = "depositorcache.RData") {
  tb <- make_treebase_associations()
  ot <- make_otol_associations()
  combined.df <- rbind(tb, ot)
  author.lumped <- stats::aggregate(combined.df$urls, by = list(author = combined.df$person), FUN = paste, collapse = ", ")
  colnames(author.lumped) <- c("person", "urls")
  author.lumped$person <- gsub("^\\s+|\\s+$", "", author.lumped$person)
  author.lumped <- author.lumped[order(author.lumped$person), ]
  depositors <- author.lumped
  if (!is.null(outputfile)) {
    save(depositors, file = outputfile)
  }
  return(author.lumped)
}
