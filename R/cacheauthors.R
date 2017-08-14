
#' Create a cache  from Open Tree of Life
#' @param outputfile Path including file name
#' @param verbose If TRUE, give status updates to the user
#' @return List containing author and curator results
#' @export
CreateContributorCache <- function(outputfile="contributorcache.rda", verbose=TRUE) {
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
  author.pretty <- CalculateOverlapTable(author.results)
  curator.pretty <- CalculateOverlapTable(curator.results)
  save(author.results, curator.results, author.pretty, curator.pretty, file=outputfile)
  return(list(author.results=author.results, curator.results=curator.results, author.pretty=author.pretty, curator.pretty=curator.pretty))
}


#' Create a cache from TreeBase
#' @param outputfile Path including file name
#' @param verbose If TRUE, give status updates to the user
#' @return List containing author and curator results
#' @export
CreateTreeBaseCache <- function(outputfile="treebasecache.rda", verbose=TRUE) {
  InvertNames <- function(author) {
    return(paste(sapply(strsplit(author, ', '), rev), collapse=" "))
  }
  all.studies <- treebase::download_metadata("", by="all")
  unlist(all.studies[[1]])[which(names( unlist(all.studies[[1]])) == "creator")]
  author.results <- data.frame()
  for (i in sequence(length(all.studies))) {
    authors <- study.id <- NULL
    try(authors <- unlist(all.studies[[i]])[which(names( unlist(all.studies[[i]])) == "creator")])
    try(study.id <- unlist(all.studies[[i]])["identifier"])
    if(!is.null(authors)) {
      authors <- sapply(authors, InvertNames)
      author.results <- rbind(author.results, expand.grid(person=authors, study=study.id, stringsAsFactors=FALSE))
    }
    if(verbose) {
      print(paste(i, "of", length(all.studies), "done"))
    }
  }
  tb.author.pretty <- CalculateOverlapTable(author.results)[,-3]
  tb.author.results <- author.results
  save(tb.author.results, tb.author.pretty, file=outputfile)
  return(list(tb.author.results=tb.author.results, tb.author.pretty=tb.author.pretty))
}

#' Create an overlap table
#' @param results.table An author.results or curator.results data.frame
#' @return Data.frame with info on people and what clades they've worked on
#' @export
CalculateOverlapTable <- function(results.table) {
  unique.person <- unique(results.table$person)
  final.table <- data.frame()
  for (person.index in sequence(length(unique.person))) {
    local.df <- subset(results.table, results.table$person==unique.person[person.index])
    clades <- table(as.character(local.df$clade))
    names(clades)[which(nchar(names(clades))==0)] <- "Unknown"
    clade.string <- paste(sort(paste0(names(clades), " (", clades, ")")), collapse=", ")
    person.split <- strsplit(unique(local.df$person), " ")[[1]]
    last.name.guess <- person.split[length(person.split)]
    final.table <- rbind(final.table, data.frame(lastname = last.name.guess, person=unique.person[person.index], study.count = nrow(local.df), clades=clade.string, stringsAsFactors=FALSE))
  }
  final.table <- final.table[order(final.table$lastname),]
  final.table <- final.table[,-1]
  return(final.table)
}
