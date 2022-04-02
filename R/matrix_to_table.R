#' Go from a patristic distance matrix to a node ages table
#'
#' @param matrix A patristic distance matrix.
#' @param reference A character vector with the study reference from where the ages come from.
#' @return A `data.frame` of "taxonA", "taxonB", and "age".
#' @export
matrix_to_table <- function(matrix, reference) {
  tt <- as.data.frame(as.table(matrix), useNA = "no", stringsAsFactors = FALSE)
  # create an empty vector length equal to number of columns
  new_colnames <- c()
  # assign column names:
  for (i in seq(ncol(tt))) {
    if (inherits(tt[[i]], "character")) {
      new_colnames <- c(new_colnames, paste0("taxon", LETTERS[i]))
    }
    if (inherits(tt[[i]], "numeric")) {
      new_colnames <- c(new_colnames, "nodeAge")
    }
  }
  colnames(tt) <- new_colnames
  # remove rows that are tips (with an age of 0)
  tt <- tt[tt$nodeAge != 0, ]
  # get the actual node age of the taxon pair by dividing the age by two
  tt$nodeAge <- tt$nodeAge/2
  # add the reference column
  tt$reference <- rep(reference, nrow(tt))
  # return the table
  tt
}

#' Go from a list of patristic distance matrix to a table of node ages
#'
#' @param matrices A names list of patristic distance matrices. Names correspond to the study reference.
#' @return A single `data.frame` of "taxonA", "taxonB", and "age".
#' @importFrom data.table rbindlist
#' @export
matrices_to_table <- function(matrices) {
  # apply matrix_to_table to a list of matrices and get a list of tables back:
  tt <- lapply(seq(length(matrices)),
               function(i) matrix_to_table(reference = names(matrices)[i],
                                           matrix = matrices[[i]]))
  # merge the list of tables and return a single table:
  data.table::rbindlist(tt)
}
