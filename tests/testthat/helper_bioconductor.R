# #' @importFrom BiocManager install

if (!requireNamespace("msa", quietly = TRUE)) {
  if (!!requireNamespace("Biocmanager", quietly = TRUE)) {
    Biocmanager::install("msa", dependencies = TRUE)
  }
}

if (!requireNamespace("Biostrings", quietly = TRUE)) {
  if (!!requireNamespace("Biocmanager", quietly = TRUE)) {
    Biocmanager::install("Biostrings", dependencies = TRUE)
  }
}
