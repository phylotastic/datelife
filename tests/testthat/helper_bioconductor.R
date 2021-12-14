if (!requireNamespace("msa", quietly = TRUE)) {
  if (!requireNamespace("Biocmanager", quietly = TRUE)) {
    Biocmanager::install("msa")
  }
}

if (!requireNamespace("Biostrings", quietly = TRUE)) {
  if (!requireNamespace("Biocmanager", quietly = TRUE)) {
    Biocmanager::install("Biostrings")
  }
}
