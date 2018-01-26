pkgname <- "datelife"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('datelife')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("datelife_search")
### * datelife_search

flush(stderr()); flush(stdout())

### Name: datelife_search
### Title: Core function to go from a vector of species, newick string, or
###   phylo object get a chronogram or dates back
### Aliases: datelife_search datelife

### ** Examples

ages <- datelife_search(c("Rhea americana", "Pterocnemia pennata", "Struthio camelus", "Mus musculus"), output.format="mrca")



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
