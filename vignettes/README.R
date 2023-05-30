# to precompile vignettes

#knit the .orig vignette file into Rmd
# this generates usually two dirs at root, cache and figure
knitr::knit("vignettes/making_bold_trees.Rmd.orig",
            output = "vignettes/making_bold_trees.Rmd")

# move figures to vignettes folder
# needed for vignette compilation
sapply(list.files("figure/"), function(x){
  file.copy(from = paste0("figure/", x),
            to = "vignettes",
            overwrite = TRUE)
})

# cache is not needed for vignette compilation, just for rendering
# sapply(list.files("cache/"), function(x){
#   file.copy(from = paste0("cache/", x), to = "vignettes")
# })

# cache dir from root is ignored for build and git
# figure dir is deleted
# tip from http://theautomatic.net/2018/07/11/manipulate-files-r/
unlink("figure", recursive = TRUE)

#######################################################################
knitr::knit("data-raw/vignettes_making_bladj_trees.Rmd",
            output = "vignettes/making_bladj_trees.Rmd")
