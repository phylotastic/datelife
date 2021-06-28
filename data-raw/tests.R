# file from Flores-Abreu 2018
# install_github("fmichonneau/phyloch")
agave <- phyloch::read.beast("data-raw/Agave5_14Nov13.TreeAnnOut_nombresPub")
names(agave)

input <- c("Panthera leo",
           "Leopardus pardalis",
           "Puma concolor")

tree <- rotl::tol_induced_subtree(rotl::ott_id(rotl::tnrs_match_names(input)), label = "name")

datelife::datelife_search(input = c("Panthera_leo", "Puma_concolor", "Leopardus_pardalis"))
