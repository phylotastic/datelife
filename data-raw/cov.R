## code to prepare `cov` dataset goes here

install.packages("covr")
library(covr)

cov <- covr::package_coverage()

usethis::use_data_raw("cov")

# usethis::use_data(cov, overwrite = TRUE)
