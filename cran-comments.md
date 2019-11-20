## Test environments
* local OS X install, R 3.6.0
* ubuntu 16.04.6 (on travis-ci), R 3.6.1
* win-builder (devel and release)

## R CMD check results
There were no ERRORs, WARNINGs or NOTEs.


## Downstream dependencies
This is the first release of this package


  #####################
  > devtools::check()
Updating datelife documentation
Warning: Version of roxygen2 last used with this package is 6.1.99.9001.  You only have version 6.1.1
Writing NAMESPACE
Loading datelife
Registered S3 method overwritten by 'compare':
  method           from
  print.comparison testthat
Registered S3 method overwritten by 'geiger':
  method            from
  unique.multiPhylo ape
Registered S3 method overwritten by 'httr':
  method                 from
  as.character.form_file crul
Warning in (function (dep_name, dep_ver = "*")  :
  Dependency package 'strap' not available.
Error: Dependency package(s) 'strap' not available.
Call `rlang::last_error()` to see a backtrace.
