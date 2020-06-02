## Resubmission
This is a resubmission. In this version I have:

* Fixed packages not on CRAN.

* Removed system requirement of mafft software, it was replaced with an R function

## Test environments
* local OS X install, R 3.6.0
* ubuntu 16.04.6 (on travis-ci), R 3.6.1
* win-builder (devel and release)

## R CMD check results
There were no ERRORs, and no WARNINGs

There was 1 NOTE:

* checking installed package size ... NOTE
    installed size is 10.2Mb
    sub-directories of 1Mb or more:
      R      5.0Mb
      data   4.7Mb

    The package is hosting a database (data/opentree_chronograms.rda), increasing the size of the data dir. This database is needed to run main functions.
    The R dir appears of a much smaller size locally (351Kb).

## Downstream dependencies

This is the first release of this package.
