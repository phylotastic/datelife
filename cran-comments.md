## Resubmission
This is a resubmission. In this version I have:

* Fixed packages not on CRAN.

* Removed system requirement of mafft software, which was replaced with an R function.

## Test environments
* local OS X install, R 4.1.0
* local ubuntu 16.04.6, R 3.6.1
* local windows, R 4.1.0
* win-builder (devel and release)

## R CMD check results
There are no ERRORs, and no WARNINGs

There are 2 NOTEs:

* checking installed package size ... NOTE
    installed size is 10.2Mb
    sub-directories of 1Mb or more:
      R      5.0Mb
      data   4.7Mb

Justification:
    The package is hosting a database (data/opentree_chronograms.rda), increasing the size of the data dir. This database is needed to run main functions.

* checking data for non-ASCII characters ... NOTE
  Note: found 2575 marked UTF-8 strings

Justification:
    We ran tools::showNonASCIIfile() on our data objects. The culprit is our
    opentree_chronograms.rda object. We then ran tools::showNonASCII() on each element of the object. non-ASCII characters are found on the opentree_chronograms$authors list, which contains names of authors of our study
    database.

## Downstream dependencies

This will be the first release of this package.
