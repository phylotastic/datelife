## Resubmission
This is a resubmission where I have:

* Refactored code to remove dependencies of packages not on CRAN.
* Removed system requirement of mafft software. It was replaced with an R function.

### Test environments:

#### R 4.1.0, Platform: x86_64-apple-darwin17.0 (64-bit), local devtools::check()

0 errors | 0 warnings | 2 notes

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
    opentree_chronograms.rda object. We then ran tools::showNonASCII() on each element of the object. non-ASCII characters are found on the opentree_chronograms$authors list, which contains names of authors of our study database.

#### R 4.1.0, Platform: , devtools::check_win_release()


Results for rhub::check_for_cran():

https://builder.r-hub.io/status/datelife_0.6.0.tar.gz-54f3cbcfda7b44c488b2402364b19949
https://builder.r-hub.io/status/datelife_0.6.0.tar.gz-f4ec0e8a04c84a828fd2986e4f925141
https://builder.r-hub.io/status/datelife_0.6.0.tar.gz-5dc6c3acd693475ab5b4130f96eabaa2

 Possibly misspelled words in DESCRIPTION:
  Chronogram (2:46)
  DateLife (22:62)
  chronograms (20:73)
  workflows (20:26)


Found the following (possibly) invalid DOIs:
  DOI: doi.org/10.1101/782094
    From: inst/CITATION
    Message: Invalid DOI
  DOI: doi.org/10.5281/zenodo.593938
    From: inst/CITATION
    Message: Invalid DOI
