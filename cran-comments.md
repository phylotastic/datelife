## Resubmission
This is a resubmission where I have:

* Refactored code to remove dependencies of packages not on CRAN.
* Removed system requirement of mafft software. It was replaced with an R function.

### Test environments:

- MacOS, Platform: x86_64-apple-darwin17.0 (64-bit)
  - R 4.1.0, local devtools::check()
  - R 4.1.1, rhub::check_for_cran(platforms = "macos-highsierra-release-cran")
- Windows, Platform x86_64-w64-mingw32 (64-bit)
  - R 4.1.2, devtools::check_win_release()
  - R Under development (2021-11-26 r81252), rhub::check_for_cran(platforms = "windows-x86_64-release")
- Linux, Platform: x86_64-pc-linux-gnu (64-bit)
  - R 4.1.2, rhub::check_for_cran("ubuntu-gcc-release")
  - R Under development (2021-11-26 r81252), rhub::check_for_cran("linux-x86_64-fedora-clang")

### Results
0 errors | 0 warnings | 3 notes

* **checking installed package size ... NOTE**

    installed size is 7.8Mb

    sub-directories of 1Mb or more:

      data         4.7Mb

      figures      2.3Mb

_Comments_:

    The package is hosting a database (data/opentree_chronograms.rda), increasing the size of the data dir. This database is needed to run main functions.

* **checking data for non-ASCII characters ... NOTE**

  Note: found 2575 marked UTF-8 strings

_Comments_:

    We ran tools::showNonASCIIfile() on our data objects. non-ASCII characters can only be found in the opentree_chronograms.rda object. We then ran tools::showNonASCII() on each element of the object. non-ASCII characters are found only in the opentree_chronograms$authors list, which contains names of authors of our study database. These are obtained from a different database, and should be conserved in their original form.

* **Possibly misspelled words in DESCRIPTION:**

  Chronogram (2:46)

  DateLife (22:62)

  chronograms (20:73)

  workflows (20:26)

_Comments_:

    To my knowledge, none of these words are mispelled.
