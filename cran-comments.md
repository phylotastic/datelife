## Resubmission
This is a resubmission where I have:

* Fixed a (possibly) invalid file URI: data-raw/hexsticker-current.R From: README.md
Fixed to https://github.com/phylotastic/datelife/blob/master/data-raw/hexsticker-current.R

* Added donrun to examples with CPU (user + system) or elapsed time > 5s

* Checked possibly misspelled words in DESCRIPTION and determined they are not misspellings.

### Test environments:

- MacOS, Platform: x86_64-apple-darwin17.0 (64-bit)
  - R 4.1.0, local `devtools::check()`
  - R 4.1.1, `rhub::check_for_cran(platforms = "macos-highsierra-release-cran")`
- Windows, Platform x86_64-w64-mingw32 (64-bit)
  - R 4.1.2, `devtools::check_win_release()`
  - R Under development (2021-11-26 r81252), `rhub::check_for_cran(platforms = "windows-x86_64-release")`
- Linux, Platform: x86_64-pc-linux-gnu (64-bit)
  - R 4.1.2, `rhub::check_for_cran("ubuntu-gcc-release")`
  - R Under development (2021-11-26 r81252), `rhub::check_for_cran("linux-x86_64-fedora-clang")`

### Results
0 errors | 0 warnings | 3 notes

* Note 1

```
checking installed package size ... NOTE
    installed size is 6.1Mb
    sub-directories of 1Mb or more:
      data         4.7Mb
```
**Comments**: <br/>
The package is hosting a database (`data/opentree_chronograms.rda`), increasing the size of the data dir. This database is needed to run main functions.

* Note 2

```
checking data for non-ASCII characters ... NOTE
  Note: found 2575 marked UTF-8 strings
```

**Comments**: <br/>
We ran `tools::showNonASCIIfile()`` on our data objects. non-ASCII characters can only be found in the opentree_chronograms.rda object. We then ran `tools::showNonASCII()`` on each element of the object. non-ASCII characters are found only in the `opentree_chronograms$authors` list, which contains names of authors of our study database. These are obtained from a different database, and should be conserved in their original form.

* Note3

```
Possibly misspelled words in DESCRIPTION:
  Chronogram (2:46)
  DateLife (22:62)
  chronograms (20:73)
  workflows (20:26)
  phylogenetic (20:49)
```

**Comments**: <br/>
To my knowledge, none of these words are misspelled.
