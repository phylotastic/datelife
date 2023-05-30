## New submission v0.6.7
This is a new submission, package was archived on CRAN after because of an unavailable dependecy, `phylocomr`.
The dependecy is available again, and I have made a couple updates to the package:

* added vignette for bold data workflow (precomputed, so cache folder is buildignored)
* added some warnings when using random age for the root


Thanks!

Luna Sanchez


## Test environments:

* local OS X install, R 4.2.2
* macOS Big Sur 10.16, x86_64-apple-darwin17.0, R 4.2.2
* Ubuntu 20.04.5, x86_64-pc-linux-gnu (github actions), R 4.2.2
* Windows Server x64 (build 20348), x86_64-w64-mingw32 (github actions), R 4.2.2
* Fedora Linux, R-devel, clang, gfortran
* Ubuntu Linux 20.04.1 LTS, R-release, GCC
* Windows Server 2022, R-devel, 64 bit
* Debian Linux, R-release, GCC


## R CMD check results

0 errors ✔ | 0 warnings ✔ | 4 notes ✖

* Note 1

```
❯ checking package dependencies ... NOTE
  Packages suggested but not available for checking: 'msa', 'Biostrings'

  Imports includes 22 non-default packages.
  Importing from so many packages makes the package vulnerable to any
  of them becoming unavailable.  Move as many as possible to Suggests
  and use conditionally.
```
**Comments**: <br/>

All the packages imported are needed for the datelife workflow.


* Note 2

```
Possibly misspelled words in DESCRIPTION:
  BLADJ (32:15)
  Britton (33:6)
  Criscuolo (24:18)
  DateLife (21:62)
  Huelsenbeck (33:71)
  O'Meara (35:9)
  PATHd (32:79)
  Ronquist (34:9)
  SDM (23:67)
  Schenk (29:50)
  al (24:31, 26:61, 31:19, 32:30, 33:17)
  chronogram (22:61)
  chronograms (19:43, 23:14, 25:26)
  congruification (30:76)
  et (24:28, 26:58, 31:16, 32:27, 33:14)
  mrBayes (33:62)
  phylogenetic (19:62, 28:13)
  treePL (34:69)
```
**Comments**: <br/>

All these words are OK.

* Note 3

```
CRAN repository db overrides:
  X-CRAN-Comment: Archived on 2023-02-24 as requires archived package
    'phylocomr'.
```

**Comments**: <br/>

Package phylocomr is available on CRAn again.

* Note 4

```
Found the following (possibly) invalid DOIs:
  DOI: 10.1111/2041-210X.12051
    From: DESCRIPTION
    Status: Forbidden
    Message: 403
```

**Comments**: <br/>

The DOI is correct.
