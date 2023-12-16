## New submission v0.6.9
This is a new submission, where I have:
- added tests
- added spelling
- updated some functions


Thanks!

Luna Sanchez


## Test environments:

* local OS X install, R 4.2.2
* github actions: macOS Big Sur 10.16, x86_64-apple-darwin17.0, R 4.2.2
* github actions: Ubuntu 20.04.5, x86_64-pc-linux-gnu, release, old release and development, R 4.2.2
* github actions: Windows Server x64 (build 20348), x86_64-w64-mingw32, R 4.2.2
* R project win builder: Windows Server 2022, 64 bit, GCC, R under development 2023-06-17 r84564 ucrt
* rub: Fedora Linux, R-devel, clang, gfortran
* rub: Debian Linux, R-release, GCC


## R CMD check results

0 errors ✔ | 0 warnings ✔ | 3 notes ✖

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
Found the following (possibly) invalid DOIs:
  DOI: 10.1111/2041-210X.12051
    From: DESCRIPTION
    Status: Forbidden
    Message: 403
```

**Comments**: <br/>

The DOI is correct.
