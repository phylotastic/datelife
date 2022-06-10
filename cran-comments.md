## New submission v 0.6.4
This is a new submission where I have:

* addressed note on package size: The package hosts a database of 4.7Mb that is needed to run main functions
* Removed unnecessary "if (interactive())" from examples
* Fixed R code problem by useing `inherits` instead of if() conditions comparing class() to string
* Fixed examples that were failing in some Linux environments.

### Test environments:

- MacOS, Platform: x86_64-apple-darwin17.0 (64-bit)
  - R 4.2.0, local `devtools::check()`
  - R 4.1.1, `rhub::check_for_cran(platforms = "macos-highsierra-release-cran")`
- Windows, Platform x86_64-w64-mingw32 (64-bit)
  - R 4.2.0, `devtools::check_win_release()`
  - R Under development (2022-03-23 r81968 ucrt)
  - `rhub::check_for_cran(platforms = "windows-x86_64-release")`
- Linux, Platform: x86_64-pc-linux-gnu (64-bit)
  - R 4.2.0 release
  - R 4.3 Under development (2022-06-08 r82470)

### Results
0 errors ✔ | 0 warnings ✔ | 4 notes ✖

* Note 1

```
checking installed package size ... NOTE
    installed size is 6.1Mb
    sub-directories of 1Mb or more:
      data         4.7Mb
```

**Comments**: <br/>
The package is hosting a database (`data/opentree_chronograms.rda`), increasing the size of the data dir.
This database is needed to run main functions.


* Note 2

```
❯ checking package dependencies ... NOTE
  Packages suggested but not available for checking: 'msa', 'Biostrings'

  Imports includes 21 non-default packages.
  Importing from so many packages makes the package vulnerable to any
  of them becoming unavailable.  Move as many as possible to Suggests
  and use conditionally.
```
**Comments**: <br/>
All the packages imported are needed for the datelife workflow.

* Note 3

```
Possibly misspelled words in DESCRIPTION:
  Ané (26:53)
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
  al (24:31, 26:60, 31:19, 32:30, 33:17)
  chronogram (22:55)
  chronograms (19:43, 23:14, 25:26)
  congruification (30:76)
  et (24:28, 26:57, 31:16, 32:27, 33:14)
  mrBayes (33:62)
  treePL (34:69)
```
**Comments**: <br/>
None of these word are misspelled

* Note 4

```
Found the following (possibly) invalid DOIs:
  DOI: 10.1111/2041-210X.12051
    From: DESCRIPTION
    Status: Service Unavailable
    Message: 503
```
**Comments**: <br/>
The DOI is functional
