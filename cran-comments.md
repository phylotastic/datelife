## New submission v0.6.7
This is a new submission where I have:

* added vignettes
*

Thanks!

Luna Sanchez


## Test environments:

* local OS X install, R 4.2.0
* macOS Big Sur 10.16, x86_64-apple-darwin17.0, R 4.2.2
* Ubuntu 20.04.5, x86_64-pc-linux-gnu (github actions), R 4.2.2
* Windows Server x64 (build 20348), x86_64-w64-mingw32 (github actions), R 4.2.2
* Fedora Linux, R-devel, clang, gfortran
* Ubuntu Linux 20.04.1 LTS, R-release, GCC
* Windows Server 2022, R-devel, 64 bit
* Debian Linux, R-release, GCC


## R CMD check results

0 errors ✔ | 0 warnings ✔ | 2 notes ✖

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
Found the following (possibly) invalid DOIs:
  DOI: 10.1111/2041-210X.12051
    From: DESCRIPTION
    Status: Service Unavailable
    Message: 503
```
**Comments**: <br/>

The DOI is functional
