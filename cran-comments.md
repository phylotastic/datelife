## New submission v 0.6.2
This is a new submission where I have:

* Fixed examples that were failing in some Linux environments.

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
0 errors ✔ | 0 warnings ✔ | 1 note ✖

* Note 1

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
