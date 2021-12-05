
<!-- README.md is generated from README.Rmd. Please edit THIS file -->

<!-- badges: start -->

<!-- [![R build status](https://github.com/phylotastic/datelife/workflows/R-CMD-check/badge.svg)](https://github.com/phylotastic/datelife/actions) -->

![GitHub R package
version](https://img.shields.io/github/r-package/v/phylotastic/datelife?color=green)
[![Build
Status](https://travis-ci.org/phylotastic/datelife.svg)](https://travis-ci.org/phylotastic/datelife)
[![codecov](https://codecov.io/gh/phylotastic/datelife/branch/master/graph/badge.svg)](https://codecov.io/gh/phylotastic/datelife)
[![Github Open
Issues](https://img.shields.io/github/issues-raw/phylotastic/datelife.svg)](https://github.com/phylotastic/datelife/issues)
[![Github Closed
Issues](https://img.shields.io/github/issues-closed-raw/phylotastic/datelife.svg)](https://github.com/phylotastic/datelife/issues?q=is%3Aissue+is%3Aclosed)
[![DOI](https://zenodo.org/badge/23036/phylotastic/datelife.svg)](https://zenodo.org/badge/latestdoi/23036/phylotastic/datelife)
[![NSF-1458603](https://img.shields.io/badge/NSF-1458603-white.svg)](https://nsf.gov/awardsearch/showAward?AWD_ID=1458603)
[![NSF-0905606](https://img.shields.io/badge/NSF-0905606-white.svg)](https://nsf.gov/awardsearch/showAward?AWD_ID=0905606)
[![NSF-1458572](https://img.shields.io/badge/NSF-1458572-white.svg)](https://nsf.gov/awardsearch/showAward?AWD_ID=1458572)

<!-- badges: end -->

# Welcome to DateLife’s GitHub repository\!

Get a phylogenetic tree with branch lengths proportional to geologic
time (aka a ***chronogram***) of any two or more lineages of interest to
you: use the `datelife` R package or go to
[www.datelife.org](http://datelife.org/query/) to make a query of
chronograms available for your lineages in the [Open Tree of
Life](https://tree.opentreeoflife.org/curator)’s tree store.

You can also date a phylogenetic tree of your own making (or choosing
one from the literature), using node ages from chronograms queried with
`datelife` and the congruification method ([Eastman et
al. 2013](http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12051/abstract))
implemented with *geiger* ([Harmon et
al. 2008](http://bioinformatics.oxfordjournals.org/content/24/1/129.short)).

`datelife` has been developed as part of the
[phylotastic](http://phylotastic.org/) (NSF-funded) project, and is
still under rapid development.

  - [Installation](#installation)
  - [Quick intro](#quick-intro)
  - [Citation](#citation)
  - [Feedback and info for developers](#feedback)
  - [License](#license)

## Installation

Get the stable version from
[CRAN](https://cran.r-project.org/web/packages/datelife/index.html):

``` r
install.packages("datelife")
```

Or, install the development version from the GitHub repository:

``` r
devtools::install_github("phylotastic/datelife")
```

## Quick intro

You can update the cached objects with the `update_all_cached()`
function.

## Citation

If you use `datelife` for a publication, please cite the R package and
the accompanying paper:

<p>

O’Meara B, Sanchez-Reyes L, Eastman J, Heath T, Wright A, Schliep K,
Chamberlain S, Midford P, Harmon L, Brown J, Pennell M, Alfaro M (2021).
<em>datelife: Go from a List of Taxa or a Tree to a Chronogram using
Open Scientific Data</em>. R package version 0.5.0.

</p>

<p>

Sanchez-Reyes L, O’Meara B (2019). “datelife: Leveraging databases and
analytical tools to reveal the dated Tree of Life.” <em>bioRxiv</em>,
<b>782094</b>. doi:
<a href="https://doi.org/10.1101/782094">10.1101/782094</a>.

</p>

You can get these citations and the bibtex entry with:

``` r
citation("datelife")
toBibtex(citation("datelife"))
```

## Feedback and Information for Developers

All comments, ideas and questions about `datelife` are encouraged. You
are welcome to post an issue
[here](https://github.com/phylotastic/datelife/issues/new), or to make a
[pull request](https://github.com/phylotastic/datelife/pulls) if you
want to contribute with code directly.

Package and function documentation was generated with
[roxygen2](https://CRAN.R-project.org/package=roxygen2):

``` r
roxygen2::roxygenise()
```

We used the package [lintr]() to check for coding style:

``` r
lintr::lint_package()
```

Code coverage was tested with the package
[covr](https://CRAN.R-project.org/package=covr):

``` r
cov <- covr::package_coverage()

usethis::use_data(cov, overwrite = TRUE)
```

You can see an interactive report of testing coverage:

``` r
covr::report(cov)
```

And, find code with zero coverage:

``` r
covr::zero_coverage(cov)
```

To create the CRAN release, first run a `R CMD check` from your
terminal. For that, change directories to the one above your package:

``` bash
cd ../
```

Generate a tar ball for your package by running `R CMD build`:

``` bash
R CMD build datelife
```

Finally, run `R CMD check` on the tar ball that you just generated:

``` bash
R CMD check datelife_0.5.0.tar.gz
```

### `datelife`s hexsticker

We used the package `hexSticker`:

```
library(hexSticker)
s <- sticker(~plot(cars, cex=.5, cex.axis=.5, mgp=c(0,.3,0), xlab="", ylab=""),
          package="datelife", p_size=20, s_x=.8, s_y=.6, s_width=1.4, s_height=1.2,
          filename="inst/figures/baseplot.png")
          
```

### Rendering the vignettes

```
knitr::knit("vignettes/Getting_started_with_datelife.Rmd")
```

## License

This package is free and open source software, licensed under GPL.
