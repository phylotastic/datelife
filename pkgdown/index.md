
<!-- index.md is generated from index.Rmd. Make sure to edit the .Rmd file and not the .md -->

<img src='https://github.com/phylotastic/datelife/raw/master/man/figures/datelife-hexsticker-ai.png' align='right' style='width:150px' />

<!--
# BADGES DEV
install.packages("badger")
library(badger)
badge_github_actions("phylotastic/datelife")
badge_cran_checks("phylotastic/datelife")
-->
<!-- badges: start -->
<!-- Stable status -->

[![CRAN
status](https://www.r-pkg.org/badges/version/datelife)](https://CRAN.R-project.org/package=datelife)
[![CRAN RStudio mirror
downloads](https://cranlogs.r-pkg.org/badges/datelife)](https://www.r-pkg.org/pkg/datelife)
[![DOI](https://zenodo.org/badge/23036/phylotastic/datelife.svg)](https://zenodo.org/badge/latestdoi/23036/phylotastic/datelife)

<!-- Development status -->

[![GitHub master branch package
version](https://img.shields.io/github/r-package/v/phylotastic/datelife/master?color=y&label=GitHub%40master)](https://github.com/phylotastic/datelife)
[![GitHub R-CMD-check
Status](https://github.com/phylotastic/datelife/workflows/R-CMD-check/badge.svg)](https://github.com/phylotastic/datelife/actions/)
[![codecov](https://codecov.io/gh/phylotastic/datelife/branch/master/graph/badge.svg)](https://app.codecov.io/gh/phylotastic/datelife)
[![Github Open
Issues](https://img.shields.io/github/issues-raw/phylotastic/datelife.svg)](https://github.com/phylotastic/datelife/issues)
[![Github Closed
Issues](https://img.shields.io/github/issues-closed-raw/phylotastic/datelife.svg)](https://github.com/phylotastic/datelife/issues?q=is%3Aissue+is%3Aclosed)

<!-- Funding -->

[![NSF-1458603](https://img.shields.io/badge/NSF-1458603-white.svg)](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1458603)
[![NSF-0905606](https://img.shields.io/badge/NSF-0905606-white.svg)](https://www.nsf.gov/awardsearch/showAward?AWD_ID=0905606)
[![NSF-1458572](https://img.shields.io/badge/NSF-1458572-white.svg)](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1458572)

<!-- badges: end -->

# Welcome to DateLife’s R package documentation!

## What is `datelife`?

`datelife` is an R package that allows researchers and the general
audience to obtain open scientific data on the age of any organism they
are interested in, by retrieving organism ages from a database of dated
phylogenetic trees (*aka* chronograms), that have been peer-reviewed and
published as part of a scientific research article, in an indexed
journal ([Open Tree of Life’s tree
store](https://tree.opentreeoflife.org/curator)). As such, organism ages
retrieved by `datelife` constitute state-of-the-art, peer-reviewed,
public scientific knowledge, that can be accessed and reused by experts
and non-experts in the field alike.

## How can you use `datelife`?

You can install the `datelife` R package on your own computer and use it
locally.

If you do not want/have time to deal with installation and R code, you
can use [DateLife’s interactive website application](). However, the
website is not up at the moment, apologies.

<!-- http://datelife.opentreeoflife.org/query/ -->

## Documentation topics:

- [Local
  installation](http://phylotastic.org/datelife/index.html#installation)
- [Getting
  started](http://phylotastic.org/datelife/articles/Getting_started_with_datelife.html)
- [Function
  documentation](http://phylotastic.org/datelife/reference/index.html)

## Local installation of the `datelife` R package

`datelife`’s most recent stable version can be installed with:

``` r
install.packages("datelife")
```

`datelife`’s previous stable versions are available for installation
from the CRAN repository. For example, to install `version 0.6.1`, you
can run:

``` r
devtools::install_version("datelife", version="0.6.1")
```

You can install `datelife`’s development version from its GitHub
repository with:

``` r
devtools::install_github("phylotastic/datelife")
```

## Citing `datelife`

If you use `datelife` for a publication, please cite the R package and
the accompanying paper:

<p>
O’Meara B, Sanchez-Reyes L, Eastman J, Heath T, Wright A, Schliep K,
Chamberlain S, Midford P, Harmon L, Brown J, Pennell M, Alfaro M (2023).
<em>datelife: Scientific Data on Time of Lineage Divergence for Your
Taxa</em>. R package version 0.6.7,
<a href="https://doi.org/10.5281/zenodo.593938">https://doi.org/10.5281/zenodo.593938</a>.
</p>
<p>
Sanchez-Reyes L, O’Meara B (2019). “datelife: Leveraging databases and
analytical tools to reveal the dated Tree of Life.” <em>bioRxiv</em>,
<b>782094</b>.
<a href="https://doi.org/10.1101/782094">https://doi.org/10.1101/782094</a>.
</p>

You can get these citations and the bibtex entry with:

``` r
citation("datelife")
toBibtex(citation("datelife"))
```

<!--.bibtex files are available-->

## License

This package is free and open source software, licensed under GPL.

## Acknowledgements

`datelife` has been developed as part of the
[phylotastic](http://phylotastic.org/) (NSF-funded) project, and is
still under development.
