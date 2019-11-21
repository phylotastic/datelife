
<!-- README.md is generated from README.Rmd. Please edit THIS file -->

[![Build
Status](https://travis-ci.org/phylotastic/datelife.svg)](https://travis-ci.org/phylotastic/datelife)
[![DOI](https://zenodo.org/badge/23036/phylotastic/datelife.svg)](https://zenodo.org/badge/latestdoi/23036/phylotastic/datelife)
[![codecov](https://codecov.io/gh/phylotastic/datelife/branch/master/graph/badge.svg)](https://codecov.io/gh/phylotastic/datelife)
[![Github Open
Issues](https://img.shields.io/github/issues-raw/phylotastic/datelife.svg)](https://github.com/phylotastic/datelife/issues)
[![Github Closed
Issues](https://img.shields.io/github/issues-closed-raw/phylotastic/datelife.svg)](https://github.com/phylotastic/datelife/issues?q=is%3Aissue+is%3Aclosed)

# `datelife`

Get a phylogenetic tree with branch lengths proportional to geologic
time (aka a ***chronogram***) of any two or more lineages of interest to
you. Use this R package or go to
[www.datelife.org](http://datelife.org/query/) to make a query of
chronograms available for your lineages in the [Open Tree of
Life](http://opentreeoflife.org)’s tree store. A phylogenetic tree of
your own making or choosing can be dated using ages from these
chronograms and the congruification method ([Eastman et
al. 2013](http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12051/abstract))
implemented in *geiger* ([Harmon et
al. 2008](http://bioinformatics.oxfordjournals.org/content/24/1/129.short)).
`datelife` was developed as part of the phylotastic (NSF-funded)
project. It is still under rapid development.

  - [Installation](#installation)
  - [Quick Intro](#quick-intro)
  - [Citing `datelife`](#citation)
  - [Feedback](#feedback)
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
function

## Citation

If you use `datelife`, please cite the R package and the accompanying
publication:

<p>

O’Meara B, Reyes LLS, Eastman J, Heath T, Wright A, Schliep K,
Chamberlain S, Midford P, Harmon L, Brown J, Pennell M, Alfaro M (2019).
<em>datelife: Go from a list of taxa or a tree to a chronogram using
open scientific data</em>. R package version 0.2.19.

</p>

<p>

Sanchez Reyes L, O’Meara B (2019). “datelife: leveraging databases and
analytical tools to reveal the dated tree of life.” <em>Bioarxiv</em>,
<b>xx</b>, xxx-xxx.

</p>

You can get these citations and the bibtex entry with:

``` r
citation("datelife")
toBibtex(citation("datelife"))
```

## Feedback

All comments, ideas and questions about \`datelife1 are encouraged. You
are welcome to post an issue
[here](https://github.com/phylotastic/datelife/issues/new) or to make a
[pull request](https://github.com/phylotastic/datelife/pulls) if you
want to contribute with code directly.

## License

This package is free and open source software, licensed under GPL.
