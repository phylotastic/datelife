
<!-- README.md is generated from README.Rmd. Make sure to edit the .Rmd file and not the .md -->

<img src='man/figures/datelife-hexsticker-ai.png' align='right' style='width:150px' />

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

# Welcome to DateLife’s R package GitHub repository!

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
locally. You can find instructions for a local installation below.

If you do not want/have time to deal with installation and R code, you
can use [DateLife’s interactive website
application](https://github.com/phylotastic/datelifedocker). Note that
the website is not live at the moment, apologies.

<!-- http://datelife.opentreeoflife.org/query/ -->
<!--
Get a phylogenetic tree with branch lengths proportional to geologic time (aka a
_**chronogram**_) of any two or more lineages of interest to you.

You can also date a phylogenetic tree of your own making (or choosing one from the literature), using node ages from chronograms found with `datelife` as secondary calibrations.
-->

To learn more, please go to [`datelife`’s documentation
website](http://phylotastic.org/datelife/index.html).

## README topics:

- [Installation](#installation)
- [Citation](#citation)
- [Feedback and info for developers](#feedback)
- [License](#license)

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

## Feedback and Information for Developers

We welcome and encourage to post a [GitHub
issue](https://github.com/phylotastic/datelife/issues/new) with any
comments, ideas and questions about `datelife`’s software and website.
If you want to contribute with code directly, we welcome and encourage
[pull requests](https://github.com/phylotastic/datelife/pulls).

#### *Function documentation*:

Package and function documentation was generated with
[roxygen2](https://CRAN.R-project.org/package=roxygen2):

``` r
roxygen2::roxygenise()
```

#### *Styling code*:

We used the package [lintr](https://CRAN.R-project.org/package=lintr) to
check for coding style:

``` r
lintr::lint_package()
```

#### *Calculating test coverage*:

Code coverage was calculated with the package
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

#### *Generating `datelife`’s hexsticker*:

Code used to generate current `datelife`’s logo hexsticker is in
[data-raw/hexsticker-current.R](https://github.com/phylotastic/datelife/blob/master/data-raw/hexsticker-current.R)

#### *Rendering the vignettes*:

Vignettes are rendered automatically upon built. However, if you wish to
see how they look rendered before releasing the package, you can do this
with `knitr::knit()`. The following command renders the vignette
`Getting_started_with_datelife` as html:

``` r
knitr::knit("vignettes/Getting_started_with_datelife.Rmd")
```

To update “pre-rendered” vignettes, follow [this
blog](https://ropensci.org/blog/2019/12/08/precompute-vignettes/#the-solution-locally-knitting-rmarkdown).
For example, to render the vignette about making trees with BOLD, do:

``` r
knitr::knit("vignettes/making_bold_trees.Rmd.orig", output = "vignettes/making_bold_trees.Rmd")
```

#### *Creating a documentation website for the package*

Using `pkgdown` for this is quite straightforward and fun:

``` r
usethis::use_pkgdown()
pkgdown::build_site()
```

#### *Preparing a CRAN release*

##### ***Updating GitHub actions R CMD check***

Run the following function from the package
[`usethis`](https://usethis.r-lib.org/reference/github_actions.html) to
update R CMD Check on GitHub:

``` r
usethis::use_github_action_check_standard()
```

This downloads the standard R CMD check workflow from [r-lib action
examples](https://github.com/r-lib/actions/blob/v2/examples/check-standard.yaml).

##### ***Local checks***

To be able to release to CRAN, the first step is to pass the checks
locally. To run a local check, you can use the command `R CMD check`
from your terminal. For that, change directories to the one above your
working clone of the `datelife` repo:

``` bash
cd ../
```

Generate a tar ball for your package by running
`R CMD build package-name`:

``` bash
R CMD build datelife
```

Finally, run `R CMD check package-tar-ball` on the tar ball that you
just generated:

``` bash
R CMD check --as-cran datelife_0.6.0.tar.gz
```

##### ***Remote checks***

If you do not have access to different OS to test your package on, the
[rhub](https://CRAN.R-project.org/package=rhub) package allows remote
testing on a variety of OS with the command:

``` r
rhub::check_for_cran()
```

For more `rhub` useful workflows, check out its
[documentation](https://r-hub.github.io/rhub/articles/rhub.html). For
example

``` r
previous_checks <- rhub::list_package_checks(".",
                                  email = "sanchez.reyes.luna@gmail.com",
                                  howmany = 4)
group_check <- rhub::get_check(previous_checks$group[1])
group_check

cran_prep <- check_for_cran()
cran_prep$cran_summary()
```

To [check for URL
validity](https://blog.r-hub.io/2020/12/01/url-checks/) and on Windows
OS, use:

``` r
devtools::check_win_release()
devtools::check_win_devel()
```

#### *Releasing to CRAN*

To submit to CRAN call:

    devtools::release()

and answer the prompted questions. If the answer to all of these is
*yes*, the package will be submitted to CRAN :rocket:

## License

This package is free and open source software, licensed under GPL.

## Acknowledgements

`datelife` has been developed as part of the
[phylotastic](http://phylotastic.org/) (NSF-funded) project, and is
still under development.
