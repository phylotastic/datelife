[![Build Status](https://travis-ci.org/phylotastic/datelife.svg)](https://travis-ci.org/phylotastic/datelife) [![DOI](https://zenodo.org/badge/23036/phylotastic/datelife.svg)](https://zenodo.org/badge/latestdoi/23036/phylotastic/datelife) [![codecov](https://codecov.io/github/bomeara/datelife/branch/master/graph/badge.svg)](https://codecov.io/github/bomeara/datelife)

This has the core functionality of DateLife, part of Phylotastic (NSF-funded). It is still under rapid development. With this package, you can get ages of clades or even chronograms using information in [OpenTree](http://opentreeoflife.org)'s tree store. Newick strings or phylo objects can also be dated using these trees and the congruifier, part of Geiger ([Eastman et al. 2013](http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12051/abstract), [Harmon et al. 2008](http://bioinformatics.oxfordjournals.org/content/24/1/129.short)).

For now, requires github version of rotl:

`devtools::install_github("ropensci/rotl")`

Then install datelife:

`devtools::install_github("phylotastic/datelife")`
