<!--TODO:
  - add congruification step to `datelife_use` functions

DONE:
-->
# datelife v0.6.9
- fix bug in `check_ott_input()`
- fix bug in `make_datelife_query()` when getting ott ids
- use spaces instead of "_" to run tnrs, avoids bug from `rotl::tnrs_match_names()` v3.0.14
- faster and more accurate way to get study ids and tree ids from opentree API on `get_opentree_chronograms()` that have branch lengths in Myrs, no relative time.
- update to chronogram database, now with 292 chronograms.
- `get_taxon_summary()` now manages case when `datelife_result` is empty; throws warning instead of criptic error.
- better `testthat` suite for `datelife_search()` inner functions.


# datelife v0.6.8
- fix bug in function `extract_calibrations_phylo()`
- update messages in `calibrations_match()`
- add faster function to retrieve descendants
- updates for new rotl version
- fix uri in DESCRIPTION

# datelife v0.6.7
  - added vignette for bold data workflow
  - `use_calibrations_bladj.matchedCalibrations()`:
    - can use a root age provided by the user
    - if there is no root age, provides a heavy warning and uses the age of the maximum calibration available plus one unit of standard deviation -- if multiple calibrations available, or plus 0.1*calibration age if only one calibration is available.
  - update links in `pkgdown/badges.Rmd` and `pkgdown/presentation.Rmd`

# datelife v0.6.6
  - when congruification fails, return `NA` and produce a warning on functions that use `geiger::congruify.phylo()`:
    - `congruify_and_mrca_phylo()`
    - `extract_calibrations_phylo()`
    - `congruify_and_check()`
  - changes to `congruify_and_mrca_multiPhylo()`:
    - filters out results of `congruify_and_mrca_phylo()` that do not inherit `data.frame`; usually this happens when congruification fails
    - correctly assigns congruified `phy` as attribute
  - `use_calibrations_bladj.matchedCalibrations()` correctly uses `calibrations` argument
  - `get_ott_children()` has error handling for when OpenTree APIs might be down, or there is no internet connection
  - various fixes to `make_datelife_query()` _et al._ functions
  - added a new vignette for making datelife query functions (allows testing various functionalities)

# datelife v0.6.5
  - data:
    - updated `opentree_chronograms` object. It now has 253 chronograms from Open Tree of Life and uses "xz" compression.
  - examples:
    - set to not test examples
  - functions:
    - added functions `matrix_to_table` and `matrices_to_table` that go from a matrix of patristic distances to a table of taxon name pairs and respective node ages.
    - bug fix on `use_calibrations_bladj` that used element $present_calibrations instead of $matched_calibrations
    - `get_otol_chronogram` is another name to call `get_opentree_chronograms`
    - added taxonomic source options argument to `make_datelife_query`. You can choose from OTT, NCBI, IRMNG and GBIF.
    - function `congruify_and_mrca`: output has congruified topology with nodelables as attribute.
    - `use_calibrations_bladj` takes an output of `congruify_and_mrca` functions.
  - added function method `congruify_and_mrca`

# datelife v0.6.1
  - Functions:
    - `get_otol_chronograms` was updated and renamed to `get_opentree_chronograms`
    - Update `match_all_calibrations`
    - Added a `summary` method for `datelifeResult` objects
  - Documentation:
    - Added a "More" section describing return value attributes (will rename to "Attributes")
    - Eliminated unnecessary examples
  - Vignettes: added case study vignette
  - Description: added `BiocManager` package to imports

# datelife v0.6.0

  - Dependencies: Bioconductor packages are used conditionally
  - Package website with `pkgdown`
  - documentation: expanded for all functions
  - examples and function files are written to tempdir()
  - function rename: `get_biggest_phylo` to `get_biggest_multiphylo`

# datelife v0.5.0

  - functions:
    - `datelife_query_check` is deprecated
    - `use_each_calibration` renamed to `use_calibrations_each`
    - plotting functions have been moved to `datelifeplot` package

# datelife v0.3.1

  - Now muscle or mafft can be used for alignment of BOLD sequences

# datelife v0.2.19

  - MrBayes and TreePL as dating methods
  - Fully ultrametric summary chronograms and datelife() outputs
  - Overlay and lineage through time plots of source chronograms
  - Groves to summarize trees
