<!--TODO:
  - add congruification step to `datelife_use` functions

DONE:
-->

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
  - varius fixes to `make_datelife_query()` _et al._ functions
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
