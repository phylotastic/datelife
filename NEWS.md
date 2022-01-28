# datelife 0.6.1
  - Functions:
    - `get_otol_chronograms` was updated and renamed to `get_opentree_chronograms`
    - Update `match_all_calibrations`
    - Added a `summary` method for `datelifeResult` objects
  - Documentation:
    - Added a "More" section describing return value attributes (will rename to "Attributes")
    - Eliminated unnecessary examples
  - Vignettes: added case study vignette
  - Description: added `BiocManager` package to imports

# datelife 0.6.0

  - Dependencies: Bioconductor packages are used conditionally
  - Package website with `pkgdown`
  - documentation: expanded for all functions
  - examples and function files are written to tempdir()

# datelife 0.5.0

  - functions:
    - `datelife_query_check` is deprecated
    - `use_each_calibration` renamed to `use_calibrations_each`
    - plotting functions have been moved to `datelifeplot` package

# datelife 0.3.1

  - Now muscle or mafft can be used for alignment of BOLD sequences

# datelife 0.2.19

  - MrBayes and TreePL as dating methods
  - Fully ultrametric summary chronograms and datelife() outputs
  - Overlay and lineage through time plots of source chronograms
  - Groves to summarize trees
