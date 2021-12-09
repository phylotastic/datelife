---
title: "Getting started with `datelife`"
output: 
  rmarkdown::html_vignette:
    keep_md: true
vignette: >
  %\VignetteIndexEntry{Getting started with `datelife`}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---




```r
library(datelife)
```

## What can you do with `datelife`?

The aim of the `datelife` R package is to allow researchers and the general audience to obtain data on ages for their organisms of interest. This age data is extracted from dated phylogenetic trees (chronograms) that have been published and peer-reviewed in association with a scientific article in an indexed journal. As such, the age data provided by `datelife` constitutes state-of-the-art, peer-reviewed, public scientific knowledge, that can be accessed by experts and non-experts in the field alike. 

Users can get age data for their organisms of interest in several forms:

- **MRCA**: They can get the age of the most recent common ancestor (MRCA) of their organisms.

- **Original chronogram**: They can get the chronograms from the original publications, pruned to contain their organisms only.

- **Biggest chronogram**: They can get the original chronogram with the most organisms of interest present as tips.

- **Summary chronogram**: They can obtain a single chronogram, summarized using data from the original chronograms. The summary chronogram can be calculated using different methods.

Users can also apply the obtained age data to date a phylogenetic tree topology of their organisms of interest.

## 1. Preparing your organism names

The function `make_datelife_query` processes the scientific names of your organisms of interest and saves them as a `datelifeQuery` object.
This allows you to easily use them later to perform a search for ages across the `datelife` chronogram database. You have two general ways of providing organism names:

### a) From a set of organism names
You can provide the scientific names as a character vector or as a single comma separated string of names.
For example, running the following:


```r
datelife::make_datelife_query(input = c("Delphinus_delphis", 
                                        "Gallus gallus", 
                                        "elephas Maximus",
                                        "felis_catus"))
#> ... Making a DateLife query.
#> ... Phylo-processing 'input'.
#> 'input' is not a phylogeny.
#> DateLife query done!
#> $cleaned_names
#> [1] "Delphinus delphis" "Gallus gallus"     "elephas Maximus"  
#> [4] "felis catus"      
#> 
#> $ott_ids
#> NULL
#> 
#> $phy
#> [1] NA
#> 
#> attr(,"class")
#> [1] "datelifeQuery"
```

Gives the same result as running:


```r
datelife::make_datelife_query(input = "Delphinus_delphis, 
                                       Gallus gallus, 
                                       elephas Maximus, 
                                       felis_catus")
#> ... Making a DateLife query.
#> ... Phylo-processing 'input'.
#> 'input' is not a phylogeny.
#> DateLife query done!
#> $cleaned_names
#> [1] "Delphinus delphis" "Gallus gallus"     "elephas Maximus"  
#> [4] "felis catus"      
#> 
#> $ott_ids
#> NULL
#> 
#> $phy
#> [1] NA
#> 
#> attr(,"class")
#> [1] "datelifeQuery"
```

Note that you can use underscores or spaces to separate binomial scientific names and `datelife` will process them.

### b) From a tree with organism names on its tips

You can also provide the scientific names of your organisms of interest as labels on the tips of any phylogenetic tree.
You can provide the tree as a newick string:


```r
newick <- "((Elephas_maximus,(Delphinus_delphis, Felis_silvestris)), Gallus_gallus);"
datelife::make_datelife_query(input = newick)
#> ... Making a DateLife query.
#> ... Phylo-processing 'input'.
#> 'input' is a phylogeny and it is correcly formatted.
#> DateLife query done!
#> $cleaned_names
#> [1] "Elephas maximus"   "Delphinus delphis" " Felis silvestris"
#> [4] " Gallus gallus"   
#> 
#> $ott_ids
#> NULL
#> 
#> $phy
#> 
#> Phylogenetic tree with 4 tips and 3 internal nodes.
#> 
#> Tip labels:
#>   Elephas_maximus, Delphinus_delphis, _Felis_silvestris, _Gallus_gallus
#> 
#> Rooted; no branch lengths.
#> 
#> attr(,"class")
#> [1] "datelifeQuery"
```

You can also provide the tree as a `phylo` object. In this example, you have to first convert your newick string to `phylo`. You can do this with the function `read.tree` from the `ape` package. If you already have a tree as a `phylo` object, you can use that directly:


```r
phylo <- ape::read.tree(text = newick)
datelife::make_datelife_query(input = phylo)
#> ... Making a DateLife query.
#> ... Phylo-processing 'input'.
#> 'input' is a phylogeny and it is correcly formatted.
#> DateLife query done!
#> $cleaned_names
#> [1] "Elephas maximus"   "Delphinus delphis" "Felis silvestris" 
#> [4] "Gallus gallus"    
#> 
#> $ott_ids
#> NULL
#> 
#> $phy
#> 
#> Phylogenetic tree with 4 tips and 3 internal nodes.
#> 
#> Tip labels:
#>   Elephas_maximus, Delphinus_delphis, Felis_silvestris, Gallus_gallus
#> 
#> Rooted; no branch lengths.
#> 
#> attr(,"class")
#> [1] "datelifeQuery"
```

Now save your `datelifeQuery` object for the next section:


```r
my_datelife_query <- datelife::make_datelife_query(input = phylo)
#> ... Making a DateLife query.
#> ... Phylo-processing 'input'.
#> 'input' is a phylogeny and it is correcly formatted.
#> DateLife query done!
```

## 2. Searching the chronogram database

The `datelife_search` function takes a `datelifeQuery` object and goes through the `datelife` chronogram database to find all chronograms matching at least two organism names. 


```r
datelife::datelife_search(input = my_datelife_query)
#> ... Running a DateLife search.
#> ... Getting a DateLife result.
#> DateLife result obtained!
#> Source chronograms from:
#> 1: Bininda-Emonds, Olaf R. P., Marcel Cardillo, Kate E. Jones, Ross D. E. MacPhee, Robin M. D. Beck, Richard Grenyer, Samantha A. Price, Rutger A. Vos, John L. Gittleman, Andy Purvis. 2007. The delayed rise of present-day mammals. Nature 446 (7135): 507-512
#> 2: Bininda-Emonds, Olaf R. P., Marcel Cardillo, Kate E. Jones, Ross D. E. MacPhee, Robin M. D. Beck, Richard Grenyer, Samantha A. Price, Rutger A. Vos, John L. Gittleman, Andy Purvis. 2007. The delayed rise of present-day mammals. Nature 446 (7135): 507-512
#> 3: Bininda-Emonds, Olaf R. P., Marcel Cardillo, Kate E. Jones, Ross D. E. MacPhee, Robin M. D. Beck, Richard Grenyer, Samantha A. Price, Rutger A. Vos, John L. Gittleman, Andy Purvis. 2007. The delayed rise of present-day mammals. Nature 446 (7135): 507-512
#> 4: Hedges, S. Blair, Julie Marin, Michael Suleski, Madeline Paymer, Sudhir Kumar. 2015. Tree of life reveals clock-like speciation and diversification. Molecular Biology and Evolution 32 (4): 835-845
#> Input taxa presence across source chronograms:
#>               taxon chronograms
#> 1   Elephas maximus         4/4
#> 2 Delphinus delphis         4/4
#> 3  Felis silvestris         4/4
#> 4     Gallus gallus         1/4
#> Input taxa completely absent from source chronograms:
#>   taxon
#> 1  None
#> DateLife search done!
#> 4 phylogenetic trees
```

`datelife_search` will also summarize the results from the chronogram database search for you.
By default, it will give you all original chronograms matching the search, as a `phylo` or `multiPhylo` object. You can specify the format you desire for the summary of results with the argument `summary_format`.
There are different summary formats to choose from:

`summary_format = "citations"` gives back the study citations from all  original chronograms matching your search:


```r
datelife::datelife_search(input = my_datelife_query,
                          summary_format = "citations")
#> ... Running a DateLife search.
#> ... Getting a DateLife result.
#> DateLife result obtained!
#> Input taxa presence across source chronograms:
#>               taxon chronograms
#> 1   Elephas maximus         4/4
#> 2 Delphinus delphis         4/4
#> 3  Felis silvestris         4/4
#> 4     Gallus gallus         1/4
#> Input taxa completely absent from source chronograms:
#>   taxon
#> 1  None
#> DateLife search done!
#> [1] "Bininda-Emonds, Olaf R. P., Marcel Cardillo, Kate E. Jones, Ross D. E. MacPhee, Robin M. D. Beck, Richard Grenyer, Samantha A. Price, Rutger A. Vos, John L. Gittleman, Andy Purvis. 2007. The delayed rise of present-day mammals. Nature 446 (7135): 507-512"
#> [2] "Bininda-Emonds, Olaf R. P., Marcel Cardillo, Kate E. Jones, Ross D. E. MacPhee, Robin M. D. Beck, Richard Grenyer, Samantha A. Price, Rutger A. Vos, John L. Gittleman, Andy Purvis. 2007. The delayed rise of present-day mammals. Nature 446 (7135): 507-512"
#> [3] "Bininda-Emonds, Olaf R. P., Marcel Cardillo, Kate E. Jones, Ross D. E. MacPhee, Robin M. D. Beck, Richard Grenyer, Samantha A. Price, Rutger A. Vos, John L. Gittleman, Andy Purvis. 2007. The delayed rise of present-day mammals. Nature 446 (7135): 507-512"
#> [4] "Hedges, S. Blair, Julie Marin, Michael Suleski, Madeline Paymer, Sudhir Kumar. 2015. Tree of life reveals clock-like speciation and diversification. Molecular Biology and Evolution 32 (4): 835-845"                                                          
#> attr(,"datelife_result")
#> $`Bininda-Emonds, Olaf R. P., Marcel Cardillo, Kate E. Jones, Ross D. E. MacPhee, Robin M. D. Beck, Richard Grenyer, Samantha A. Price, Rutger A. Vos, John L. Gittleman, Andy Purvis. 2007. The delayed rise of present-day mammals. Nature 446 (7135): 507-512`
#>                   Elephas maximus Delphinus delphis
#> Elephas maximus               0.0             202.6
#> Delphinus delphis           202.6               0.0
#> Felis silvestris            202.6             177.0
#>                   Felis silvestris
#> Elephas maximus              202.6
#> Delphinus delphis            177.0
#> Felis silvestris               0.0
#> 
#> $`Bininda-Emonds, Olaf R. P., Marcel Cardillo, Kate E. Jones, Ross D. E. MacPhee, Robin M. D. Beck, Richard Grenyer, Samantha A. Price, Rutger A. Vos, John L. Gittleman, Andy Purvis. 2007. The delayed rise of present-day mammals. Nature 446 (7135): 507-512`
#>                   Elephas maximus Felis silvestris
#> Elephas maximus               0.0            192.2
#> Felis silvestris            192.2              0.0
#> Delphinus delphis           192.2            173.6
#>                   Delphinus delphis
#> Elephas maximus               192.2
#> Felis silvestris              173.6
#> Delphinus delphis               0.0
#> 
#> $`Bininda-Emonds, Olaf R. P., Marcel Cardillo, Kate E. Jones, Ross D. E. MacPhee, Robin M. D. Beck, Richard Grenyer, Samantha A. Price, Rutger A. Vos, John L. Gittleman, Andy Purvis. 2007. The delayed rise of present-day mammals. Nature 446 (7135): 507-512`
#>                   Delphinus delphis Felis silvestris
#> Delphinus delphis               0.0            179.6
#> Felis silvestris              179.6              0.0
#> Elephas maximus               221.6            221.6
#>                   Elephas maximus
#> Delphinus delphis           221.6
#> Felis silvestris            221.6
#> Elephas maximus               0.0
#> 
#> $`Hedges, S. Blair, Julie Marin, Michael Suleski, Madeline Paymer, Sudhir Kumar. 2015. Tree of life reveals clock-like speciation and diversification. Molecular Biology and Evolution 32 (4): 835-845`
#>                   Felis silvestris Delphinus delphis
#> Felis silvestris            0.0000          156.7850
#> Delphinus delphis         156.7850            0.0000
#> Elephas maximus           210.0193          210.0193
#> Gallus gallus             641.0476          641.0476
#>                   Elephas maximus Gallus gallus
#> Felis silvestris         210.0193      641.0476
#> Delphinus delphis        210.0193      641.0476
#> Elephas maximus            0.0000      641.0476
#> Gallus gallus            641.0476        0.0000
#> 
#> attr(,"class")
#> [1] "datelifeResult"
#> attr(,"datelife_query")
#> $cleaned_names
#> [1] "Elephas maximus"   "Delphinus delphis" "Felis silvestris" 
#> [4] "Gallus gallus"    
#> 
#> $ott_ids
#> NULL
#> 
#> $phy
#> 
#> Phylogenetic tree with 4 tips and 3 internal nodes.
#> 
#> Tip labels:
#>   Elephas_maximus, Delphinus_delphis, Felis_silvestris, Gallus_gallus
#> 
#> Rooted; no branch lengths.
#> 
#> attr(,"class")
#> [1] "datelifeQuery"
```

`summary_format = "mrca"` gives back the age of the MRCA of your organisms of interest:


```r
datelife::datelife_search(input = my_datelife_query,
                          summary_format = "mrca")
#> ... Running a DateLife search.
#> ... Getting a DateLife result.
#> DateLife result obtained!
#> Source chronograms from:
#> 1: Bininda-Emonds, Olaf R. P., Marcel Cardillo, Kate E. Jones, Ross D. E. MacPhee, Robin M. D. Beck, Richard Grenyer, Samantha A. Price, Rutger A. Vos, John L. Gittleman, Andy Purvis. 2007. The delayed rise of present-day mammals. Nature 446 (7135): 507-512
#> 2: Bininda-Emonds, Olaf R. P., Marcel Cardillo, Kate E. Jones, Ross D. E. MacPhee, Robin M. D. Beck, Richard Grenyer, Samantha A. Price, Rutger A. Vos, John L. Gittleman, Andy Purvis. 2007. The delayed rise of present-day mammals. Nature 446 (7135): 507-512
#> 3: Bininda-Emonds, Olaf R. P., Marcel Cardillo, Kate E. Jones, Ross D. E. MacPhee, Robin M. D. Beck, Richard Grenyer, Samantha A. Price, Rutger A. Vos, John L. Gittleman, Andy Purvis. 2007. The delayed rise of present-day mammals. Nature 446 (7135): 507-512
#> 4: Hedges, S. Blair, Julie Marin, Michael Suleski, Madeline Paymer, Sudhir Kumar. 2015. Tree of life reveals clock-like speciation and diversification. Molecular Biology and Evolution 32 (4): 835-845
#> Input taxa presence across source chronograms:
#>               taxon chronograms
#> 1   Elephas maximus         4/4
#> 2 Delphinus delphis         4/4
#> 3  Felis silvestris         4/4
#> 4     Gallus gallus         1/4
#> Input taxa completely absent from source chronograms:
#>   taxon
#> 1  None
#> DateLife search done!
#> Bininda-Emonds, Olaf R. P., Marcel Cardillo, Kate E. Jones, Ross D. E. MacPhee, Robin M. D. Beck, Richard Grenyer, Samantha A. Price, Rutger A. Vos, John L. Gittleman, Andy Purvis. 2007. The delayed rise of present-day mammals. Nature 446 (7135): 507-512 
#>                                                                                                                                                                                                                                                       101.3000 
#> Bininda-Emonds, Olaf R. P., Marcel Cardillo, Kate E. Jones, Ross D. E. MacPhee, Robin M. D. Beck, Richard Grenyer, Samantha A. Price, Rutger A. Vos, John L. Gittleman, Andy Purvis. 2007. The delayed rise of present-day mammals. Nature 446 (7135): 507-512 
#>                                                                                                                                                                                                                                                        96.1000 
#> Bininda-Emonds, Olaf R. P., Marcel Cardillo, Kate E. Jones, Ross D. E. MacPhee, Robin M. D. Beck, Richard Grenyer, Samantha A. Price, Rutger A. Vos, John L. Gittleman, Andy Purvis. 2007. The delayed rise of present-day mammals. Nature 446 (7135): 507-512 
#>                                                                                                                                                                                                                                                       110.8000 
#>                                                           Hedges, S. Blair, Julie Marin, Michael Suleski, Madeline Paymer, Sudhir Kumar. 2015. Tree of life reveals clock-like speciation and diversification. Molecular Biology and Evolution 32 (4): 835-845 
#>                                                                                                                                                                                                                                                       320.5238 
#> attr(,"datelife_result")
#> $`Bininda-Emonds, Olaf R. P., Marcel Cardillo, Kate E. Jones, Ross D. E. MacPhee, Robin M. D. Beck, Richard Grenyer, Samantha A. Price, Rutger A. Vos, John L. Gittleman, Andy Purvis. 2007. The delayed rise of present-day mammals. Nature 446 (7135): 507-512`
#>                   Elephas maximus Delphinus delphis
#> Elephas maximus               0.0             202.6
#> Delphinus delphis           202.6               0.0
#> Felis silvestris            202.6             177.0
#>                   Felis silvestris
#> Elephas maximus              202.6
#> Delphinus delphis            177.0
#> Felis silvestris               0.0
#> 
#> $`Bininda-Emonds, Olaf R. P., Marcel Cardillo, Kate E. Jones, Ross D. E. MacPhee, Robin M. D. Beck, Richard Grenyer, Samantha A. Price, Rutger A. Vos, John L. Gittleman, Andy Purvis. 2007. The delayed rise of present-day mammals. Nature 446 (7135): 507-512`
#>                   Elephas maximus Felis silvestris
#> Elephas maximus               0.0            192.2
#> Felis silvestris            192.2              0.0
#> Delphinus delphis           192.2            173.6
#>                   Delphinus delphis
#> Elephas maximus               192.2
#> Felis silvestris              173.6
#> Delphinus delphis               0.0
#> 
#> $`Bininda-Emonds, Olaf R. P., Marcel Cardillo, Kate E. Jones, Ross D. E. MacPhee, Robin M. D. Beck, Richard Grenyer, Samantha A. Price, Rutger A. Vos, John L. Gittleman, Andy Purvis. 2007. The delayed rise of present-day mammals. Nature 446 (7135): 507-512`
#>                   Delphinus delphis Felis silvestris
#> Delphinus delphis               0.0            179.6
#> Felis silvestris              179.6              0.0
#> Elephas maximus               221.6            221.6
#>                   Elephas maximus
#> Delphinus delphis           221.6
#> Felis silvestris            221.6
#> Elephas maximus               0.0
#> 
#> $`Hedges, S. Blair, Julie Marin, Michael Suleski, Madeline Paymer, Sudhir Kumar. 2015. Tree of life reveals clock-like speciation and diversification. Molecular Biology and Evolution 32 (4): 835-845`
#>                   Felis silvestris Delphinus delphis
#> Felis silvestris            0.0000          156.7850
#> Delphinus delphis         156.7850            0.0000
#> Elephas maximus           210.0193          210.0193
#> Gallus gallus             641.0476          641.0476
#>                   Elephas maximus Gallus gallus
#> Felis silvestris         210.0193      641.0476
#> Delphinus delphis        210.0193      641.0476
#> Elephas maximus            0.0000      641.0476
#> Gallus gallus            641.0476        0.0000
#> 
#> attr(,"class")
#> [1] "datelifeResult"
#> attr(,"datelife_query")
#> $cleaned_names
#> [1] "Elephas maximus"   "Delphinus delphis" "Felis silvestris" 
#> [4] "Gallus gallus"    
#> 
#> $ott_ids
#> NULL
#> 
#> $phy
#> 
#> Phylogenetic tree with 4 tips and 3 internal nodes.
#> 
#> Tip labels:
#>   Elephas_maximus, Delphinus_delphis, Felis_silvestris, Gallus_gallus
#> 
#> Rooted; no branch lengths.
#> 
#> attr(,"class")
#> [1] "datelifeQuery"
```

`summary_format = "newick_all"`


```r
datelife::datelife_search(input = my_datelife_query,
                          summary_format = "newick_all")
#> ... Running a DateLife search.
#> ... Getting a DateLife result.
#> DateLife result obtained!
#> Source chronograms from:
#> 1: Bininda-Emonds, Olaf R. P., Marcel Cardillo, Kate E. Jones, Ross D. E. MacPhee, Robin M. D. Beck, Richard Grenyer, Samantha A. Price, Rutger A. Vos, John L. Gittleman, Andy Purvis. 2007. The delayed rise of present-day mammals. Nature 446 (7135): 507-512
#> 2: Bininda-Emonds, Olaf R. P., Marcel Cardillo, Kate E. Jones, Ross D. E. MacPhee, Robin M. D. Beck, Richard Grenyer, Samantha A. Price, Rutger A. Vos, John L. Gittleman, Andy Purvis. 2007. The delayed rise of present-day mammals. Nature 446 (7135): 507-512
#> 3: Bininda-Emonds, Olaf R. P., Marcel Cardillo, Kate E. Jones, Ross D. E. MacPhee, Robin M. D. Beck, Richard Grenyer, Samantha A. Price, Rutger A. Vos, John L. Gittleman, Andy Purvis. 2007. The delayed rise of present-day mammals. Nature 446 (7135): 507-512
#> 4: Hedges, S. Blair, Julie Marin, Michael Suleski, Madeline Paymer, Sudhir Kumar. 2015. Tree of life reveals clock-like speciation and diversification. Molecular Biology and Evolution 32 (4): 835-845
#> Input taxa presence across source chronograms:
#>               taxon chronograms
#> 1   Elephas maximus         4/4
#> 2 Delphinus delphis         4/4
#> 3  Felis silvestris         4/4
#> 4     Gallus gallus         1/4
#> Input taxa completely absent from source chronograms:
#>   taxon
#> 1  None
#> DateLife search done!
#> Bininda-Emonds, Olaf R. P., Marcel Cardillo, Kate E. Jones, Ross D. E. MacPhee, Robin M. D. Beck, Richard Grenyer, Samantha A. Price, Rutger A. Vos, John L. Gittleman, Andy Purvis. 2007. The delayed rise of present-day mammals. Nature 446 (7135): 507-512 
#>                                                                                                                                                                                 "((Felis_silvestris:88.5,Delphinus_delphis:88.5):12.8,Elephas_maximus:101.3);" 
#> Bininda-Emonds, Olaf R. P., Marcel Cardillo, Kate E. Jones, Ross D. E. MacPhee, Robin M. D. Beck, Richard Grenyer, Samantha A. Price, Rutger A. Vos, John L. Gittleman, Andy Purvis. 2007. The delayed rise of present-day mammals. Nature 446 (7135): 507-512 
#>                                                                                                                                                                                   "((Delphinus_delphis:86.8,Felis_silvestris:86.8):9.3,Elephas_maximus:96.1);" 
#> Bininda-Emonds, Olaf R. P., Marcel Cardillo, Kate E. Jones, Ross D. E. MacPhee, Robin M. D. Beck, Richard Grenyer, Samantha A. Price, Rutger A. Vos, John L. Gittleman, Andy Purvis. 2007. The delayed rise of present-day mammals. Nature 446 (7135): 507-512 
#>                                                                                                                                                                                   "((Felis_silvestris:89.8,Delphinus_delphis:89.8):21,Elephas_maximus:110.8);" 
#>                                                           Hedges, S. Blair, Julie Marin, Michael Suleski, Madeline Paymer, Sudhir Kumar. 2015. Tree of life reveals clock-like speciation and diversification. Molecular Biology and Evolution 32 (4): 835-845 
#>                                                                                                                     "((Elephas_maximus:105.009658,(Felis_silvestris:78.392518,Delphinus_delphis:78.392517):26.617141):215.5141645,Gallus_gallus:320.5238235);" 
#> attr(,"datelife_result")
#> $`Bininda-Emonds, Olaf R. P., Marcel Cardillo, Kate E. Jones, Ross D. E. MacPhee, Robin M. D. Beck, Richard Grenyer, Samantha A. Price, Rutger A. Vos, John L. Gittleman, Andy Purvis. 2007. The delayed rise of present-day mammals. Nature 446 (7135): 507-512`
#>                   Elephas maximus Delphinus delphis
#> Elephas maximus               0.0             202.6
#> Delphinus delphis           202.6               0.0
#> Felis silvestris            202.6             177.0
#>                   Felis silvestris
#> Elephas maximus              202.6
#> Delphinus delphis            177.0
#> Felis silvestris               0.0
#> 
#> $`Bininda-Emonds, Olaf R. P., Marcel Cardillo, Kate E. Jones, Ross D. E. MacPhee, Robin M. D. Beck, Richard Grenyer, Samantha A. Price, Rutger A. Vos, John L. Gittleman, Andy Purvis. 2007. The delayed rise of present-day mammals. Nature 446 (7135): 507-512`
#>                   Elephas maximus Felis silvestris
#> Elephas maximus               0.0            192.2
#> Felis silvestris            192.2              0.0
#> Delphinus delphis           192.2            173.6
#>                   Delphinus delphis
#> Elephas maximus               192.2
#> Felis silvestris              173.6
#> Delphinus delphis               0.0
#> 
#> $`Bininda-Emonds, Olaf R. P., Marcel Cardillo, Kate E. Jones, Ross D. E. MacPhee, Robin M. D. Beck, Richard Grenyer, Samantha A. Price, Rutger A. Vos, John L. Gittleman, Andy Purvis. 2007. The delayed rise of present-day mammals. Nature 446 (7135): 507-512`
#>                   Delphinus delphis Felis silvestris
#> Delphinus delphis               0.0            179.6
#> Felis silvestris              179.6              0.0
#> Elephas maximus               221.6            221.6
#>                   Elephas maximus
#> Delphinus delphis           221.6
#> Felis silvestris            221.6
#> Elephas maximus               0.0
#> 
#> $`Hedges, S. Blair, Julie Marin, Michael Suleski, Madeline Paymer, Sudhir Kumar. 2015. Tree of life reveals clock-like speciation and diversification. Molecular Biology and Evolution 32 (4): 835-845`
#>                   Felis silvestris Delphinus delphis
#> Felis silvestris            0.0000          156.7850
#> Delphinus delphis         156.7850            0.0000
#> Elephas maximus           210.0193          210.0193
#> Gallus gallus             641.0476          641.0476
#>                   Elephas maximus Gallus gallus
#> Felis silvestris         210.0193      641.0476
#> Delphinus delphis        210.0193      641.0476
#> Elephas maximus            0.0000      641.0476
#> Gallus gallus            641.0476        0.0000
#> 
#> attr(,"class")
#> [1] "datelifeResult"
#> attr(,"datelife_query")
#> $cleaned_names
#> [1] "Elephas maximus"   "Delphinus delphis" "Felis silvestris" 
#> [4] "Gallus gallus"    
#> 
#> $ott_ids
#> NULL
#> 
#> $phy
#> 
#> Phylogenetic tree with 4 tips and 3 internal nodes.
#> 
#> Tip labels:
#>   Elephas_maximus, Delphinus_delphis, Felis_silvestris, Gallus_gallus
#> 
#> Rooted; no branch lengths.
#> 
#> attr(,"class")
#> [1] "datelifeQuery"
```

`summary_format = "newick_sdm"`


```r
datelife::datelife_search(input = my_datelife_query,
                          summary_format = "newick_sdm")
```

`summary_format = "newick_median"`


```r
datelife::datelife_search(input = my_datelife_query,
                          summary_format = "newick_median")
#> ... Running a DateLife search.
#> ... Getting a DateLife result.
#> DateLife result obtained!
#> Trying with overlap = 2
#> Registered S3 method overwritten by 'compare':
#>   method           from    
#>   print.comparison testthat
#> Success!
#> ... Calculating a median summary chronogram.
#> Source chronograms from:
#> 1: Bininda-Emonds, Olaf R. P., Marcel Cardillo, Kate E. Jones, Ross D. E. MacPhee, Robin M. D. Beck, Richard Grenyer, Samantha A. Price, Rutger A. Vos, John L. Gittleman, Andy Purvis. 2007. The delayed rise of present-day mammals. Nature 446 (7135): 507-512
#> 2: Bininda-Emonds, Olaf R. P., Marcel Cardillo, Kate E. Jones, Ross D. E. MacPhee, Robin M. D. Beck, Richard Grenyer, Samantha A. Price, Rutger A. Vos, John L. Gittleman, Andy Purvis. 2007. The delayed rise of present-day mammals. Nature 446 (7135): 507-512
#> 3: Bininda-Emonds, Olaf R. P., Marcel Cardillo, Kate E. Jones, Ross D. E. MacPhee, Robin M. D. Beck, Richard Grenyer, Samantha A. Price, Rutger A. Vos, John L. Gittleman, Andy Purvis. 2007. The delayed rise of present-day mammals. Nature 446 (7135): 507-512
#> 4: Hedges, S. Blair, Julie Marin, Michael Suleski, Madeline Paymer, Sudhir Kumar. 2015. Tree of life reveals clock-like speciation and diversification. Molecular Biology and Evolution 32 (4): 835-845
#> Input taxa presence across source chronograms:
#>               taxon chronograms
#> 1   Elephas maximus         4/4
#> 2 Delphinus delphis         4/4
#> 3  Felis silvestris         4/4
#> 4     Gallus gallus         1/4
#> Input taxa completely absent from source chronograms:
#>   taxon
#> 1  None
#> DateLife search done!
#> [1] "((Elephas_maximus:103.154831,(Delphinus_delphis:87.650002,Felis_silvestris:87.650002)n3:15.504829)n2:217.369003,Gallus_gallus:320.523834)n1:1;"
#> attr(,"datelife_result")
#> $`Bininda-Emonds, Olaf R. P., Marcel Cardillo, Kate E. Jones, Ross D. E. MacPhee, Robin M. D. Beck, Richard Grenyer, Samantha A. Price, Rutger A. Vos, John L. Gittleman, Andy Purvis. 2007. The delayed rise of present-day mammals. Nature 446 (7135): 507-512`
#>                   Elephas maximus Delphinus delphis
#> Elephas maximus               0.0             202.6
#> Delphinus delphis           202.6               0.0
#> Felis silvestris            202.6             177.0
#>                   Felis silvestris
#> Elephas maximus              202.6
#> Delphinus delphis            177.0
#> Felis silvestris               0.0
#> 
#> $`Bininda-Emonds, Olaf R. P., Marcel Cardillo, Kate E. Jones, Ross D. E. MacPhee, Robin M. D. Beck, Richard Grenyer, Samantha A. Price, Rutger A. Vos, John L. Gittleman, Andy Purvis. 2007. The delayed rise of present-day mammals. Nature 446 (7135): 507-512`
#>                   Elephas maximus Felis silvestris
#> Elephas maximus               0.0            192.2
#> Felis silvestris            192.2              0.0
#> Delphinus delphis           192.2            173.6
#>                   Delphinus delphis
#> Elephas maximus               192.2
#> Felis silvestris              173.6
#> Delphinus delphis               0.0
#> 
#> $`Bininda-Emonds, Olaf R. P., Marcel Cardillo, Kate E. Jones, Ross D. E. MacPhee, Robin M. D. Beck, Richard Grenyer, Samantha A. Price, Rutger A. Vos, John L. Gittleman, Andy Purvis. 2007. The delayed rise of present-day mammals. Nature 446 (7135): 507-512`
#>                   Delphinus delphis Felis silvestris
#> Delphinus delphis               0.0            179.6
#> Felis silvestris              179.6              0.0
#> Elephas maximus               221.6            221.6
#>                   Elephas maximus
#> Delphinus delphis           221.6
#> Felis silvestris            221.6
#> Elephas maximus               0.0
#> 
#> $`Hedges, S. Blair, Julie Marin, Michael Suleski, Madeline Paymer, Sudhir Kumar. 2015. Tree of life reveals clock-like speciation and diversification. Molecular Biology and Evolution 32 (4): 835-845`
#>                   Felis silvestris Delphinus delphis
#> Felis silvestris            0.0000          156.7850
#> Delphinus delphis         156.7850            0.0000
#> Elephas maximus           210.0193          210.0193
#> Gallus gallus             641.0476          641.0476
#>                   Elephas maximus Gallus gallus
#> Felis silvestris         210.0193      641.0476
#> Delphinus delphis        210.0193      641.0476
#> Elephas maximus            0.0000      641.0476
#> Gallus gallus            641.0476        0.0000
#> 
#> attr(,"class")
#> [1] "datelifeResult"
#> attr(,"datelife_query")
#> $cleaned_names
#> [1] "Elephas maximus"   "Delphinus delphis" "Felis silvestris" 
#> [4] "Gallus gallus"    
#> 
#> $ott_ids
#> NULL
#> 
#> $phy
#> 
#> Phylogenetic tree with 4 tips and 3 internal nodes.
#> 
#> Tip labels:
#>   Elephas_maximus, Delphinus_delphis, Felis_silvestris, Gallus_gallus
#> 
#> Rooted; no branch lengths.
#> 
#> attr(,"class")
#> [1] "datelifeQuery"
```

`summary_format = "phylo_sdm"`


```r
datelife::datelife_search(input = my_datelife_query,
                          summary_format = "phylo_sdm")
```

`summary_format = "phylo_median"`


```r
datelife::datelife_search(input = my_datelife_query,
                          summary_format = "phylo_median")
#> ... Running a DateLife search.
#> ... Getting a DateLife result.
#> DateLife result obtained!
#> Trying with overlap = 2
#> Success!
#> ... Calculating a median summary chronogram.
#> Source chronograms from:
#> 1: Bininda-Emonds, Olaf R. P., Marcel Cardillo, Kate E. Jones, Ross D. E. MacPhee, Robin M. D. Beck, Richard Grenyer, Samantha A. Price, Rutger A. Vos, John L. Gittleman, Andy Purvis. 2007. The delayed rise of present-day mammals. Nature 446 (7135): 507-512
#> 2: Bininda-Emonds, Olaf R. P., Marcel Cardillo, Kate E. Jones, Ross D. E. MacPhee, Robin M. D. Beck, Richard Grenyer, Samantha A. Price, Rutger A. Vos, John L. Gittleman, Andy Purvis. 2007. The delayed rise of present-day mammals. Nature 446 (7135): 507-512
#> 3: Bininda-Emonds, Olaf R. P., Marcel Cardillo, Kate E. Jones, Ross D. E. MacPhee, Robin M. D. Beck, Richard Grenyer, Samantha A. Price, Rutger A. Vos, John L. Gittleman, Andy Purvis. 2007. The delayed rise of present-day mammals. Nature 446 (7135): 507-512
#> 4: Hedges, S. Blair, Julie Marin, Michael Suleski, Madeline Paymer, Sudhir Kumar. 2015. Tree of life reveals clock-like speciation and diversification. Molecular Biology and Evolution 32 (4): 835-845
#> Input taxa presence across source chronograms:
#>               taxon chronograms
#> 1   Elephas maximus         4/4
#> 2 Delphinus delphis         4/4
#> 3  Felis silvestris         4/4
#> 4     Gallus gallus         1/4
#> Input taxa completely absent from source chronograms:
#>   taxon
#> 1  None
#> DateLife search done!
#> 
#> Phylogenetic tree with 4 tips and 3 internal nodes.
#> 
#> Tip labels:
#>   Elephas_maximus, Delphinus_delphis, Felis_silvestris, Gallus_gallus
#> Node labels:
#>   n1, n2, n3
#> 
#> Rooted; includes branch lengths.
```

`summary_format = "phylo_all"`


```r

datelife::datelife_search(input = my_datelife_query,
                          summary_format = "phylo_all")
#> ... Running a DateLife search.
#> ... Getting a DateLife result.
#> DateLife result obtained!
#> Source chronograms from:
#> 1: Bininda-Emonds, Olaf R. P., Marcel Cardillo, Kate E. Jones, Ross D. E. MacPhee, Robin M. D. Beck, Richard Grenyer, Samantha A. Price, Rutger A. Vos, John L. Gittleman, Andy Purvis. 2007. The delayed rise of present-day mammals. Nature 446 (7135): 507-512
#> 2: Bininda-Emonds, Olaf R. P., Marcel Cardillo, Kate E. Jones, Ross D. E. MacPhee, Robin M. D. Beck, Richard Grenyer, Samantha A. Price, Rutger A. Vos, John L. Gittleman, Andy Purvis. 2007. The delayed rise of present-day mammals. Nature 446 (7135): 507-512
#> 3: Bininda-Emonds, Olaf R. P., Marcel Cardillo, Kate E. Jones, Ross D. E. MacPhee, Robin M. D. Beck, Richard Grenyer, Samantha A. Price, Rutger A. Vos, John L. Gittleman, Andy Purvis. 2007. The delayed rise of present-day mammals. Nature 446 (7135): 507-512
#> 4: Hedges, S. Blair, Julie Marin, Michael Suleski, Madeline Paymer, Sudhir Kumar. 2015. Tree of life reveals clock-like speciation and diversification. Molecular Biology and Evolution 32 (4): 835-845
#> Input taxa presence across source chronograms:
#>               taxon chronograms
#> 1   Elephas maximus         4/4
#> 2 Delphinus delphis         4/4
#> 3  Felis silvestris         4/4
#> 4     Gallus gallus         1/4
#> Input taxa completely absent from source chronograms:
#>   taxon
#> 1  None
#> DateLife search done!
#> 4 phylogenetic trees
```

`summary_format = "phylo_biggest"`


```r

datelife::datelife_search(input = my_datelife_query,
                          summary_format = "phylo_biggest")
#> ... Running a DateLife search.
#> ... Getting a DateLife result.
#> DateLife result obtained!
#> Source chronograms from:
#> 1: Bininda-Emonds, Olaf R. P., Marcel Cardillo, Kate E. Jones, Ross D. E. MacPhee, Robin M. D. Beck, Richard Grenyer, Samantha A. Price, Rutger A. Vos, John L. Gittleman, Andy Purvis. 2007. The delayed rise of present-day mammals. Nature 446 (7135): 507-512
#> 2: Bininda-Emonds, Olaf R. P., Marcel Cardillo, Kate E. Jones, Ross D. E. MacPhee, Robin M. D. Beck, Richard Grenyer, Samantha A. Price, Rutger A. Vos, John L. Gittleman, Andy Purvis. 2007. The delayed rise of present-day mammals. Nature 446 (7135): 507-512
#> 3: Bininda-Emonds, Olaf R. P., Marcel Cardillo, Kate E. Jones, Ross D. E. MacPhee, Robin M. D. Beck, Richard Grenyer, Samantha A. Price, Rutger A. Vos, John L. Gittleman, Andy Purvis. 2007. The delayed rise of present-day mammals. Nature 446 (7135): 507-512
#> 4: Hedges, S. Blair, Julie Marin, Michael Suleski, Madeline Paymer, Sudhir Kumar. 2015. Tree of life reveals clock-like speciation and diversification. Molecular Biology and Evolution 32 (4): 835-845
#> Input taxa presence across source chronograms:
#>               taxon chronograms
#> 1   Elephas maximus         4/4
#> 2 Delphinus delphis         4/4
#> 3  Felis silvestris         4/4
#> 4     Gallus gallus         1/4
#> Input taxa completely absent from source chronograms:
#>   taxon
#> 1  None
#> DateLife search done!
#> 
#> Phylogenetic tree with 4 tips and 3 internal nodes.
#> 
#> Tip labels:
#>   Gallus_gallus, Felis_silvestris, Delphinus_delphis, Elephas_maximus
#> 
#> Rooted; includes branch lengths.
```

`summary_format = "html"`


```r

datelife::datelife_search(input = my_datelife_query,
                          summary_format = "html")
#> ... Running a DateLife search.
#> ... Getting a DateLife result.
#> DateLife result obtained!
#> Input taxa presence across source chronograms:
#>               taxon chronograms
#> 1   Elephas maximus         4/4
#> 2 Delphinus delphis         4/4
#> 3  Felis silvestris         4/4
#> 4     Gallus gallus         1/4
#> Input taxa completely absent from source chronograms:
#>   taxon
#> 1  None
#> DateLife search done!
#> [1] "<table border='1'><tr><th>MRCA Age (MY)</th><th>Ntax</th><th>Citation</th><th>Newick</th></tr><tr><td>101.3</td><td>3</td><td>Bininda-Emonds, Olaf R. P., Marcel Cardillo, Kate E. Jones, Ross D. E. MacPhee, Robin M. D. Beck, Richard Grenyer, Samantha A. Price, Rutger A. Vos, John L. Gittleman, Andy Purvis. 2007. The delayed rise of present-day mammals. Nature 446 (7135): 507-512</td><td>((Felis_silvestris:88.5,Delphinus_delphis:88.5):12.8,Elephas_maximus:101.3);</td></tr><tr><td>96.1</td><td>3</td><td>Bininda-Emonds, Olaf R. P., Marcel Cardillo, Kate E. Jones, Ross D. E. MacPhee, Robin M. D. Beck, Richard Grenyer, Samantha A. Price, Rutger A. Vos, John L. Gittleman, Andy Purvis. 2007. The delayed rise of present-day mammals. Nature 446 (7135): 507-512</td><td>((Delphinus_delphis:86.8,Felis_silvestris:86.8):9.3,Elephas_maximus:96.1);</td></tr><tr><td>110.8</td><td>3</td><td>Bininda-Emonds, Olaf R. P., Marcel Cardillo, Kate E. Jones, Ross D. E. MacPhee, Robin M. D. Beck, Richard Grenyer, Samantha A. Price, Rutger A. Vos, John L. Gittleman, Andy Purvis. 2007. The delayed rise of present-day mammals. Nature 446 (7135): 507-512</td><td>((Felis_silvestris:89.8,Delphinus_delphis:89.8):21,Elephas_maximus:110.8);</td></tr><tr><td>320.5238235</td><td>4</td><td>Hedges, S. Blair, Julie Marin, Michael Suleski, Madeline Paymer, Sudhir Kumar. 2015. Tree of life reveals clock-like speciation and diversification. Molecular Biology and Evolution 32 (4): 835-845</td><td>((Elephas_maximus:105.009658,(Felis_silvestris:78.392518,Delphinus_delphis:78.392517):26.617141):215.5141645,Gallus_gallus:320.5238235);</td></tr></table>"
#> attr(,"datelife_result")
#> $`Bininda-Emonds, Olaf R. P., Marcel Cardillo, Kate E. Jones, Ross D. E. MacPhee, Robin M. D. Beck, Richard Grenyer, Samantha A. Price, Rutger A. Vos, John L. Gittleman, Andy Purvis. 2007. The delayed rise of present-day mammals. Nature 446 (7135): 507-512`
#>                   Elephas maximus Delphinus delphis
#> Elephas maximus               0.0             202.6
#> Delphinus delphis           202.6               0.0
#> Felis silvestris            202.6             177.0
#>                   Felis silvestris
#> Elephas maximus              202.6
#> Delphinus delphis            177.0
#> Felis silvestris               0.0
#> 
#> $`Bininda-Emonds, Olaf R. P., Marcel Cardillo, Kate E. Jones, Ross D. E. MacPhee, Robin M. D. Beck, Richard Grenyer, Samantha A. Price, Rutger A. Vos, John L. Gittleman, Andy Purvis. 2007. The delayed rise of present-day mammals. Nature 446 (7135): 507-512`
#>                   Elephas maximus Felis silvestris
#> Elephas maximus               0.0            192.2
#> Felis silvestris            192.2              0.0
#> Delphinus delphis           192.2            173.6
#>                   Delphinus delphis
#> Elephas maximus               192.2
#> Felis silvestris              173.6
#> Delphinus delphis               0.0
#> 
#> $`Bininda-Emonds, Olaf R. P., Marcel Cardillo, Kate E. Jones, Ross D. E. MacPhee, Robin M. D. Beck, Richard Grenyer, Samantha A. Price, Rutger A. Vos, John L. Gittleman, Andy Purvis. 2007. The delayed rise of present-day mammals. Nature 446 (7135): 507-512`
#>                   Delphinus delphis Felis silvestris
#> Delphinus delphis               0.0            179.6
#> Felis silvestris              179.6              0.0
#> Elephas maximus               221.6            221.6
#>                   Elephas maximus
#> Delphinus delphis           221.6
#> Felis silvestris            221.6
#> Elephas maximus               0.0
#> 
#> $`Hedges, S. Blair, Julie Marin, Michael Suleski, Madeline Paymer, Sudhir Kumar. 2015. Tree of life reveals clock-like speciation and diversification. Molecular Biology and Evolution 32 (4): 835-845`
#>                   Felis silvestris Delphinus delphis
#> Felis silvestris            0.0000          156.7850
#> Delphinus delphis         156.7850            0.0000
#> Elephas maximus           210.0193          210.0193
#> Gallus gallus             641.0476          641.0476
#>                   Elephas maximus Gallus gallus
#> Felis silvestris         210.0193      641.0476
#> Delphinus delphis        210.0193      641.0476
#> Elephas maximus            0.0000      641.0476
#> Gallus gallus            641.0476        0.0000
#> 
#> attr(,"class")
#> [1] "datelifeResult"
#> attr(,"datelife_query")
#> $cleaned_names
#> [1] "Elephas maximus"   "Delphinus delphis" "Felis silvestris" 
#> [4] "Gallus gallus"    
#> 
#> $ott_ids
#> NULL
#> 
#> $phy
#> 
#> Phylogenetic tree with 4 tips and 3 internal nodes.
#> 
#> Tip labels:
#>   Elephas_maximus, Delphinus_delphis, Felis_silvestris, Gallus_gallus
#> 
#> Rooted; no branch lengths.
#> 
#> attr(,"class")
#> [1] "datelifeQuery"
```

`summary_format = "data_frame"`


```r

datelife::datelife_search(input = my_datelife_query,
                          summary_format = "data_frame")
#> ... Running a DateLife search.
#> ... Getting a DateLife result.
#> DateLife result obtained!
#> Input taxa presence across source chronograms:
#>               taxon chronograms
#> 1   Elephas maximus         4/4
#> 2 Delphinus delphis         4/4
#> 3  Felis silvestris         4/4
#> 4     Gallus gallus         1/4
#> Input taxa completely absent from source chronograms:
#>   taxon
#> 1  None
#> DateLife search done!
#>        Age Ntax
#> 1 101.3000    3
#> 2  96.1000    3
#> 3 110.8000    3
#> 4 320.5238    4
#>                                                                                                                                                                                                                                                         Citation
#> 1 Bininda-Emonds, Olaf R. P., Marcel Cardillo, Kate E. Jones, Ross D. E. MacPhee, Robin M. D. Beck, Richard Grenyer, Samantha A. Price, Rutger A. Vos, John L. Gittleman, Andy Purvis. 2007. The delayed rise of present-day mammals. Nature 446 (7135): 507-512
#> 2 Bininda-Emonds, Olaf R. P., Marcel Cardillo, Kate E. Jones, Ross D. E. MacPhee, Robin M. D. Beck, Richard Grenyer, Samantha A. Price, Rutger A. Vos, John L. Gittleman, Andy Purvis. 2007. The delayed rise of present-day mammals. Nature 446 (7135): 507-512
#> 3 Bininda-Emonds, Olaf R. P., Marcel Cardillo, Kate E. Jones, Ross D. E. MacPhee, Robin M. D. Beck, Richard Grenyer, Samantha A. Price, Rutger A. Vos, John L. Gittleman, Andy Purvis. 2007. The delayed rise of present-day mammals. Nature 446 (7135): 507-512
#> 4                                                           Hedges, S. Blair, Julie Marin, Michael Suleski, Madeline Paymer, Sudhir Kumar. 2015. Tree of life reveals clock-like speciation and diversification. Molecular Biology and Evolution 32 (4): 835-845
#>                                                                                                                                     Newick
#> 1                                                             ((Felis_silvestris:88.5,Delphinus_delphis:88.5):12.8,Elephas_maximus:101.3);
#> 2                                                               ((Delphinus_delphis:86.8,Felis_silvestris:86.8):9.3,Elephas_maximus:96.1);
#> 3                                                               ((Felis_silvestris:89.8,Delphinus_delphis:89.8):21,Elephas_maximus:110.8);
#> 4 ((Elephas_maximus:105.009658,(Felis_silvestris:78.392518,Delphinus_delphis:78.392517):26.617141):215.5141645,Gallus_gallus:320.5238235);
```








