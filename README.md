
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RHermes

<!-- badges: start -->

[![Codecov test
coverage](https://codecov.io/gh/RogerGinBer/RHermes/branch/master/graph/badge.svg?token=HL73R4GHFJ)](https://codecov.io/gh/RogerGinBer/RHermes?branch=master)
[![R-CMD-check](https://github.com/RogerGinBer/RHermes/actions/workflows/check-bioc.yml/badge.svg)](https://github.com/RogerGinBer/RHermes/actions/workflows/check-bioc.yml)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

<img align = "right" style = "padding-left: 10%; padding-bottom: 10%; padding-top: 2%" 
width = "250px" height = "285px" src = "https://i.postimg.cc/Tw9SvJ11/sticker-No-Bioc.png">

`RHermes` is a broad-scoped targeted metabolomics package designed to
analyse LC-MS and LC-MS/MS data to identify compounds in biological or
environmental samples.

The `RHermes` workflow works with both Orbitrap and qTOF instrument data
and comes with an interactive GUI to process and visualize the data.

Whether it’s natural products, biomedical or enviormental samples,
`RHermes` can help you improve your matrix characterization. By
selecting an appropiate formula database, you can focus on just the
compounds you are interested in and improve your coverage.

With `RHermes` you can **see** the metabolites in your data and much
more. There’s no need to manually calculate m/z’s to plot the XIC of
different adducts: with the GUI you are just one click away from a
metabolite-centric plot.

For more info, check out the documentation
[here](https://rogerginber.github.io/RHermes/) and our article
[here](https://www.nature.com/articles/s41592-021-01307-z).

## Installation

You can download the development version from
[GitHub](https://github.com/RogerGinBer/RHermes) with:

``` r
if(!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("RogerGinBer/RHermes")
```

## Analyzing LC-MS data with RHermes

Once installed, you can use `RHermes` programmatically like this:

``` r
library(RHermes)
#Generate a Exp object
example <- RHermesExp()

#Set your formula and adduct database
example <- setDB(example, db = "hmdb")

#Process your MS1 files
example <- processMS1(example,
                        system.file("extdata", "MS1TestData.mzML",
                        package = "RHermes"))
#Generate SOIs
example <- findSOI(example, getSOIpar(), 1)

#Generate an IL (Inclusion List)
example <- generateIL(example, 1, ILParam())
```

With the generated inclusion list, you can export it and run a Parallel
Reaction Monitoring (PRM) MS2 experiment to reveal coeluting isomers or
use any other MS2 mode you see fit.

Or start the interactive GUI typing:

``` r
RHermesGUI()
```

In the GUI you will find abundant help pages to guide you along the
processing :+1:

Please check the User Guide
[vignette](https://rogerginber.github.io/RHermes/articles/RHermes_UserGuide.html)
for more detailed info and real examples.

## Database availability

We have compiled some molecular formula open databases ready to be used
with RHermes for all sorts of samples: HMDB, ChEBI, NORMAN, LipidMaps,
COCONUT, etc. They are freely available at [this Zenodo
repository](https://zenodo.org/record/5025560).

## Bug reporting

Suggestions and bug reports are more than welcome at:
<https://github.com/RogerGinBer/RHermes/issues>

## Citation

Please cite this package as:

Giné, R., Capellades, J., Badia, J.M. et al. HERMES: a
molecular-formula-oriented method to target the metabolome. Nat Methods
18, 1370–1376 (2021). <https://doi.org/10.1038/s41592-021-01307-z>
