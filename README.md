
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RHermes

<!-- badges: start -->

[![Codecov test
coverage](https://codecov.io/gh/RogerGinBer/RHermes/branch/master/graph/badge.svg?token=HL73R4GHFJ)](https://codecov.io/gh/RogerGinBer/RHermes?branch=master)

<!-- badges: end -->

<img align = "right" style = "padding-left: 10%; padding-bottom: 10%; padding-top: 2%" 
width = "250px" height = "285px" src = "https://i.postimg.cc/Tw9SvJ11/sticker-No-Bioc.png">

`RHermes` is a broad-scoped targeted metabolomics software designed to
analyse LC-MS and LC-MS/MS data to identify compounds in biological or
environmental samples.

The `RHermes` workflow works with both Orbitrap and q-TOF instrument
data and comes with an easy to use GUI that will guide you every step of
the way.

You are in control of your metabolites: whether it’s natural products,
biomedical or enviormental samples, `RHermes` has you covered. By
restricting the formula database, you can focus on just the compounds
you are interested in and achieve greater metabolome coverage depth.

Have you ever wished you could just **see** the metabolites in your
data? With `RHermes` you can do that and much more. Say goodbye to
manually calculating m/z’s and plotting XIC of different adducts: with
our GUI you are just one click away from a metabolite-centric plot.

For more info, check out the documentation
[here](https://rogerginber.github.io/RHermes/) and our article preprint
[here](https://www.biorxiv.org/content/10.1101/2021.03.08.434466v1.full.pdf)

## System requirements

The recommended system requirements are:

-   At least 8-16 GB of RAM
-   An internet connection to perform KEGG queries

## Installation

You can download the development version from
[GitHub](https://github.com/RogerGinBer/RHermes) with:

``` r
if(!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("RogerGinBer/RHermes")
```

## Setup

`RHermes` can perform almost all its functions after installation, but
the SOI Blank Substraction step requires a valid `keras` and
`tensorflow` installation (which rely on Python).

### Option 1: Default installation

For most users, both `keras` and `tensorflow` can be configured with:

``` r
reticulate::install_miniconda()
keras::install_keras()
tensorflow::install_tensorflow()
```

After which you can check the following:

``` r
tensorflow::tf_config()
model <- keras::load_model_hdf5(system.file("extdata", "ImprovedModel.h5",
                                            package = "RHermes"))
is(model, "python.builtin.object") #Gives TRUE if the loading is successful.
```

If both commands don’t give any error (the “Your CPU supports …” warning
is fine) the installation has been successful. If it fails (which can
happen in some users with previous Python installations, try Option 2).

### Option 2: Manual installation

First install Miniconda:

``` r
reticulate::install_miniconda()
```

Now find the Miniconda Prompt in your computer. Instead of relying on
the default r-reticulate environment, type the following to create a new
environment:

`conda create -n newenv python=3.6 tensorflow keras`

When finished, type in R:

``` r
reticulate::use_condaenv("newenv", required = TRUE)
tensorflow::tf_config()
model <- keras::load_model_hdf5(system.file("extdata", "ImprovedModel.h5",
                                            package = "RHermes"))
```

Everything should run smoothly. If not, try manually installing Anaconda
from their website and letting `reticulate` know where to find the
environment.

Also check out
[Keras](https://tensorflow.rstudio.com/tutorials/beginners/basic-ml/)
and [Tensorflow](https://tensorflow.rstudio.com/tutorials/) R tutorials.

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

## Bug reporting

Suggestions and bug reports are more than welcome at:
<https://github.com/RogerGinBer/RHermes/issues>

## Citation

Please cite this package as:

HERMES: a molecular formula-oriented method to target the metabolome
Roger Giné, Jordi Capellades, Josep M. Badia, Dennis Vughs, Michaela
Schwaiger-Haber, Maria Vinaixa, Andrea M. Brunner, Gary J. Patti, Oscar
Yanes bioRxiv 2021.03.08.434466; doi:
<https://doi.org/10.1101/2021.03.08.434466>
