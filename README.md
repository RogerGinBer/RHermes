
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RHermes

<!-- badges: start -->

[![Codecov test
coverage](https://codecov.io/gh/RogerGinBer/Proves/branch/master/graph/badge.svg)](https://codecov.io/gh/RogerGinBer/Proves?branch=master)
<!-- badges: end -->

The goal of RHermes is to analyse LC-MS and LC-MS/MS files to identify
compounds in biological or environmental samples. The RHermes workflow
can work with either Orbitrap or q-TOF instrument data, but the former
will yield best results.

## System requirements

RHermes is quite heavy on the CPU and memory loads of your system, as it
has been developed to use almost all available cores. For that reason
you will need:

  - At least a 4 core processor
  - 16 GB of RAM or more
  - An internet connection to perform KEGG queries

## Installation

You can install the released version of RHermes from
[Bioconductor](https://bioconductor.org/) with:

``` r
# install.packages("BiocManager")
BiocManager::install("RHermes")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("RogerGinBer/Proves")
```

(REMOVE LATER)

There’s yet another way to install RHermes. Clone the RHermes repo and
enter the .RProj file. Build the package with Ctrl+Shift+B and manually
install the dependencies with either install.packages() for CRAN or
BiocManager::install() for Bioconductor packages.

Once you’ve installed all dependencies, build the package again and it
will be installed.

## Setup

RHermes can perform almost all its routines out of the box, but the SOI
Blank Substraction step requires a valid Keras and Tensorflow
installation, which rely on Python.

### Option 1: Default installation

In theory both Keras and Tensorflow can be configured with:

``` r
reticulate::install_miniconda()
keras::install_keras()
tensorflow::install_tensorflow()
```

After which you can check the following:

``` r
tensorflow::tf_config()
model <- keras::load_model_tf(system.file("extdata", "model", package = "RHermes"))
```

If both commands don’t yield any error (the “Your CPU supports …”
warning is fine) you are set to go. If it fails (which can happen in
some devices, try Option 2).

### Option 2: Dealing with the Installation yourself

First install Miniconda:

``` r
reticulate::install_miniconda()
```

Now find the Miniconda Prompt in your computer. Instead of relying on
the default r-reticulate environment, type the following to create a new
environment:

conda create -n newenv python=3.6 tensorflow keras

When finished, type in R:

``` r
reticulate::use_condaenv("newenv", required = TRUE)
tensorflow::tf_config()
model <- keras::load_model_tf(system.file("extdata", "model", package = "RHermes"))
```

Everything should run smoothly. If not, try manually installing Anaconda
from their website and letting reticulate know where to find the
environment.

Also check out Keras and Tensorflow R tutorials.

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(RHermes)
#> Warning: replacing previous import 'igraph::groups' by 'plotly::groups' when
#> loading 'RHermes'
#> Warning: replacing previous import 'ggplot2::last_plot' by 'plotly::last_plot'
#> when loading 'RHermes'
#> Warning: replacing previous import 'data.table::between' by 'dplyr::between'
#> when loading 'RHermes'
## basic example code
```
