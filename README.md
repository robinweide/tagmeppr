
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tagMeppr <img src="vignettes/logo_tagMeppr.png" align="right" alt="" width="120" />

![GitHub](https://img.shields.io/github/license/robinweide/tagMeppr)
[![Build
Status](https://travis-ci.org/robinweide/tagmeppr.svg?branch=master)](https://travis-ci.org/robinweide/tagmeppr)
[![Project Status: WIP – Initial development is in progress, but there
has not yet been a stable, usable release suitable for the
public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![minimal R
version](https://img.shields.io/badge/R%3E%3D-3.4.4-blue.svg)](https://cran.r-project.org/)

*A computational pipeline to map tagmap-insertions*

TagMap is a very useful methods for transposon mapping (Stern 2017),
enabling researchers to map the insertion sites with ease and generate
long sequencing reads. However, there is little to none automatisation
and downstream analysis software available for these reads. TagMeppr is
an easy to use, memory efficient fastq-to-figure package written in R.

## Installation

``` r
# Install development version from GitHub
devtools::install_github("robinweide/tagmeppr")
```

## Usage

To align the tagMeppr-sample to the index, just use the
`align()`-function.

``` r
align(exp = C9, 
      ref = reference_hg19_PB, 
      cores = 30, 
      empericalCentre = T)
```

## Code of conduct

Please note that this project is released with a [Contributor Code of
Conduct](.github/CODE_OF_CONDUCT.md). By participating in this project
you agree to abide by its terms.

<div id="refs" class="references">

<div id="ref-Stern037762">

Stern, David L. 2017. “Tagmentation-Based Mapping (Tagmap) of Mobile Dna
Genomic Insertion Sites.” *bioRxiv*. Cold Spring Harbor Laboratory.
<https://doi.org/10.1101/037762>.

</div>

</div>
