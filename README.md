
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

*A computational pipeline to map TagMap-reads*

TagMap is a very useful method for transposon mapping (Stern 2017),
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

The basic usage of tagMapper revolves around three clear steps:

1.  index: a tagMapper-index is made once for a specific genome and
    protocol (e.g. hg19 and PigyBac).
2.  align: a tagMapperSample-object is made and aligned to the index
3.  analyse: determine and plot highly likely integraton-sites

See the
[vignette](https://raw.githubusercontent.com/robinweide/tagmeppr/master/vignettes/tagmeppr.pdf)
for a more in-depth coverage of all things tagMeppr\!

## Code of conduct

Please note that this project is released with a [Contributor Code of
Conduct](.github/CODE_OF_CONDUCT.md). By participating in this project
you agree to abide by its terms.

-----

<div id="refs" class="references">

<div id="ref-Stern037762">

Stern, David L. 2017. “Tagmentation-Based Mapping (Tagmap) of Mobile Dna
Genomic Insertion Sites.” *bioRxiv*. Cold Spring Harbor Laboratory.
<https://doi.org/10.1101/037762>.

</div>

</div>
