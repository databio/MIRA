# MIRA (Methylation-based Inference of Regulatory Activity)
[![Build Status](https://travis-ci.org/databio/MIRA.svg?branch=master)](https://travis-ci.org/databio/MIRA)

The MIRA package aggregates DNA methylation across the genome for instances of a genomic feature like histone ChIP peaks, transcription factor ChIP peaks, or open chromatin regions in order to give a single signature and score for that feature, which may be used to infer regulatory activity.


### Installing MIRA

MIRA may be installed from Github:

```{r}
devtools::install_github("databio/MIRA")
```

or locally after downloading/cloning the source code:

```{r}
install.packages("path/to/MIRA/directory", repos=NULL, type="source")
```


### Learning how to use MIRA

A couple of vignettes are included with the package. For an overview of the package see the "Getting Started" vignette. A more realistic example of how to use MIRA may be found in the "BiologicalApplication" vignette.
