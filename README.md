# SpotClean: a computational method to adjust for spot swapping in spatial transcriptomics data

SpotClean is a computational method to adjust for spot swapping in spatial transcriptomics data. Recent spatial transcriptomics experiments utilize slides containing thousands of spots with spot-specific barcodes that bind mRNA. Ideally, unique molecular identifiers at a spot measure spot-specific expression, but this is often not the case due to bleed from nearby spots, an artifact we refer to as spot swapping. SpotClean is able to estimate the contamination rate in observed data and decontaminate the spot swapping effect, thus increase the sensitivity and precision of downstream analyses.

### Introduction

Spatial transcriptomics (ST), named [Method of the Year 2020](https://www.nature.com/articles/s41592-020-01033-y) by *Nature Methods* in 2020, is a powerful and widely-used experimental method for profiling genome-wide gene expression across a tissue. In a typical ST experiment, fresh-frozen (or FFPE) tissue is sectioned and placed onto a slide containing spots, with each spot containing millions of capture oligonucleotides with spatial barcodes unique to that spot. The tissue is imaged, typically via Hematoxylin and Eosin (H&E) staining. Following imaging, the tissue is permeabilized to release mRNA which then binds to the capture oligonucleotides, generating a cDNA library consisting of transcripts bound by barcodes that preserve spatial information. Data from an ST experiment consists of the tissue image coupled with RNA-sequencing data collected from each spot. A first step in processing ST data is tissue detection, where spots on the slide containing tissue are distinguished from background spots without tissue. Unique molecular identifier (UMI) counts at each spot containing tissue are then used in downstream analyses.

Ideally, a gene-specific UMI at a given spot would represent expression of that gene at that spot, and spots without tissue would show no (or few) UMIs. This is not the case in practice. Messenger RNA bleed from nearby spots causes substantial contamination of UMI counts, an artifact we refer to as spot swapping. On average, we observe that more than 30% of UMIs at a tissue spot did not originate from this spot, but from other spots contaminating it. Spot swapping confounds downstream inferences including normalization, marker gene-based annotation, differential expression and cell type decomposition.

We developed **SpotClean** to adjust for the effects of spot swapping in ST experiments. SpotClean is able to measure the per-spot contamination rates in observed data and decontaminate gene expression levels, thus increases the sensitivity and precision of downstream analyses. Our package `SpotClean` is built based on 10x Visium spatial transcriptomics experiments, currently the most widely-used commercial protocol, providing functions to load raw spatial transcriptomics data from 10x Space Ranger outputs, decontaminate the spot swapping effect, estimate contamination levels, visualize expression profiles and spot labels on the slide, and connect with other widely-used packages for further analyses. SpotClean can be potentially extended to other spatial transcriptomics data as long as the gene expression data in both tissue and background regions are available.


### Installation

Install the GitHub version:

```{r}
if(!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("zijianni/SpotClean", build_manual = TRUE, build_vignettes = TRUE)

```

Install the Bioconductor version:

```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("SpotClean")

```

Load package after installation:

```{r}
library(SpotClean)
```

### Tutorial

After installing the package, access the vignette by running

```{r}
vignette("SpotClean")
```

### Citation

We appreciate it if you could cite our work when using `SpotClean`:

Ni, Z., Prasad, A., Chen, S. et al. SpotClean adjusts for spot swapping in spatial transcriptomics data. *Nat Commun* 13, 2971 (2022). https://doi.org/10.1038/s41467-022-30587-y

A BibTeX entry for LaTeX users can be found by running 

```{r}
citation("SpotClean")
```
