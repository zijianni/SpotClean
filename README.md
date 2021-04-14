# TBD

Droplet-based single cell RNA-seq technologies provide a novel insight in transcriptome profiles of individual cells with a much larger cell numbers and cheaper cost compared with microfluidic-based single cell RNA-seq. During library preparation, each cell is expected to be captured by one droplet. The number of droplets are usually much more than the number of cells, thus most droplets do not contain real cells. However, there are always free-floating RNA fragments in the library due to broken cells or library contamination. Empty droplets will capture them and have non-zero expression values after sequencing.

**CB2** is a cluster-based approach for distinguishing true cells from background barcodes in droplet-based single cell RNA-seq experiments (especially for 10X Chromium output), while `scCB2` is its corresponding R package. It is based on clustering similar barcodes and calculating Monte-Carlo p-value for each cluster to test against background distribution. This cluster-level test outperforms single-barcode-level tests not only for high count barcodes, but also in dealing with low count barcodes and homogeneous sequencing library, while keeping FDR well controlled.

This package has been accepted by Bioconductor. For more instructions, see [here](https://bioconductor.org/packages/3.12/bioc/html/scCB2.html).

## Installation

Install via Bioconductor: 

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version='devel')
BiocManager::install("scCB2")
```


Alternatively, install via Github using `devtools`:

```
if(!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("zijianni/scCB2", build_manual = TRUE, build_vignettes = TRUE)
```

Note: This may take a few minutes for building vignettes. If you don't need vignettes (which rarely happens), set `build_vignettes = FALSE`.

After installing, you will find package vignettes by `vignette("scCB2")`.


## Vignettes 
[Link to Bioconductor](https://bioconductor.org/packages/3.12/bioc/vignettes/scCB2/inst/doc/scCB2.html)



## Citation

Ni, Z., Chen, S., Brown, J., & Kendziorski, C. (2020). CB2 improves power of cell detection in droplet-based single-cell RNA sequencing data. Genome Biology, 21(1), 137. https://doi.org/10.1186/s13059-020-02054-8
