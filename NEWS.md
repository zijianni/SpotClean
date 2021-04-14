# Updates

## scCB2 0.99.0

---------------------

* Created initial github page.
* Created package.
* Created vignettes.
* Checked package... passed R CMD check and BiocCheck.
* Ready to submit to Bioconductor.

## scCB2 0.99.11

---------------------

* Customized knee point calculation function.
* CB2FindCell will print out retain threshold in standard output.
* Minor bug fixes when there is no cluster.

## scCB2 0.99.12

---------------------

* Minor changes on vignette.
* Baseline (2 * background) clustering threshold will be printed out.

## scCB2 0.99.13

---------------------

* Rounded baseline threshold.
* Minor changes on vignette.
* Bug fix when retain threshold is larger than any barcode.

## scCB2 0.99.14

---------------------

* Change parameter names to match those in the paper: retain -> upper. background -> lower.

## scCB2 0.99.15

---------------------

* When dividing barcodes into groups with similar barcode counts, the last group will be combined with second last group if the number of barcodes in the last group is less than half of that in the second last group. 

## scCB2 0.99.16

---------------------

* Change function name `Read10X` to be `Read10xRaw`, `Read10Xh5` to be `Read10xRawH5`.
* Change parameter name `PrintProg` to be `verbose`.
* Minor edits on parameter description of function `GetCellMat`.
* Added a quick function `QuickCB2` by combining all necessary functions into one. Input: directory of raw data. Output: filtered matrix, or a Seurat object containing filtered matrix.

## scCB2 0.99.17

---------------------

* Initial version preparing to submit to Bioconductor.

## scCB2 0.99.18

---------------------

* Minor bug fix.
* Minor edits on package vignettes.

## scCB2 0.99.19

---------------------

* Minor edits on package vignettes.

## scCB2 0.99.20

---------------------

* Updated citation.

## scCB2 0.99.21

---------------------

* Updated testthat scripts.
* Minor bug fix in Read10xRaw.R. Added as.numeric when building sparse matrix.

## scCB2 0.99.22

---------------------

* Received comments from Bioconductor reviewer.
* Moved `BiocStyle` from `Imports:` to `Suggests:`.
* Modified vignettes. Deleted GitHub installation instruction.
* Changed default value of `Ncores` to be 2.
* Changed output of `CB2FindCell` to be an `SummarizedExperiment` object. Changed `GetCellMat` accordingly. Changed examples and testthat accordingly.

## scCB2 0.99.23

---------------------

* Minor edit of README.md

## scCB2 0.99.24

---------------------

* Removed direct download files in vignettes to speed up building.
* Changed R dependencies so that package can be built in older R versions.
* Changed vignette title to match the citation.

## scCB2 0.99.25

---------------------

* Changed R dependency back to 4.0.0.

## scCB2 0.99.29

---------------------

* Updated functions of calculating Pearson correlation in sparse matrix. New functions are more accurate.

## scCB2 0.99.30

---------------------

* Added thresholding of number of top genes to use in testing steps. This avoids high number of false positives in ultra-high dimensional datasets, e.g. 10x barnyard data. 

## scCB2 0.99.31

---------------------

* Added extra filtering for barnyard data (a mixture of multiple species). 
