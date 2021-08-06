# Updates

## SpotClean 0.99.0

---------------------

* Created initial github page.
* Created package.
* Created vignettes.
* Checked package... passed R CMD check and BiocCheck.
* Ready to submit.

## SpotClean 0.99.1

---------------------

* Fixed a bug due to genes with zero expression in tissue spots but positive expression in background spots.
* As a result, when creating the slide object, genes will be filtered by average expression in tissue spots, not all spots.
