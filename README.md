# EasyCellType
EasyCellType is a R package developed to automatically examine the input marker
lists obtained from existing software such as Seurat over the cell marker databases.
Two quantification approaches are provided to annotate cell types: Gene set 
enrichment analysis (GSEA) and a fisher exact test combining with expression 
significance. The package generates straight-forward visualization of the test results, 
which help the researchers to understand the results and pick the optimal annotations.

# Installation
```
library(devtools)
devtools::install_github("rx-li/EasyCellType")
```
