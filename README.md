
# scDNS

<!-- badges: start -->
<!-- badges: end -->

scDNS is a computational framework that quantifies gene-specific functional perturbations at single-cell resolution by computing information-theoretic divergence to compare gene network configurations across biological conditions

## Installation

You can install the development version of scDNS like so:
```r
devtools::install_github('xiaolab-xjtu/scDNS')
```

``` r
# FILL THIS IN! HOW CAN PEOPLE INSTALL YOUR DEV PACKAGE?
```

## An example for scDNS

This is a basic example which shows you how to solve a common problem:
### Load data and preprocessing
``` r
library(scDNS)
library(patchwork)
library(dplyr)
library(reticulate)
library(Seurat)
library(Rmagic)
use_condaenv('r-reticulate')
options(future.globals.maxSize = 2 * 1024^3)
#02 load scRNA-seq data with Seurat object----------

load(file = './GEM24h_scDNS.RData')

PANC1GEM24H_SOB_3 <- rmNotExpressedGene(PANC1GEM24H_SOB_3)

sob <- PANC1GEM24H_SOB_3

rm(PANC1GEM24H_SOB_3)
gc()

DefaultAssay(sob) <- 'RNA'
group_col <- 'Type'  # conditions
replicates_col <- 'Type'

#impute gene expression
if (!'MAGIC_RNA'%in%names(sob@assays)){
  sob <- Magic4MultipleData(sob,split.by = group_col)
}
```

### 1 Creat scDNS object
``` r
scDNSob <- seurat2scDNSObj(sob = sob,
                           GroupBy = group_col,
                           imputedAssay = 'MAGIC_RNA',
                           parallel.sz = 10,
                           loop.size=6000)
```
						   
#### 2 Calculate network divergence
``` r
scDNSob <- scDNS_1_CalDivs(scDNSob)
```

#### 3 Creat network ensmeble modle to combine  network divergence

``` r
scDNSob <- scDNS_2_creatNEAModel_v2(scDNSobjcet = scDNSob,
                                    n.randNet = 20000)
```
#### 4 Construct context-adaptive GINs using a GAT based model

#### 5 Calculate gene perpturbation score (Gene Z-score)
``` r
scDNSob <- scDNS_2_creatNEAModel_v2(scDNSobjcet = scDNSob,
                                    n.randNet = 20000)
```
#### 6 Calculate single cell perptrubation score
``` r
scDNSob <- scDNS_4_scContribution(scDNSob)
```

