
# scDNS

<!-- badges: start -->
<!-- badges: end -->

The goal of scDNS2 is to ...

## Installation

You can install the development version of scDNS like so:



``` r
# FILL THIS IN! HOW CAN PEOPLE INSTALL YOUR DEV PACKAGE?
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(scDNS)
library(patchwork)
library(dplyr)
library(reticulate)
library(Seurat)
library(Rmagic)
use_condaenv('r-reticulate')
options(future.globals.maxSize = 2 * 1024^3)
#02 加载单细胞数据Seurat对象----------

load(file = '/home/rstudio/Projects/scDDG/scDNS_application/GEM24h_scDNS.RData')

PANC1GEM24H_SOB_3 <- rmNotExpressedGene(PANC1GEM24H_SOB_3)

sob <- PANC1GEM24H_SOB_3

rm(PANC1GEM24H_SOB_3)
gc()
DefaultAssay(sob) <- 'RNA'
# 03设置项目相关的变量--------------
ProjectName <- 'PADC_GEM_CellLine'  # 项目名称
group_col <- 'Type'  # 比较组的列名
cellType_col <- 'seurat_clusters'  # 细胞类型的列名
replicates_col <- 'Type'
savePath <- '/home/rstudio/Projects/scDDG/Application/GEM/'  # 结果保存路径

# sob <- Magic4MultipleData(sob,split.by = 'Type')

# # # 03scDNS--------------
# #
# 0.0 impute gene expression
if (!'MAGIC_RNA'%in%names(sob@assays)){
  sob <- Magic4MultipleData(sob,split.by = group_col)
}
# 1 creat scDNS object
scDNSob <- seurat2scDNSObj(sob = sob,
                           GroupBy = group_col,
                           imputedAssay = 'MAGIC_RNA',
                           parallel.sz = 10,
                           loop.size=6000)
# 2 calculate network divergence
scDNSob <- scDNS_1_CalDivs(scDNSob)
# 3 creat network ensmeble modle to combine  network divergence
scDNSob <- scDNS_2_creatNEAModel_v2(scDNSobjcet = scDNSob,
                                    n.randNet = 20000)
# 4 construct context-adaptive GINs 

# 5 calculate gene perpturbation score (Gene Z-score)
scDNSob <- scDNS_3_GeneZscore_v2(scDNSob)

# 6 calculate single cell perptrubation score
scDNSob <- scDNS_4_scContribution(scDNSob)


```

