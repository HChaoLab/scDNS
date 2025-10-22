# scDNS

**scDNS** is a computational framework that quantifies **gene-specific functional perturbations at single-cell resolution**.
By leveraging **information-theoretic divergence**, scDNS systematically compares **gene network configurations** across different biological conditions, enabling precise characterization of cell-state-specific functional alterations.

---

## 🧩 Installation

Install the development version of **scDNS** from GitHub using:

```r
devtools::install_github("xiaolab-xjtu/scDNS")
```

---

## 🚀 Quick Start Example

This example demonstrates a typical analysis workflow using **scDNS** on scRNA-seq data.

### 1. Load and preprocess data

```r
library(scDNS)
library(patchwork)
library(dplyr)
library(reticulate)
library(Seurat)
library(Rmagic)

# Activate Python environment
use_condaenv("r-reticulate")

# Increase global object size limit
options(future.globals.maxSize = 2 * 1024^3)

# Load example Seurat object
load(file = "./GEM24h_scDNS.RData")

# Remove non-expressed genes
PANC1GEM24H_SOB_3 <- rmNotExpressedGene(PANC1GEM24H_SOB_3)
sob <- PANC1GEM24H_SOB_3
rm(PANC1GEM24H_SOB_3)
gc()

DefaultAssay(sob) <- "RNA"
group_col <- "Type"       # experimental condition
replicates_col <- "Type"  # replicate label

# Impute gene expression using MAGIC if not already available
if (!"MAGIC_RNA" %in% names(sob@assays)) {
  sob <- Magic4MultipleData(sob, split.by = group_col)
}
```

---

### 2. Create a `scDNS` object

```r
scDNSob <- seurat2scDNSObj(
  sob = sob,
  GroupBy = group_col,
  imputedAssay = "MAGIC_RNA",
  parallel.sz = 10,
  loop.size = 6000
)
```

---

### 3. Compute network divergence

```r
scDNSob <- scDNS_1_CalDivs(scDNSob)
```

This step quantifies the **information-theoretic divergence** between gene networks across conditions.

---

### 4. Build a network ensemble model

```r
scDNSob <- scDNS_2_creatNEAModel_v2(
  scDNSobject = scDNSob,
  n.randNet = 20000
)
```

A **network ensemble model** is constructed to integrate divergence values and estimate background distributions using random network permutations.

---

### 5. Construct context-adaptive GINs (optional)

This step builds **context-adaptive gene interaction networks (GINs)** using a **graph attention network (GAT)** model, capturing nonlinear dependencies between genes across conditions.

---

### 6. Calculate gene-level perturbation scores

```r
scDNSob <- scDNS_3_GeneZscore_v2(scDNSobject = scDNSob)
```

This step computes **gene perturbation Z-scores**, representing the degree of network rewiring for each gene.

---

### 7. Calculate cell-level perturbation scores

```r
scDNSob <- scDNS_4_scContribution(
  scDNSobject,
  topGene = 100,
  InterstingGene = NULL,
  q.th = 0.1
)
```

This step summarizes gene-level perturbations into **single-cell perturbation scores (scZscores)**, allowing visualization of functional heterogeneity across cells.

---

## 🧠 Citation

If you use **scDNS** in your research, please cite:

>  *scDNS: Characterizing Gene Perturbations in Single Cells via Network Divergence Analysis.* (2025)

---

## 📘 License

This package is released under the MIT License.

---
