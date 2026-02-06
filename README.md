# Single-Cell RNA-Seq Analysis of Tumor vs. Tumor Microenvironment 

single cell RNA-seq project 2022; updated 2026

## Overview

This project analyzes single-cell RNA sequencing data from tumor (T) and tumor microenvironment (TME) samples across three patients (JMY, WHG, XDL) to investigate inter-patient tumor heterogeneity and characterize immunosuppression in the TME. The analysis follows the standard Seurat workflow with SCTransform normalization.

## Data

Six samples from three patients, each contributing a tumor and TME sample:

| Sample | Condition | Patient |
|--------|-----------|---------|
| JMY_Ca | T         | JMY     |
| JMY_N  | TME       | JMY     |
| WHG_Ca | T         | WHG     |
| WHG_N  | TME       | WHG     |
| XDL_Ca | T         | XDL     |
| XDL_N  | TME       | XDL     |

Raw count matrices (`.csv`) are private and available upon request. Place them in a `data/` directory.

## Project Structure

```
├── data/                        # Raw count matrices (not included)
├── outputs/                     # Saved Seurat objects (.rds)
├── Raw_Data_Import.Rmd          # Step 1: Import and create Seurat object
├── QC_to_Clustering_Seurat.Rmd  # Step 2: QC, normalization, clustering
├── Seurat_DifferentialExpression_Classification.Rmd  # Step 3: DE, annotation, composition
└── README.md
```

## Pipeline

### 1. Raw Data Import (`Raw_Data_Import.Rmd`)

- Reads `.csv` count matrices and transposes them into gene × cell format
- Creates individual Seurat objects per sample (`min.cells = 3`, `min.features = 200`)
- Merges into a single object with cell ID prefixes
- Adds `condition` (T/TME) and `patient_id` metadata columns
- Saves: `outputs/merged_Seurat_obj.rds`

### 2. QC, Normalization, and Clustering (`QC_to_Clustering_Seurat.Rmd`)

- Calculates percent mitochondrial (`percent.mt`) using `^MT\\.` pattern
- Visualizes QC metrics (nCount, nFeature, percent.mt) per sample
- Applies filtering: `nCount_RNA > 500` and `percent.mt < 30`
- Normalizes with **SCTransform** (regressing out `percent.mt`)
- Runs PCA (30 PCs for neighbors), FindNeighbors, FindClusters (resolution 0.1), and UMAP (30 dims)
- Saves: `outputs/obj_filt_sctran_clust0.1.rds`

### 3. Differential Expression, Annotation, and Composition (`Seurat_DifferentialExpression_Classification.Rmd`)

**Differential Expression:**
- `FindAllMarkers` (Wilcoxon) across clusters to identify cluster marker genes
- `FindMarkers` (Wilcoxon) between T vs. TME conditions
- Pseudobulk DE via `AggregateExpression` + DESeq2 to address p-value inflation
- Visualization with DoHeatmap, DotPlot, VlnPlot, and FeaturePlot

**Cell Type Annotation (SingleR):**
- Reference: `celldex::HumanPrimaryCellAtlasData()`
- Cell-level annotation → `humanRNASeq.main`
- Cluster-level annotation via `AverageExpression` → `humanRNASeq.main.clust`

**Differential Composition:**
- Computes relative cell type proportions per sample
- Compares composition between T vs. TME conditions
- Compares composition across patients (JMY, WHG, XDL)
- Saves: `outputs/obj_final_Seurat.rds`

## Key Dependencies

| Package      | Source        | Purpose                          |
|--------------|---------------|----------------------------------|
| Seurat       | CRAN          | Core scRNA-seq toolkit           |
| tidyverse    | CRAN          | Data wrangling and plotting      |
| glmGamPoi    | Bioconductor  | SCTransform speed improvement    |
| presto       | GitHub        | Fast Wilcoxon DE                 |
| SingleR      | Bioconductor  | Reference-based cell annotation  |
| celldex      | Bioconductor  | Annotation reference datasets    |
| MAST         | Bioconductor  | Hurdle model DE                  |
| patchwork    | CRAN          | Composite plot layouts           |
| hdf5r        | CRAN          | HDF5 data import                 |

## Notes

- The `XDL_Ca` sample shows a different count distribution from other tumor samples; sample-specific QC thresholds are discussed but a universal threshold is applied.
- Memory allocation may need adjustment for `PrepSCTFindMarkers` — use `mem.maxVSize(vsize=16384*4)` if needed.
- The integrated/batch-corrected assay should **not** be used for differential expression; use the `SCT` assay instead.
