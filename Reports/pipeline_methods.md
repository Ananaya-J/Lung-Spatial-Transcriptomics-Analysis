# Xenium Spatial Transcriptomics Pipeline: Methods and Reproducibility

## Overview
This document describes the step-by-step methods used to process, analyze, and visualize spatial transcriptomics data from the 10x Genomics Xenium platform, using Python, Scanpy, Squidpy, and Matplotlib. The pipeline is designed for robust, reproducible analysis of custom spatial transcriptomics datasets.

---

## 1. Data Loading and Integrity Checks
- **Gene expression matrix** loaded from Xenium output (`cell_feature_matrix.h5`) using Scanpy.
- **Cell metadata** (`cells.csv.gz`) loaded and aligned to gene expression data.
- Data integrity verified by printing shapes and matching barcodes.

## 2. Quality Control and Filtering
- **Cells** filtered for minimal expression (≥1 gene).
- **Genes** filtered for minimal presence (≥1 cell).
- AnnData shape checked after filtering.

## 3. Normalization and Transformation
- **Total counts normalization** per cell.
- **Log1p transformation** for variance stabilization.
- Matrix statistics checked before and after normalization.

## 4. Feature Selection and Dimensionality Reduction
- **Highly variable genes** identified.
- **Scaling** and **Principal Component Analysis (PCA)** performed.

## 5. Clustering and Embedding
- **Neighbor graph** computed.
- **Leiden clustering** performed.
- **UMAP** embedding computed and plotted.

## 6. Spatial Visualization
- **Spatial coordinates** (`x_centroid`, `y_centroid`) checked in metadata.
- **Matplotlib** used for robust spatial scatter plots colored by cluster.

## 7. Marker Gene Analysis
- **Lowly expressed genes** filtered out (expressed in <3 cells).
- **Differential expression** via `rank_genes_groups` to identify cluster markers.
- Marker genes exported to CSV and plotted.

## 8. Advanced Spatial Analysis (Squidpy)
- **Spatial neighbors** and **spatial autocorrelation (Moran's I)** computed.
- **Top 10 spatially variable genes** visualized by Moran's I statistic.

## 9. Cell Type Annotation
- **Canonical lung marker genes** defined (AT1, AT2, Club, Ciliated, Endothelial, Macrophage).
- **UMAP and spatial heatmaps** plotted for each marker gene.
- **Special handling for LTBP2:** UMAP and spatial heatmap plotted if present.

## 10. Export and Output Organization
- All results (plots, CSVs, processed AnnData) saved in a dedicated `bulk_data` directory for reproducibility and easy access.

---

## Design Rationale
- **Matplotlib** is used for spatial plots to ensure compatibility with custom Xenium data.
- **Squidpy** is leveraged for spatial statistics, not for plotting (which requires Visium-style metadata).
- **Robust error handling** and **clear output organization** ensure reproducibility and easy troubleshooting.
- **Extensible marker gene logic** allows for easy customization.

---

## Reproducibility Notes
- All code is version-controlled and documented in `xenium_analysis.py`.
- The pipeline is modular and can be extended for further analyses (e.g., region-specific DE, integration with histology).
- Outputs are automatically organized by analysis step and saved with descriptive filenames.

---

## Dependencies
- Python ≥3.8
- scanpy
- squidpy
- matplotlib
- pandas
- numpy

---

## Citation
If you use this pipeline in your research, please cite:
- Wolf et al., "Scanpy: large-scale single-cell gene expression data analysis." Genome Biology (2018).
- Palla et al., "Squidpy: a scalable framework for spatial single cell analysis." Nature Methods (2022).
- 10x Genomics Xenium platform documentation.

---
