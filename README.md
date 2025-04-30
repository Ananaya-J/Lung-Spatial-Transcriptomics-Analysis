# 🔬 Xenium Lung Spatial Transcriptomics Analysis


## ✨ Project Overview

This repository contains a complete workflow for analyzing 10x Genomics Xenium spatial transcriptomics data, focusing on lung tissue. The pipeline processes raw Xenium data through preprocessing, quality control, clustering, spatial analysis, and cell type annotation based on lung-specific marker genes.

## 🧪 Features

- **Comprehensive Data Processing** | Raw data loading, QC, normalization, dimensionality reduction
- **Spatial Analysis** | Leverages Squidpy for neighborhood analysis and visualization
- **Cell Type Annotation** | Automated identification using lung-specific marker genes
- **Advanced Visualization** | UMAP embeddings, spatial plots, gene expression heatmaps
- **Differential Expression** | Analysis between clusters and spatial regions
- **Containerized Environment** | Docker setup for reproducible analysis

## 📁 Repository Structure

```
├── data/
│   └── raw/                    # Raw Xenium data
├── bulk_data/                  # Output images and results
├── scripts/
│   ├── xenium_analysis.py      # Main analysis script
│   └── utils/                  # Utility functions
├── docker/
│   └── Dockerfile.squidpy      # Python environment Dockerfile
├── docker-compose.yml          # Docker Compose configuration
├── reports/                    # Analysis reports
```

## 🔎 Key Findings

| Discovery | Description |
|-----------|-------------|
| **Cell Populations** | Identified distinct cell clusters in lung tissue |
| **Spatial Relationships** | Characterized neighborhood interactions between cell types |
| **Spatially Variable Genes** | Discovered genes like LTBP2 with spatial expression patterns |
| **Cell Type Atlas** | Annotated major lung cells (alveolar, fibroblasts, immune cells) |
| **Region Analysis** | Differential expression between tumor and stroma |

## 🛠️ Technologies Used

`Scanpy` • `Squidpy` • `Matplotlib` • `Docker` • `Python` • `Jupyter`

## 🚀 Getting Started

### Prerequisites

- Docker and Docker Compose
- Data(https://www.10xgenomics.com/datasets/preview-data-ffpe-human-lung-cancer-with-xenium-multimodal-cell-segmentation-1-standard)

### Quick Start

1. **Launch containers:**
   ```bash
   docker-compose up python-env
   ```

2. **Access Jupyter:**
   - Navigate to `http://localhost:8888`
   - Open `scripts/xenium_analysis.py` to run the analysis

## 📊 Results

The pipeline generates:

- Clustering results on UMAP embeddings
- Spatial plots of cell distributions
- Differential expression analyses
- Cell type annotations
- Neighborhood enrichment visualizations

## 🔮 Future Directions

- [ ] Integration with other Xenium datasets
- [ ] Cell-cell communication inference via CellPhoneDB
- [ ] Trajectory analysis for developmental insights
- [ ] Expansion to additional tissue types
- [ ] Deep learning for spatial pattern recognition

## 👏 Acknowledgments

Thanks to 10x Genomics for the Xenium platform, the Scanpy/Squidpy development teams, and the entire spatial transcriptomics community.

---

<p align="center">
  <i>Uncovering spatial biology, one cell at a time.</i>
</p>
