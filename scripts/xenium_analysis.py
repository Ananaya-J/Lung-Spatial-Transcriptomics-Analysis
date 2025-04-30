import scanpy as sc
import pandas as pd
import numpy as np
import squidpy as sq
import matplotlib.pyplot as plt
import h5py
import warnings
import os

# Suppress runtime and performance warnings for cleaner output
warnings.filterwarnings('ignore')

# --- Load gene expression matrix ---
adata = sc.read_10x_h5("data/raw/xenium_lung_ffpe/cell_feature_matrix.h5")
print("AnnData shape after loading:", adata.shape)

# --- Load cell metadata and align ---
cell_meta = pd.read_csv("data/raw/xenium_lung_ffpe/cells.csv.gz", index_col=0)
print("cell_meta shape:", cell_meta.shape)
print("Number of matching barcodes:", sum(cell_meta.index.isin(adata.obs_names)))
adata.obs = cell_meta.reindex(adata.obs_names)
print("Number of cells with all-NaN metadata:", adata.obs.isnull().all(axis=1).sum())

# --- Check AnnData shape before filtering ---
print("AnnData shape before filtering:", adata.shape)

# Minimal filtering (or skip if already QC'd)
sc.pp.filter_cells(adata, min_genes=1)
sc.pp.filter_genes(adata, min_cells=1)
print("AnnData shape after filtering:", adata.shape)

# --- Check matrix stats before normalization ---
try:
    X = adata.X.A if hasattr(adata.X, 'A') else adata.X  # handle sparse
    print("Before normalization: min =", np.min(X), "max =", np.max(X), "NaNs:", np.isnan(X).sum())
except Exception as e:
    print("Could not compute matrix stats before normalization:", e)

if adata.shape[0] > 0 and adata.shape[1] > 0:
    # Normalize and log-transform
    sc.pp.normalize_total(adata)
    try:
        X = adata.X.A if hasattr(adata.X, 'A') else adata.X
        print("After normalization: min =", np.min(X), "max =", np.max(X), "NaNs:", np.isnan(X).sum())
    except Exception as e:
        print("Could not compute matrix stats after normalization:", e)
    sc.pp.log1p(adata)
    try:
        X = adata.X.A if hasattr(adata.X, 'A') else adata.X
        print("After log1p: min =", np.min(X), "max =", np.max(X), "NaNs:", np.isnan(X).sum())
    except Exception as e:
        print("Could not compute matrix stats after log1p:", e)

    if np.max(X) > 0 and np.isnan(X).sum() == 0:
        # Find highly variable genes
        sc.pp.highly_variable_genes(adata)
        # Scale and run PCA
        sc.pp.scale(adata)
        sc.tl.pca(adata)
        sc.pp.neighbors(adata)
        sc.tl.leiden(adata)
        sc.tl.umap(adata)
        sc.pl.umap(adata, color=["leiden"], show=True, save="_xenium_umap.png")

        # Spatial plot
        if {'x_centroid', 'y_centroid'}.issubset(adata.obs.columns):
            plt.figure(figsize=(6, 6))
            plt.scatter(
                adata.obs['x_centroid'],
                adata.obs['y_centroid'],
                c=adata.obs['leiden'].astype(int),
                cmap='tab20',
                s=1
            )
            plt.xlabel('x_centroid')
            plt.ylabel('y_centroid')
            plt.title('Spatial plot colored by Leiden clusters')
            plt.gca().set_aspect('equal')
            plt.show()

        # Filter genes before marker analysis
        sc.pp.filter_genes(adata, min_cells=3)
        sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
        sc.pl.rank_genes_groups(adata, n_genes=10, sharey=False, show=True, save="_xenium_marker_genes.png")

        # --- 1. Spatial Neighborhood Analysis (Squidpy) ---
        if {'x_centroid', 'y_centroid'}.issubset(adata.obs.columns):
            adata.obsm['spatial'] = adata.obs[['x_centroid', 'y_centroid']].values
            sq.gr.spatial_neighbors(adata, coord_type='generic')
            sq.gr.spatial_autocorr(adata, mode='moran')

            output_dir = os.path.abspath(os.path.join(os.getcwd(), 'bulk_data'))
            os.makedirs(output_dir, exist_ok=True)
            print("Saving outputs to:", output_dir)

            plt.figure(figsize=(6, 6))
            plt.scatter(
                adata.obs['x_centroid'],
                adata.obs['y_centroid'],
                c=adata.obs['leiden'].astype(int),
                cmap='tab20', s=1
            )
            plt.xlabel('x_centroid')
            plt.ylabel('y_centroid')
            plt.title('Spatial plot colored by Leiden clusters')
            plt.gca().set_aspect('equal')
            plt.savefig(f"{output_dir}/spatial_clusters.png")
            plt.show()

            # Moran's I visualization
            if 'moranI' in adata.uns:
                moran_df = adata.uns['moranI']
                top_moran = moran_df.sort_values('I', ascending=False).head(10)
                plt.figure(figsize=(8, 4))
                plt.bar(top_moran.index, top_moran['I'])
                plt.ylabel("Moran's I")
                plt.xlabel("Gene")
                plt.title("Top 10 spatially variable genes (Moran's I)")
                plt.xticks(rotation=45, ha='right')
                plt.tight_layout()
                plt.savefig(f"{output_dir}/top10_moransI_genes.png")
                plt.show()

        # --- 2a. Region Annotation ---
        if 'region' not in adata.obs.columns:
            if {'x_centroid', 'y_centroid'}.issubset(adata.obs.columns):
                adata.obs['region'] = np.where(adata.obs['x_centroid'] > adata.obs['x_centroid'].median(), 'Tumor', 'Stroma')
                print('Mock region annotation added (Tumor vs. Stroma)')

        # --- 2b. Region-Specific DE ---
        if 'region' in adata.obs.columns:
            sc.tl.rank_genes_groups(adata, 'region', method='t-test')
            sc.pl.rank_genes_groups(adata, n_genes=10, sharey=False, show=True, save="_region_marker_genes.png")
            region_de = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(10)
            region_de.to_csv(f"{output_dir}/region_marker_genes.csv")
            print('Region-specific DE analysis completed and saved.')

        # --- 2c. Neighborhood Enrichment ---
        if 'leiden' in adata.obs.columns:
            print('Running Squidpy neighborhood enrichment...')
            sq.gr.nhood_enrichment(adata, cluster_key='leiden')
            sq.pl.nhood_enrichment(adata, cluster_key='leiden', show=True, save="_nhood_enrichment.png")

        # --- 2d. Cell-Cell Communication (Placeholder) ---
        print('To perform cell-cell communication inference, export cluster assignments and expression matrix for use with CellPhoneDB or NicheNet.')
        adata.obs[['leiden']].to_csv(f"{output_dir}/cell_clusters_for_communication.csv")
        adata.to_df().to_csv(f"{output_dir}/expression_matrix_for_communication.csv")
        print('Exported files for ligand-receptor analysis.')

        # --- LTBP2 Visualization ---
        if 'LTBP2' in adata.var_names:
            print('Plotting UMAP and spatial heatmap for LTBP2')
            sc.pl.umap(adata, color='LTBP2', show=True, save="_LTBP2_umap.png")
            if {'x_centroid', 'y_centroid'}.issubset(adata.obs.columns):
                plt.figure(figsize=(6, 6))
                plt.scatter(
                    adata.obs['x_centroid'],
                    adata.obs['y_centroid'],
                    c=adata[:, 'LTBP2'].X.toarray().flatten() if hasattr(adata[:, 'LTBP2'].X, 'toarray') else adata[:, 'LTBP2'].X.flatten(),
                    cmap='hot', s=3
                )
                plt.title(f"Spatial heatmap: LTBP2")
                plt.xlabel('x_centroid')
                plt.ylabel('y_centroid')
                plt.colorbar(label='Expression')
                plt.savefig(f"{output_dir}/spatial_heatmap_LTBP2.png")
                plt.show()
        else:
            print('LTBP2 not found in gene list.')

        # --- 3. Differential Expression Between Clusters ---
        import collections
        if isinstance(adata.uns['rank_genes_groups']['names'], np.ndarray):
            de_df = pd.DataFrame(adata.uns['rank_genes_groups']['names'])
        else:
            de_df = pd.DataFrame({
                group: adata.uns['rank_genes_groups']['names'][group]
                for group in adata.uns['rank_genes_groups']['names'].dtype.names
            })
        de_df.to_csv(f"{output_dir}/differential_expression_by_cluster.csv")

        # --- 4. Custom Spatial Heatmap for Top Marker ---
        try:
            top_gene = adata.uns['rank_genes_groups']['names']['0'][0]
        except Exception:
            top_gene = None

        if top_gene and top_gene in adata.var_names and {'x_centroid', 'y_centroid'}.issubset(adata.obs.columns):
            plt.figure(figsize=(6, 6))
            plt.scatter(
                adata.obs['x_centroid'],
                adata.obs['y_centroid'],
                c=adata[:, top_gene].X.toarray().flatten() if hasattr(adata[:, top_gene].X, 'toarray') else adata[:, top_gene].X.flatten(),
                cmap='hot', s=3
            )
            plt.title(f"Spatial heatmap: {top_gene}")
            plt.xlabel('x_centroid')
            plt.ylabel('y_centroid')
            plt.colorbar(label='Expression')
            plt.savefig(f"{output_dir}/spatial_heatmap_{top_gene}.png")
            plt.show()
    else:
        print("Matrix is all zeros or contains NaN after normalization/log1p. Cannot proceed with downstream analysis.")
else:
    print("AnnData object is empty after filtering. Please check your data or filtering thresholds.")

import scanpy as sc
import pandas as pd
import numpy as np
from sklearn.preprocessing import MinMaxScaler

# 1. Define lung cell type markers
gene_sets = {
    'Airway epithelial cells': ['ADH7', 'AQP1', 'CDH1', 'SEC14L3'],
    'Airway goblet cells': ['AGR2', 'AQP5', 'CEACAM1', 'DMBT1', 'DUSP4'],
    'Mesothelial cells': ['C2', 'CALB2', 'CD44', 'CDH1', 'DES'],
    'Fibroblasts': ['COL3A1', 'COL5A2', 'DPT', 'FN1', 'GSN'],
    'Basal cells (Airway progenitor cells)': ['ABI3BP', 'AQP3', 'DAPL1', 'GSTM2', 'HPGD'],
    'Alveolar macrophages': ['ABCG1', 'CCL3', 'CD36', 'CLEC7A', 'CSF2RB'],
    'Ciliated cells': ['APPL2', 'ATP5MD', 'CCDC153', 'CCDC17', 'CCDC181'],
    'Clara cells': ['AHR', 'ALDH1A1', 'BPIFA1', 'CTSE', 'CYP2E1'],
    'Immune system cells': ['ADGRE1', 'ARG1', 'BIRC5', 'CCL17', 'CCL18'],
    'Endothelial cell': ['CD34', 'EGFL7', 'EMCN', 'ESAM', 'FLT1'],
    'Epithelial cells': ['ANPEP', 'EPCAM', 'IL10', 'IL6R'],
    'Ionocytes': ['CFTR', 'CLCNKB', 'FOXI1', 'KCNMA1', 'SCGB1A1'],
    'Pulmonary alveolar type I cells': ['AGER', 'AKAP5', 'AQP3', 'AQP5', 'CCN2'],
    'Pulmonary alveolar type II cells': ['ABCA3', 'ADGRF5', 'AGER', 'CD36', 'CD3G'],
    'Secretory cell': ['MUC5B', 'PIGR', 'SCGB1A1'],
    'Cancer stem cells': ['ABCG2', 'ALCAM', 'ALDH1A1', 'BMI1', 'CD24']
}

# 2. Initialize scoring matrix (fixed parentheses)
scoring_matrix = pd.DataFrame(
    np.zeros((adata.n_obs, len(gene_sets))),
    index=adata.obs_names,
    columns=list(gene_sets.keys())
)

# 3. Calculate average expression for each cell type's markers
for ct, genes in gene_sets.items():
    valid_genes = [g for g in genes if g in adata.var_names]
    if valid_genes:
        scoring_matrix[ct] = adata[:, valid_genes].X.mean(axis=1)

# 4. Normalize scores (0-1 range)
scaler = MinMaxScaler()
scoring_matrix[:] = scaler.fit_transform(scoring_matrix)

# 5. Assign cell types with confidence
adata.obs['cell_type'] = scoring_matrix.idxmax(axis=1)
adata.obs['confidence'] = scoring_matrix.max(axis=1)

# 6. Compute UMAP if missing
if 'X_umap' not in adata.obsm:
    sc.pp.pca(adata)
    sc.pp.neighbors(adata, n_pcs=30, n_neighbors=15)
    sc.tl.umap(adata)

# 7. Visualize results with legend beside the plot (not on data)
sc.pl.umap(adata, 
           color='cell_type', 
           legend_loc='best',  # Changed from 'on data' to 'right'
           legend_fontsize=6,   # Adjust legend font size as needed
           frameon=False,
           title='Cell Type Annotation')

# 8. Visualize spatial distribution based on available coordinates
# Check which spatial coordinates are available
if 'spatial' in adata.uns:
    # Standard Visium/10X format
    sc.pl.spatial(adata, 
                  color='cell_type', 
                  spot_size=20,
                  title='Spatial Distribution')
elif {'x_centroid', 'y_centroid'}.issubset(adata.obs.columns):
    # Custom spatial coordinates in obs
    sc.pl.embedding(adata, 
                   basis='spatial',
                   color='cell_type',
                   gene_symbols=None,
                   use_raw=False,
                   title='Spatial Distribution')
elif {'x', 'y'}.issubset(adata.obs.columns):
    # Alternative coordinate naming
    sc.pl.embedding(adata, 
                   basis='spatial',
                   color='cell_type',
                   gene_symbols=None,
                   use_raw=False,
                   title='Spatial Distribution')
    
# If trying to use squidpy, handle it separately with proper error checking
def plot_with_squidpy():
    try:
        import squidpy as sq
        # Check if spatial data exists in the expected format
        if 'spatial' in adata.uns:
            sq.pl.spatial_scatter(adata, 
                                color='cell_type',
                                spot_size=20, 
                                title='Spatial Distribution')
        else:
            print("Cannot use squidpy.pl.spatial_scatter: 'spatial' key not found in adata.uns")
            # Fall back to basic plotting of coordinates if available
            if {'x_centroid', 'y_centroid'}.issubset(adata.obs.columns):
                # Create spatial basis from centroid coordinates
                adata.obsm['spatial'] = adata.obs[['x_centroid', 'y_centroid']].values
                sc.pl.embedding(adata, 
                              basis='spatial',
                              color='cell_type',
                              title='Spatial Distribution')
    except ImportError:
        print("Squidpy not available. Using scanpy for visualization.")
