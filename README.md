# Loupy
Export clusters, annotations and projections from ScanPy's AnnData object to Loupe Cell Browser's .cloupe file format.


Uses LoupeR from 10X: [https://github.com/10XGenomics/loupeR](https://github.com/10XGenomics/loupeR)

# Usage

## Export all projections and categorical annotations
 - `python loupy/loupy.py pbmc3k_analysed.h5ad pbmc3k.cloupe --layer=counts`

## Export specific projections and categorical annotations
- `python loupy/loupy.py pbmc3k_analysed.h5ad pbmc3k.cloupe --layer=counts --projections=tsne,umap --categoricals=cell_type`

# Environment

## Python
    - ScanPy
    - rpy2

## R
    - Seurat
    - SeuratData
    - SeuratDisk
    - loupeR
    - rhdf5
