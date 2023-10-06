"# Master's Thesis Poject" 
# Developing methods and tools to enrich regions of gene of interest on the UMAP for single and multipe genes for a single cell RNA-seq data. Linking the enrichment to pathway information which will in turn help us biologically interpret the genes and cell types in the Single Cell data-set.
# Tools can be a good modification or integration to the Seurat pipeline used for analysis of Single-cell Data.


Tools used:
# Kernel Smoothing: KDE estimation #
Package: ks
Version: 1.14.0
Date: 2022-11-24
Title: Kernel Smoothing
Author: Tarn Duong [aut, cre] (<https://orcid.org/0000-0002-1198-3482>),
  Matt Wand [ctb] (<https://orcid.org/0000-0003-2555-896X>),
  Jose Chacon [ctb],
  Artur Gramacki [ctb] (<https://orcid.org/0000-0002-1610-9743>)
Repository: CRAN

# Spatial Statistics : Spatial Autocorrelation tests # 
Local Moran's I statistics
https://cran.r-project.org/web/packages/spdep/spdep.pdf

# Pathway Infomation aka scores for each cell : AUCell #
https://github.com/aertslab/AUCell
