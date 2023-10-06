"# Master's Thesis Poject" 
# Developing methods and tools to enrich regions of gene of interest on the UMAP for single and multipe genes for a single cell RNA-seq data. Linking the enrichment to pathway information which will in turn help us biologically interpret the genes and cell types in the data-set.
# Tools can be a good modification or integration to the Seurat pipeline used for analysis of Single-cell Data.
Single-cell RNA sequencing (sc-RNA-seq) is a powerful technique that allows researchers to study the gene expression profiles of individual cells. This has led to a better understanding of the heterogeneity of cell populations and the complex biological processes that take place within various types of tissues. The Seurat vignette provides tools to plot sc-RNA-seq data on to the UMAP (Uniform Manifold Approximation and projection) and interpret gene expressions across a population cells. However, the existing tools do not provide a method to visualize and interpret single or multiple genes enrichment (according to a particular score or threshold) in the UMAP space. Additionally, they do not link or show relevance to the pathway information active across the cells in the dataset. This study hypothesizes that visualization and interpretation of KDE (Kernel Density Estimation) estimates coupled with spatial statistic tool i.e., local Moran’s I statistics will result in single and multiple gene enrichments on the UMAP space. Additionally, visualizations backed with pathway information via scores resulted from AUCell tool can be used to interpret pathways active in cell and sub-cell types in the single cell data. These methods were performed on a publicly available dataset. The results of this study have proven that the proposed methods have a potential to provide single, multiple gene enrichment and pathway information for the single cell data. This could be a valuable tool for researchers to better understand the biology and cause of diseases in various organisms.




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

# Pathway Infomation aka scores for each cell on the UMAP : AUCell #
https://github.com/aertslab/AUCell
