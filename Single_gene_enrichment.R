# loading all the libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(ks)
library(spdep)
library(wesanderson)
library(AUCell)
library(GSEABase)
library(devEMF)
# loading all PBMC data set 

pbmc.data = Read10X(data.dir = "/path_to_to_the_file/filtered_matrices_genes_bc")


pbmc = CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

# QC step
# filtering out cells which have unique feature counts over 2,500 or less than 200
# filtering cells that have >5% mitochondrial counts


pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc = subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt <5)

# normalizing data 

pbmc = NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# feature selection

pbmc = FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Scaling data
all.genes <- rownames(pbmc)
pbmc = ScaleData(pbmc, features = all.genes)

# performing PCA on scaled data

pbmc = RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Finding neighbors

pbmc = FindNeighbors(pbmc,return.neighbor = TRUE)


# Cluster Cells

pbmc = FindNeighbors(pbmc, dims = 1:10)
pbmc = FindClusters(pbmc, resolution = 0.5)

# running UMAP

pbmc = RunUMAP(pbmc, dims = 1:10)


# Finding Variable features 
pbmc = FindVariableFeatures(pbmc, selection.method = "disp", nfeatures = 2000)



# find markers for every cluster compared to all remaining cells, report only the positive
# ones

pbmc.markers = FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

# Making a custom theme fot the UMAP plots

custom_theme <- function(base_size = 12, base_family = "") {
  theme_minimal() +
    theme(
      text = element_text(family = base_family),
      axis.title = element_text(size = rel(1.2)),
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.title.y = element_text(margin = margin(r = 10)),
      plot.title = element_text(size = rel(1.5), face = "bold", hjust = 0.5),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      panel.grid.major = element_line(color = "gray80"),  # Add grid lines
      panel.grid.minor = element_blank(),
      strip.text = element_text(size = rel(1.1)),
      legend.title = element_text(size = rel(1.1)),
      legend.text = element_text(size = rel(1)),
      legend.key = element_blank(),
      legend.background = element_rect(fill = "white", color = "black"),
      legend.position = "right",
      legend.spacing.y = unit(0.2, "cm"),
      legend.box.spacing = unit(0.2, "cm"),
      plot.margin = margin(1, 1, 1, 1, "cm")
    )
}


##### function to create a scatter plot for a single gene having colored by density values using ggplot2 #####

single_gene_enrichment = function(gene_s,data,densityplot,LISA_logfc,LISA_log10p)
{
  # extracting gene expression data
  expression_data = FetchData(data,vars = gene_s)
  
  # extracting cell embedding from the Seurat object
  umap_coords = data@reductions$umap@cell.embeddings
  
  # creating data frame with expression data and cluster identities
  df = data.frame(Expression = expression_data,UMAP1 = umap_coords[, 1], UMAP2 = umap_coords[, 2])
  
  # changing name of the column to maintain the uniformity
  names(df)[names(df) == gene_s] = "Expression"
  
  # Calculating density values using KDE function
  density_values = kde(x=umap_coords,w=df$Expression,compute.cont = TRUE,binned = FALSE,eval.points = umap_coords)
  
  # Adding the density estimates to the main data-frame
  df$density = density_values$estimate
  
  # Creating a neighbors list using UMAP coordinates from k-nn object
  nb_list = knn2nb(knearneigh(umap_coords,k=20))
  
  # converting the neighbors list into spatial weights matrix
  w = nb2listw(nb_list,style = "W")
  
  # Local Moran's I statistics calculation to signify spatial enrichment 
  LISA_statistics = localmoran(as.vector(df$density),listw = w,conditional = TRUE,alternative = "two.sided",spChk = NULL,adjust.x = FALSE)
  
  # converting the stats results into a data frame
  LISA_statistics = as.data.frame(LISA_statistics)
  
  # Calculating -log10p and log2fold changes using the Moran's I statistical test values
  df$mlog10p = -log10(LISA_statistics$`Pr(z != E(Ii))`)
  df$logfc = log2(LISA_statistics$Z.Ii)

  
  if (densityplot==TRUE)
  {
    # plotting the UMAP using for density Visualizations
    
    # Color palette for the UMAP plots
    pal = pal = wes_palette("Rushmore1", 100, type = "continuous")
    
    plot1 = ggplot(df, aes(x = UMAP1, y = UMAP2,color=df$density)) +
      geom_point(size=1.5)+
      labs(x = "UMAP 1", y = "UMAP 2", title =gene_s) +
      custom_theme()+
      scale_color_gradientn(name = "Density", colors = pal)
    return(plot1)
  }
  else if(LISA_logfc==TRUE)
  {
    # Color palette for the UMAP plots
    pal = wes_palette("Rushmore1", 100, type = "continuous")
    
    plot2 = ggplot(df, aes(x = UMAP1, y = UMAP2)) +
      geom_point(data = df,aes(colour = df$logfc),size = 1)+
      labs(x = "UMAP 1", y = "UMAP 2", title = gene_s) +
      theme_minimal()+
      scale_color_gradientn(name = "Logfc:(log(z-scores))", colors = pal)
    return(plot2)
  }
  
  else if (LISA_log10p == TRUE)
    
    # Color palette for the UMAP plots
    pal = wes_palette("Rushmore1", 100, type = "continuous")
  
  plot3 = ggplot(df, aes(x = UMAP1, y = UMAP2)) +
    geom_point(data = df,aes(colour = df$mlog10p),size = 1)+
    labs(x = "UMAP 1", y = "UMAP 2", title =gene_s)+
    custom_theme()+
    scale_color_gradientn(name = "-log10p", colors = pal)
  return(plot3)
  
}

d_plot = single_gene_enrichment("MS4A1",pbmc,densityplot = FALSE,LISA_logfc =FALSE,LISA_log10p=TRUE)
print(d_plot)


# function for saving plots to disk 

saving_plots = function(path,ggp)
{
  emf(path,height=10,width=14)
  print(ggp)
  dev.off()
}

saving_plots("/path_to_to_the_file/Enrichment_MS4A1.emf",d_plot)

