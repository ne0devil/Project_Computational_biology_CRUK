
# loading all the libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(ks)
library(spdep)
library(AUCell)
library(GSEABase)
library(wesanderson)
library(devEMF)
# loading all PBMC data set 

pbmc.data = Read10X(data.dir = "D:/AbroadStudies/Uofg/Project/filtered_matrices_genes_bc")


pbmc = CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc


# QC step
# filtering out cells which have unique feature counts over 2,500 or less than 200
# filtering cells that have >5% mitochondrial counts


pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc = subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt <5)

# normalizing data 

pbmc = NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc = NormalizeData(pbmc)

# feature selection

pbmc = FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Scaling data 

all.genes = rownames(pbmc)
pbmc = ScaleData(pbmc, features = all.genes)

# performing PCA on scaled data

pbmc = RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Finding neighbors

pbmc = FindNeighbors(pbmc,features = VariableFeatures(object = pbmc),return.neighbor = TRUE)


# Cluster Cells

pbmc = FindNeighbors(pbmc, dims = 1:10)
pbmc = FindClusters(pbmc, resolution = 0.5)

# running UMAP

pbmc = RunUMAP(pbmc, dims = 1:10)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones

pbmc.markers = FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)


# Making a custom theme for the UMAP plots

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

# loading the RDS object
# pbmc = readRDS(file = "pbmc_seuratobject.rds")

# extracting cell embedding from the Seurat object
umap_coords = pbmc@reductions$umap@cell.embeddings
View(umap_coords)

##### function to extract the expression data and add it to a data-frame #####


expression_profiles = function(gene_s,data){
  
  # extracting gene expression data
  expression_data = FetchData(data,vars = gene_s)
  
  # creating data frame with expression data and cluster identities
  df = as.data.frame(expression_data)
  return(df)
  
}

expression_df = expression_profiles(c("CD8A","IL7R"),pbmc)


##### function for calculation of KDE estimates and appending it to a data-frame ######


KDE_estimation = function(df,umapcoords)
{
  #making a variable that takes in list of values
  density_values = list()
  # iterating through the columns to calculate the densities for the expression values 
  for (i in 1:ncol(df)) {
    colname = paste("density", i, sep = "")
    density_values[[i]] = kde(x = umapcoords, w = df[,i], compute.cont = TRUE, binned = FALSE, eval.points = umapcoords)
    
    # Renaming the column names with respect to gene no.s
    df[[colname]] = density_values[[i]]$estimate
  }
  
  return(df)
}

global_df= KDE_estimation(expression_df,umap_coords)

# Extracting Column names which have densities in it
density_column_names = grep("density", colnames(global_df), value = TRUE)

# Product of all the densities to get probabilities of specific gene expression independently.  
product_of_densities = apply(global_df[,density_column_names],1,prod)

global_df$product_d = product_of_densities




####### Function to plot the estimated densities onto the UMAP #######

density_plot = function(df, umaap_coords) {
  
  #adding UMAP coordinates to the data-frame
  df$UMAP1 = umaap_coords[,1]
  df$UMAP2 = umaap_coords[,2]
  
  # Getting Gene names
  density_gene_names = colnames(expression_df)
  
  # plotting the UMAP
  # Color palette for the UMAP plots
  pal = wes_palette("Rushmore1", 100, type = "continuous")
  
  plot = ggplot(df, aes(x = UMAP1, y = UMAP2,color=df$product_d)) +
    geom_point(size=1.5)+
    labs(x = "UMAP 1", y = "UMAP 2", title =paste("Joint densities:",paste(density_gene_names, collapse = ", "))) +
    custom_theme()+
    scale_color_gradientn(name = "Density", colors = pal)
  return(plot)
  
  
}

jd_plot = density_plot(global_df,umap_coords)
print(jd_plot)


####### Function for spatial statistical test #######

LISA_statistics = function(df,umaap_coords){
  
  # Creating a neighbors list using UMAP coordinates from k-nn object
  nb_list = knn2nb(knearneigh(umaap_coords))
  
  # converting the neighbors list into spatial weights matrix
  w = nb2listw(nb_list,style = "W",glist = df$product_d)
  
  # Local Moran's I statistics calculation to signify spatial enrichment 
  statistics = localmoran(as.vector(df$product_d),listw = w,conditional = TRUE,alternative = "two.sided",spChk = NULL,adjust.x = FALSE)
  
  # converting the stats results into a data frame
  statistics = as.data.frame(statistics)
  
 return(statistics)
  
}

LISA_statistics = LISA_statistics(global_df,umap_coords)
print(LISA_statistics)


####### Function to plot enrichment plots using statistical inferences #######

plot_with_statistics = function(df,umaap_coords,statistics,LISA_logfc,LISA_log10p){
  
  #adding UMAP coordinates to the data-frame
  df$UMAP1 = umaap_coords[,1]
  df$UMAP2 = umaap_coords[,2]
  
  # Calculating -log10p and log2fold changes using the Moran's I statistical test values
  df$mlog10p = -log10(statistics$`Pr(z != E(Ii))`)
  df$logfc = log2(statistics$Z.Ii)
  
  
  # Plotting the UMAP using ggplot2
  # Getting Gene names
  density_gene_names = colnames(expression_df)

  
  if(LISA_logfc==TRUE)
  {
    # Color palette for the UMAP plots
    pal = wes_palette("Rushmore1", 100, type = "continuous")
    
    plot1 = ggplot(df, aes(x = UMAP1, y = UMAP2)) +
      geom_point(data = df,aes(colour = df$logfc),size = 1)+
      labs(x = "UMAP 1", y = "UMAP 2", title =density_gene_names)+
      custom_theme()+
      scale_color_gradientn(name = "logfc", colors=pal)
    return(plot1)
  }
  
  
  else if(LISA_log10p==TRUE)
  {
    # Color palette for the UMAP plots
    pal = wes_palette("Rushmore1", 100, type = "continuous")
    
    plot2 = ggplot(df, aes(x = UMAP1, y = UMAP2)) +
      geom_point(data = df,aes(colour = df$mlog10p),size = 1)+
      # contour(df$mlog10p,col=pal)+
      labs(x = "UMAP 1", y = "UMAP 2", title =paste("Joint Enrichment:",paste(density_gene_names, collapse = ", ")))+
      custom_theme()+
      scale_color_gradientn(name = "-log10p",colors=pal)
    return(plot2)
    
  }
  
  }

jd_plot_with_statistics = plot_with_statistics(global_df,umap_coords,LISA_statistics,LISA_logfc=FALSE,LISA_log10p =TRUE)
print(jd_plot_with_statistics)

####### Function to extract pathway information and plot in UMAP space #######

# Constructing a gene Set
set.seed(123)
broad_set1 = getBroadSets(uri = c("D:/AbroadStudies/Uofg/Project/Pathways/BIOCARTA_BLYMPHOCYTE_PATHWAY.v2023.1.Hs.xml","D:/AbroadStudies/Uofg/Project/Pathways/BIOCARTA_BCR_PATHWAY.v2023.1.Hs.xml","D:/AbroadStudies/Uofg/Project/Pathways/HAY_BONE_MARROW_NAIVE_T_CELL.v2023.1.Hs.xml","D:/AbroadStudies/Uofg/Project/Pathways/HAY_BONE_MARROW_NK_CELLS.v2023.1.Hs.xml","D:/AbroadStudies/Uofg/Project/Pathways/HAY_BONE_MARROW_CD8_T_CELL.v2023.1.Hs.xml","D:/AbroadStudies/Uofg/Project/Pathways/HAY_BONE_MARROW_MONOCYTE.v2023.1.Hs.xml","D:/AbroadStudies/Uofg/Project/Pathways/HAY_BONE_MARROW_PLATELET.v2023.1.Hs.xml","D:/AbroadStudies/Uofg/Project/Pathways/HAY_BONE_MARROW_FOLLICULAR_B_CELL.v2023.1.Hs.xml","D:/AbroadStudies/Uofg/Project/Pathways/HAY_BONE_MARROW_IMMATURE_NEUTROPHIL.v2023.1.Hs.xml"),membersId ="MEMBERS_SYMBOLIZED")

# extracting expression matrix from the Seurat object
expression_mat = pbmc@assays[["RNA"]]@counts
# ranking the genes in the data-set
cells_rankings = AUCell_buildRankings(expression_mat, plotStats=FALSE)

# Calculating AUC activity levels for the gene-set in each of the cell
cells_AUC = AUCell_calcAUC(broad_set1, cells_rankings,normAUC = TRUE)

# Extracting AU-cell scores to a data-frame
aucell_scores = as.data.frame(as.matrix(cells_AUC@assays@data@listData[["AUC"]]))
aucell_scores = data.frame(t(aucell_scores))
View(aucell_scores)
pathway_enrichment=function(scoresdf,df,pathway_name,umaap_coords)
{
  
  #adding UMAP coordinates to the data-frame
  df$UMAP1 = umaap_coords[,1]
  df$UMAP2 = umaap_coords[,2]
  
  # Extracting column names for specific pathway name
  pathway_column_names = grep(pathway_name, colnames(scoresdf), value = TRUE)
  
  #Adding the pathway information into the main data frame
  df$pathway = scoresdf[,pathway_column_names]
  
  # Creating a neighbors list using UMAP coordinates from k-nn object
  nb_list = knn2nb(knearneigh(umaap_coords))
  
  # converting the neighbors list into spatial weights matrix
  w = nb2listw(nb_list,style = "W",glist = df$pathway)
  
  # Local Moran's I statistics calculation to signify spatial enrichment 
  statistics = localmoran(as.vector(df$pathway),listw = w,conditional = TRUE,alternative = "two.sided",spChk = NULL,adjust.x = FALSE)
  
  # converting the stats results into a data frame
  statistics = as.data.frame(statistics)
  
  # Calculating -log10p using the Moran's I statistical test values
  df$mlog10p_pathway = -log10(statistics$`Pr(z != E(Ii))`)
  
  # Color palette for the UMAP plots
  pal = wes_palette("Rushmore1", 100, type = "continuous")
  
  plot2 = ggplot(df, aes(x = UMAP1, y = UMAP2)) +
    geom_point(data = df,aes(colour = df$mlog10p_pathway),size = 1)+
    labs(x = "UMAP 1", y = "UMAP 2", title =paste("Pathway_enrichment:",pathway_name))+
    custom_theme()+
    scale_color_gradientn(name = "-log10p", colors = pal)
  return(plot2)
}

pe_plot = pathway_enrichment(scoresdf =aucell_scores,df = global_df,"NEUTROPHIL",umaap_coords = umap_coords)
print(pe_plot)




####### Function to extract multiple pathway information and plot in UMAP space #######

multiple_pathway_enrichment=function(scoresdf,df,pathway_name,umaap_coords)
{
  
  #adding UMAP coordinates to the data-frame
  df$UMAP1 = umaap_coords[,1]
  df$UMAP2 = umaap_coords[,2]
  
  # Extracting column names for specific pathway name
  pathway_column_names = list()
  for (i in 1:length(pathway_name)){
    pathway_column_names[[i]] = grep(pathway_name[i], colnames(scoresdf), value = TRUE)
  }
  pathway_column_names = unlist(pathway_column_names)
  
  # Product of the pathway scores to get probabilities of specific pathway enrichment independently.  
  product_of_pathways = apply(scoresdf[,pathway_column_names],1,prod)
  
  #Adding the product of pathway information into the main data frame
  df$pathway = product_of_pathways
  
  # Creating a neighbors list using UMAP coordinates from k-nn object
  nb_list = knn2nb(knearneigh(umaap_coords))
  
  # converting the neighbors list into spatial weights matrix
  w = nb2listw(nb_list,style = "W",glist = df$pathway)
  
  # Local Moran's I statistics calculation to signify spatial enrichment 
  statistics = localmoran(as.vector(df$pathway),listw = w,conditional = TRUE,alternative = "two.sided",spChk = NULL,adjust.x = FALSE)
  
  # converting the stats results into a data frame
  statistics = as.data.frame(statistics)
  
  # Calculating -log10p using the Moran's I statistical test values
  df$mlog10p_pathway = -log10(statistics$`Pr(z != E(Ii))`)
  
  # Color palette for the UMAP plots
  pal = wes_palette("Rushmore1", 100, type = "continuous")
  
  plot2 = ggplot(df, aes(x = UMAP1, y = UMAP2)) +
    geom_point(data = df,aes(colour = df$mlog10p_pathway),size = 1)+
    labs(x = "UMAP 1", y = "UMAP 2", title =paste("Multiple_pathway_enrichment:",paste(pathway_name, collapse = ", ")))+
    custom_theme()+
    scale_color_gradientn(name = "-log10p", colors = pal)
  return(plot2)
}

jp_pe_plot = multiple_pathway_enrichment(scoresdf =aucell_scores,df = global_df,pathway_name = c("NEUTROPHIL","MONOCYTE"),umaap_coords = umap_coords)
print(jp_pe_plot)



 
# function for saving plots to disk 

saving_plots = function(path,ggp)
{
  emf(path,height=10,width=14)
  print(ggp)
  dev.off()
}

saving_plots("D:/AbroadStudies/Uofg/Project/plots/Pathway(multiple_Neutrophil_monnocyte)_enrichment.emf",jp_pe_plot)


