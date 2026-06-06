## Single cell RNAseq analyses using Seurat
## - Includes dimensionality reduction, clustering and cell type identification.
## Reva Shenwai
## 19-May-2024
## ---------------------------------------------------

# Set parameters
install_pkg <- FALSE    #do packages need to be installed?

# Install packages if required
if (install_pkg==TRUE) {
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("Seurat")
  install.packages("dplyr")
  install.packages("patchwork")
}

# Load libraries
library(dplyr)
library(Seurat)
library(patchwork)

## Load data
pbmc.data <- Seurat::Read10X(data.dir="~/Bioinformatics/R/data/filtered_gene_bc_matrices/hg19/")
# Initialize Seurat obj (raw data)
pbmc <- Seurat::CreateSeuratObject(counts=pbmc.data, project="Seurat_project",
                                   min.cells=3, min.features=200)
## QC
pbmc[["percent.mt"]] <- Seurat::PercentageFeatureSet(pbmc,pattern="^MT-")
# Visualize QC metrics
Seurat::VlnPlot(pbmc,features=c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3)
# Feature-feature relationships
plot1 <- FeatureScatter(pbmc, feature1="nCount_RNA",feature2="percent.mt")
plot2 <- FeatureScatter(pbmc, feature1="nCount_RNA",feature2="nFeature_RNA")
plot1 + plot2

## Normalize
pbmc <- Seurat::NormalizeData(pbmc)

## Analysis
# Feature selection
pbmc <- Seurat::FindVariableFeatures(pbmc,selection.method="vst",nfeatures=2000)
top10 <- head(Seurat::VariableFeatures(pbmc),10)  # top 10 DEGs
# Visualize
Seurat::LabelPoints(plot=Seurat::VariableFeaturePlot(pbmc),
                    points=top10,repel=TRUE)

# Scale
all.genes <- rownames(pbmc)
pbmc <- Seurat::ScaleData(pbmc, features=all.genes)

# Linear dimensional reduction
pbmc <- Seurat::RunPCA(pbmc,features=Seurat::VariableFeatures(pbmc))
# Visualize PCA results
print(pbmc[["pca"]],dims=1:5,nfeatures=5)
Seurat::VizDimLoadings(pbmc,dims=1:2,reduction="pca")
Seurat::DimPlot(pbmc,reduction="pca") + Seurat::NoLegend()
Seurat::DimHeatmap(pbmc,dims=1,cells=500,balanced=TRUE)
Seurat::DimHeatmap(pbmc, dims=1:15, cells=500, balanced=TRUE)

# Determine dimensionality of dataset
Seurat::ElbowPlot(pbmc)
pbmc <- Seurat::FindNeighbors(pbmc,dims=1:10)
pbmc <- Seurat::FindClusters(pbmc,resolution=0.5)

# Non-linear dimensional reduction
pbmc <- Seurat::RunUMAP(pbmc,dims=1:10)
Seurat::DimPlot(pbmc,reduction="umap")

# Identify cluster biomarkers
pbmc.markers <- Seurat::FindAllMarkers(pbmc,only.pos=TRUE)
# view
pbmc.markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::filter(avg_log2FC>1)

## Visualization of GEXP across clusters
# Useful for identifying cell type clusters
# Violin plot
Seurat::VlnPlot(pbmc, features=c("FCGR3A","MS4A7"))  # CD16+ monocyte markers
# Overlay gexp on UMAP
Seurat::FeaturePlot(pbmc,features=c("FCGR3A","MS4A7"))
# Heatmap showing top features per cluster
top10 <- pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC>1) %>%
  dplyr::slice_head(n = 10) %>%
  dplyr::ungroup()
Seurat::DoHeatmap(pbmc, features=top10$gene) + Seurat::NoLegend()

# Assign cell types to clusters based on biomarker genes
new.cluster.ids <- c("Naive CD4 T","CD14+ Mono","Memory CD4 T","B","CD8 T", 
                     "CD16+ Mono","NK","DC","Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- Seurat::RenameIdents(pbmc, new.cluster.ids)

# Plot UMAP with cell cluster assignments
Seurat::DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()









