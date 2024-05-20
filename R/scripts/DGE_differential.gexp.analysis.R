## Differential gene expression analysis on bulk RNAseq data
## Reva Shenwai
## 19-May-2024
## ---------------------------------------------------

# Set parameters
install_pkg <- FALSE    #do packages need to be installed?

# Install packages if required
if (install_pkg==TRUE) {
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("edgeR")
  install.packages("ggplot2")
  install.packages("ggrepel")
  install.packages("dplyr")
  BiocManager::install("ComplexHeatmap")
  install.packages("circlize")
}
  
# Load libraries
library(edgeR)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(ComplexHeatmap)
library(circlize)

## Input data
# sample gene expression matrix
gexp_mat <- read.csv("~/Bioinformatics/R/data/trt.vs.ctr_gexp.mat.csv",
                     row.names = 1)

# meta data
meta_data <- c("ctr", "ctr", "ctr", "trt", "trt", "trt")

# DGElist data object
dge <- edgeR::DGEList(counts=gexp_mat,
                      group=factor(meta_data))

## Analysis: single factor design
# Filter genes with low counts
keep <- edgeR:: filterByExpr(y=dge)
dge <- dge[keep, , keep.lib.sizes=FALSE]

# Normalization using TMM (trimmed Mean of M-values)
dge <- edgeR::calcNormFactors(object=dge)

# Model fitting & estimating dispersions
dge <- edgeR::estimateDisp(y = dge)

# Testing for diff exp
et <- edgeR::exactTest(object=dge)

# Extract table with FDR
top_degs <- edgeR::topTags(object=et,
                           n="Inf")

# Summary of DGE analysis results
# +ve logFC indicates gexp is higher in trt condition compared to ctr.
summary(limma::decideTests(object=et,
                           lfc=1))    # min logFC=1

## Results
# Write to CSV
write.csv(as.data.frame(top_degs), 
          file="~/Bioinformatics/R/results/trt_vs_ctr_dge.csv")

## Visualizations ---------
# Volcano plot
# format DGE table
diffexp_df <- top_degs$table
diffexp_df$Gene <- row.names(diffexp_df)
diffexp_df <- diffexp_df %>%
  dplyr::relocate(Gene, .before=logFC)
# Is each gene differentially expressed?
diffexp_df$significant <- "NO"
diffexp_df$significant[diffexp_df$logFC > 0.6 &
                         diffexp_df$FDR < 0.05] <- "UP"
diffexp_df$significant[diffexp_df$logFC < -0.6 &
                         diffexp_df$FDR < 0.05] <- "DOWN"
# Names of differentially expressed genes
diffexp_df$delabel <- NA
diffexp_df$delabel[diffexp_df$significant!="NO"] <-
  diffexp_df$Gene[diffexp_df$significant!="NO"]

# plot
ggplot2::ggplot(data=diffexp_df, 
                ggplot2::aes(x=logFC, 
                             y=-log10(FDR), 
                             col=significant, 
                             label=delabel)) +
  ggplot2::geom_point() + 
  ggplot2::theme_minimal() +
  ggrepel::geom_text_repel() +
  ggplot2::scale_color_manual(values=c("blue", "black", "red")) +
  ggplot2::geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  ggplot2::geom_hline(yintercept=-log10(0.05), col="red")


# Heatmap
# select top and bottom 20 significant genes with |logFC|>0.6
selected_genes <- diffexp_df %>%
  dplyr::arrange(desc(logFC)) %>%
  dplyr::filter(significant!="NO") %>%
  dplyr::pull(Gene)
selected_genes <- c(selected_genes[1:20],
                    selected_genes[(length(selected_genes)-20+1):length(selected_genes)])
# gexp matrix subset
sub_mat <- gexp_mat[row.names(gexp_mat) %in% selected_genes, ]  #subset
sub_mat <- sub_mat[match(selected_genes,row.names(sub_mat)), ]  #reorder
# normalize & scale gexp matrix
col_names <- colnames(sub_mat)
sub_mat <- t(apply(sub_mat,1,scale))
colnames(sub_mat) <- col_names

# plot
ComplexHeatmap::Heatmap(matrix=sub_mat,
                        col=circlize::colorRamp2(c(-2,0,2),
                                                 c("blue","white","red")),
                        name="Relative gene expression",
                        cluster_rows=FALSE,
                        cluster_columns=FALSE,
                        column_title="Samples",
                        row_title="Genes",
                        show_column_names=TRUE,
                        show_row_names=TRUE,
                        top_annotation=ComplexHeatmap::HeatmapAnnotation(df=dge$samples["group"]),
                        column_split = c("","",""," "," "," ")
                        )

