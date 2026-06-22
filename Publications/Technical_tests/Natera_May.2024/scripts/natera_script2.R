# R script

# First tab, "Samples", contains the information about 16 samples used in this study, 
# first column has sample ID, second is the cell line and third is the drug used. 
# Second tab, "Counts", contains RNA-seq data, with the first two columns indicating 
# the genes, and the next 16 columns are gene counts from corresponding samples.
# Using R or Python, do the following:

# 1. Run unsupervised clustering for this dataset, ideally with two different approaches. Visualize the results.
# 2. For each cell line, run differential gene expression analysis between two conditions - Untreated versus Dex. 
# Visualize the results.

# You see this task is not described in much detail - just use your best judgement and common sense!

# libraries
library(openxlsx)
library(edgeR)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(ComplexHeatmap)
library(circlize)

# load data
samples_df <- openxlsx::read.xlsx("~/Bioinformatics/Interviews/Natera_test2_data.xlsx",
                                  sheet="Samples")
counts_df <- openxlsx::read.xlsx("~/Bioinformatics/Interviews/Natera_test2_data.xlsx",
                                   sheet="Counts")


# Unsupervised clustering
# approach 1:



# DGE analysis
unique(samples_df$Cells)
unique(samples_df$Drug)

## Input data
# sample gene expression matrix
gexp_mat <- counts_df
rownames(gexp_mat) <- counts_df$gene_id
gexp_mat <- gexp_mat[,3:ncol(gexp_mat)]

# meta data
meta_data <- samples_df
rownames(meta_data) <- samples_df$Sample
meta_data <- meta_data[,2:3]

# match metadata & gexp samples order
meta_data <- meta_data[match(colnames(gexp_mat),rownames(meta_data)),]
identical(colnames(gexp_mat),rownames(meta_data))

for(cell_type in unique(samples_df$Cells)) {
  # TODO: subset gexp & metadata DFs to 1 cell type
  # Reassign meta_data as meta_data$Drug
  
  # DGElist data object
  dge <- edgeR::DGEList(counts=gexp_mat,
                        group=meta_data)
  
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
}