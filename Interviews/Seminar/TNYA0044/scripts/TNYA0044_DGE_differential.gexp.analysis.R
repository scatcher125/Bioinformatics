## TNYA0044 Differential gene expression analysis on bulk RNAseq data
## - Script is based on https://web.stanford.edu/class/bios221/labs/rnaseq/lab_4_rnaseq.html
## Reva Shenwai
## 15-Jun-2024
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
gexp_mat <- read.delim("~/Bioinformatics/Interviews/Seminar/TNYA0044/data/GSE253225_TNYA0044.counts.txt",
                       row.names=1)
gexp_mat <- gexp_mat[-1,]

# meta data
meta_data <- read.csv(file="~/Bioinformatics/Interviews/Seminar/TNYA0044/data/TNYA0044_metadata.csv")
meta_data <- meta_data[,-1]
row.names(meta_data) <- meta_data$Sample_ID

identical(colnames(gexp_mat),rownames(meta_data))   # check

# DGElist data object
dge <- edgeR::DGEList(counts=gexp_mat,
                      group=factor(meta_data$Group,
                                   levels=c("WT_HBSS","PKP2_HBSS","PKP2_m3e13","PKP2_m1e14","PKP2_m1e14_2.5wk")))

# group colors
group_colors <- c(WT_HBSS="black",
                  PKP2_HBSS="blue",
                  PKP2_m3e13="red",
                  PKP2_m1e14="red4",
                  PKP2_m1e14_2.5wk="forestgreen"
                  )

## Analysis: single factor design
# Filter genes with low counts
keep <- edgeR::filterByExpr(y=dge)
dge <- dge[keep, , keep.lib.sizes=FALSE]

# Normalization using TMM (trimmed Mean of M-values)
dge <- edgeR::calcNormFactors(object=dge)

# Multi-dimensional scaling plot
# - Do samples from the same group cluster together?
limma::plotMDS(dge, method="bcv", 
               col=group_colors[dge$samples$group])
legend("topright", as.character(unique(dge$samples$group)),
       col=group_colors, pch=20, pt.cex=2)

# Model fitting & estimating dispersions
dge <- edgeR::estimateDisp(dge)
names(dge)

# plot tagwise biological coefficient of variation
plotBCV(dge)

# Testing for diff exp
et.WT_HBSS.vs.PKP2_HBSS <- edgeR::exactTest(object=dge, pair=c("WT_HBSS","PKP2_HBSS"))
et.PKP2_HBSS.vs.PKP2_m3e13 <- edgeR::exactTest(object=dge, pair=c("PKP2_HBSS","PKP2_m3e13"))
et.PKP2_HBSS.vs.PKP2_m1e14 <- edgeR::exactTest(object=dge, pair=c("PKP2_HBSS","PKP2_m1e14"))
et.PKP2_HBSS.vs.PKP2_m1e14_2.5wk <- edgeR::exactTest(object=dge, pair=c("PKP2_HBSS","PKP2_m1e14_2.5wk"))

edgeR::topTags(et.WT_HBSS.vs.PKP2_HBSS, n=10)

# Extract table with FDR
top.degs_WT_HBSS.vs.PKP2_HBSS <- edgeR::topTags(object=et.WT_HBSS.vs.PKP2_HBSS, n="Inf")
top.degs_PKP2_HBSS.vs.PKP2_m3e13 <- edgeR::topTags(object=et.PKP2_HBSS.vs.PKP2_m3e13, n="Inf")
top.degs_PKP2_HBSS.vs.PKP2_m1e14 <- edgeR::topTags(object=et.PKP2_HBSS.vs.PKP2_m1e14, n="Inf")
top.degs_PKP2_HBSS.vs.PKP2_m1e14_2.5wk <- edgeR::topTags(object=et.PKP2_HBSS.vs.PKP2_m1e14_2.5wk, n="Inf")

# Summary of DGE analysis results
# +ve logFC indicates gexp is higher in trt condition compared to ctr.
summary(limma::decideTests(object=et.WT_HBSS.vs.PKP2_HBSS))
summary(limma::decideTests(object=et.PKP2_HBSS.vs.PKP2_m3e13))
summary(limma::decideTests(object=et.PKP2_HBSS.vs.PKP2_m1e14))
summary(limma::decideTests(object=et.PKP2_HBSS.vs.PKP2_m1e14_2.5wk))

## Results
# Write to CSV
write.csv(as.data.frame(top.degs_WT_HBSS.vs.PKP2_HBSS), 
          file="~/Bioinformatics/Interviews/Seminar/TNYA0044/results/TNYA0044.WT_HBSS.vs.PKP2_HBSS.csv")
write.csv(as.data.frame(top.degs_PKP2_HBSS.vs.PKP2_m3e13), 
          file="~/Bioinformatics/Interviews/Seminar/TNYA0044/results/TNYA0044.PKP2_HBSS.vs.PKP2_m3e13.csv")
write.csv(as.data.frame(top.degs_PKP2_HBSS.vs.PKP2_m1e14), 
          file="~/Bioinformatics/Interviews/Seminar/TNYA0044/results/TNYA0044.PKP2_HBSS.vs.PKP2_m1e14.csv")
write.csv(as.data.frame(top.degs_PKP2_HBSS.vs.PKP2_m1e14_2.5wk), 
          file="~/Bioinformatics/Interviews/Seminar/TNYA0044/results/TNYA0044.PKP2_HBSS.vs.PKP2_m1e14_2.5wk.csv")

## Visualizations ---------
## Volcano plot function
make_volcano_plot <- function(edgeR_topTags) {
  # format DGE table
  diffexp_df <- edgeR_topTags$table
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
    ggplot2::geom_vline(xintercept=c(-0.6, 0.6), col="red4") +
    ggplot2::geom_hline(yintercept=-log10(0.05), col="red4")
}

# WT_HBSS.vs.PKP2_HBSS
make_volcano_plot(top.degs_WT_HBSS.vs.PKP2_HBSS)

# PKP2_HBSS.vs.PKP2_m3e13
make_volcano_plot(top.degs_PKP2_HBSS.vs.PKP2_m3e13)

# PKP2_HBSS.vs.PKP2_m1e14
make_volcano_plot(top.degs_PKP2_HBSS.vs.PKP2_m1e14)

# PKP2_HBSS.vs.PKP2_m1e14_2.5wk
make_volcano_plot(top.degs_PKP2_HBSS.vs.PKP2_m1e14_2.5wk)


## Heatmap
# Select top & bottom 20 genes from WT vs KO
#format DGE table
diffexp_df <- top.degs_WT_HBSS.vs.PKP2_HBSS$table
diffexp_df$Gene <- row.names(diffexp_df)
diffexp_df <- diffexp_df %>%
  dplyr::relocate(Gene, .before=logFC)
#Is each gene differentially expressed?
diffexp_df$significant <- "NO"
diffexp_df$significant[diffexp_df$logFC > 0.6 &
                         diffexp_df$FDR < 0.05] <- "UP"
diffexp_df$significant[diffexp_df$logFC < -0.6 &
                         diffexp_df$FDR < 0.05] <- "DOWN"
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
                        column_title="Top 20 DEGs for PKP2_HBSS vs WT_HBSS",
                        row_title="Genes",
                        show_column_names=TRUE,
                        show_row_names=TRUE,
                        top_annotation=ComplexHeatmap::HeatmapAnnotation(
                          df=dge$samples["group"],
                          col=list(group=group_colors)),
                        column_split = factor(meta_data$Group),
                        row_split=c(rep("up",20),rep("down",20)),
                        row_names_gp = gpar(fontsize = 10)
                        )

