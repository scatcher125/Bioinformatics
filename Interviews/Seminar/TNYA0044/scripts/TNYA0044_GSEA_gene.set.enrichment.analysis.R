## Gene set enrichment analysis (GSEA) on bulk RNAseq data for TNYA0044 dataset
## Reva Shenwai
## 21-Jun-2024
## ---------------------------------------------------

  # Set parameters
install_pkg <- FALSE    #do packages need to be installed?

# Install packages if required
if (install_pkg==TRUE) {
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("clusterProfiler")
  install.packages("DOSE")
  install.packages("msigdbr")
  install.packages("dplyr")
  install.packages("enrichplot")
  install.packages("ggupset")
  install.packages("ggridges")
  BiocManager::install("org.Mm.eg.db")
}

# Load libraries
library(clusterProfiler)
library(DOSE)
library(msigdbr)
library(dplyr)
library(enrichplot)
library(ggupset)
library(ggridges)
library(org.Mm.eg.db)

## Load data
comparisons <- c("TNYA0044.WT_HBSS.vs.PKP2_HBSS",
                 "TNYA0044.PKP2_HBSS.vs.PKP2_m3e13",
                 "TNYA0044.PKP2_HBSS.vs.PKP2_m1e14",
                 "TNYA0044.PKP2_HBSS.vs.PKP2_m1e14_2.5wk")
gene.rank_list <- lapply(comparisons,
                         function(file_nm_base) {
                           # read from file
                           df <- read.csv(paste0("~/Bioinformatics/Interviews/Seminar/TNYA0044/results/",
                                                 file_nm_base,".csv"))
                           # format
                           df <- df %>% dplyr::rename(Gene="X")
                           # make rank value
                           df$rank_stat <- sign(df$logFC) * -log10(df$FDR)
                           # convert gene names to ENTREZID
                           Gene_ENTREZID <- clusterProfiler::bitr(df$Gene,
                                                                  fromType="SYMBOL",
                                                                  toType="ENTREZID",
                                                                  OrgDb=org.Mm.eg.db)
                           # subset to mapped genes
                           temp_df <- subset(na.omit(df[,c("Gene","rank_stat")]),
                                       Gene %in% Gene_ENTREZID$SYMBOL)
                           # remove duplicated gene symbols
                           Gene_ENTREZID <- distinct(Gene_ENTREZID, SYMBOL, .keep_all=TRUE)
                           # add rank_stat to ENTREZID
                           Gene_ENTREZID$rank_stat <- temp_df$rank_stat
                           rm(temp_df)  # remove old data
                           # sort by decreasing order of rank_stat
                           Gene_ENTREZID <- Gene_ENTREZID %>%
                             dplyr::arrange(desc(rank_stat))
                           
                           # gene rank list
                           rank_vec <- Gene_ENTREZID$rank_stat
                           names(rank_vec) <- Gene_ENTREZID$ENTREZID
                           # return ranked gene list
                           return(rank_vec)
                         })
names(gene.rank_list) <- comparisons
                           
#gene <- names(geneList)[abs(geneList) > 2]

## GSEA analysis
# get canonical pathway genesets from MSigDB
c2.cp_genesets <- msigdbr::msigdbr(species="Mus musculus",
                             category="C2",
                             subcategory="CP") %>%
  dplyr::select(gs_name,entrez_gene)

# run analysis for each comparison
gsea.results_list <- lapply(gene.rank_list,
                            function(rank_vec) {
                              gsea_results <- clusterProfiler::GSEA(geneList=rank_vec,
                                                                    TERM2GENE=c2.cp_genesets)
                            })
names(gsea.results_list) <- comparisons

## Visualizations
# loop
lapply(gsea.results_list, function(gsea_results) {
  # GSEA rug plot
  enrichplot::gseaplot2(gsea_results,
                        geneSetID=5,
                        title=gsea_results@result$Description[5]
                        )
})

# loop
lapply(gsea.results_list, function(gsea_results) {
  # UpSet plot
  enrichplot::upsetplot(gsea_results)
})

# loop
lapply(gsea.results_list, function(gsea_results) {
  # Ridgeplot
  enrichplot::ridgeplot(gsea_results)
})

# NES plot showing top 20 Up & down regulated genesets with word clouds


