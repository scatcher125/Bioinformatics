## Gene set enrichment analysis (GSEA) on bulk RNAseq data
## Reva Shenwai
## 19-May-2024
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
}

# Load libraries
library(clusterProfiler)
library(DOSE)
library(msigdbr)
library(dplyr)
library(enrichplot)
library(ggupset)
library(ggridges)

## Load data
data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]

## GSEA analysis
# get canonical pathway genesets from MSigDB
c2.cp_genesets <- msigdbr::msigdbr(species="Homo sapiens",
                             category="C2",
                             subcategory="CP") %>%
  dplyr::select(gs_name,entrez_gene)
# run analysis
gsea_results <- clusterProfiler::GSEA(geneList,
                            TERM2GENE=c2.cp_genesets)

## Visualizations
# GSEA rug plot
enrichplot::gseaplot2(gsea_results,
                      geneSetID=5,
                      title=gsea_results@result$Description[5]
)

# UpSet plot
enrichplot::upsetplot(gsea_results)

# Ridgeplot
enrichplot::ridgeplot(gsea_results)

