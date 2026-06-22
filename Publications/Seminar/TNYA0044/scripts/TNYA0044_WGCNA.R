## TNYA0044 Weighted gene co-expression network analysis (WGCNA) on bulk RNAseq data
## Reva Shenwai
## 24-Jun-2024
## ---------------------------------------------------

# Set parameters
install_pkg <- FALSE    #do packages need to be installed?
wd <- "~/Bioinformatics/R/results"
minModuleSize=15

# Install packages if required
if (install_pkg==TRUE) {
  install.packages("WGCNA")
  install.packages("tidyverse")
  install.packages("magrittr")
  BiocManager::install("impute")
  BiocManager::install("preprocessCore")
  install.packages("corrplot")
  install.packages("pals")
}

# Load libraries
library(WGCNA)
library(tidyverse)
library(magrittr)
library(impute)
library(preprocessCore)
library(corrplot)
library(pals)
library(dynamicTreeCut)

# Load data
datExpr <- read.delim(paste0("~/Bioinformatics/Interviews/Seminar/TNYA0044/",
                             "data/GSE253225_TNYA0044.TPM.txt"),
                      row.names=1)
datExpr <- datExpr[-1,]
datExpr <- t(datExpr)   # format for WGCNA tool

datTraits <- read.csv(paste0("~/Bioinformatics/Interviews/Seminar/TNYA0044/",
                             "data/TNYA0044_metadata.csv"),
                      row.names=1)
#datTraits <- datTraits[,2:3]  # no numeric traits
corr_method <- "spearman"
network_plot_threshold <- 0.1

# settings
options(stringsAsFactors=FALSE)
WGCNA::enableWGCNAThreads()

wd <- "~/Bioinformatics/Interviews/Seminar/TNYA0044"
setwd(wd)

# Cluster samples to find outliers
sampleTree <- stats::hclust(dist(datExpr), method="average")
# plot
plot(sampleTree, main="Sample clustering to detect outliers",
     sub="",xlab="",cex.lab=1.5,cex.axis=1.5,cex.main=2)
# ask for user input to remove outgroups from the above plot
cat("Base on the 'Sample clustering to detect outliers' plot indicate which SampleID's\n",
    "to remove. Separate multiple ID's with commas; leave blank if not removing any ID's.")
outlier_samples <- readline(prompt="Indicate outlier sample ID's:")
# remove outlier samples
outlier_samples <- unlist(strsplit(outlier_samples,","))
outlier_samples <- gsub(" ","",outlier_samples)
if(length(outlier_samples)>0){
  datExpr <- datExpr[,-outlier_samples]
}

gene.names=colnames(datExpr)

powers <- c(c(1:100), seq(from=12,to=20,by=2))
sft <- WGCNA::pickSoftThreshold(datExpr,dataIsExpr=TRUE,powerVector=powers,corFnc=cor,
                                corOption=list(use='p'),networkType="unsigned")

# Plot results
sizeGrWindow(9,5)
par(mfrow=c(1,2))
cex1=0.9
# Scale-free topology fit index as fxn of soft-thresholding power
plot(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft threshold (power)", ylab="Scale free topology model fit, signed R^2",
     type="n",main="Scale independence")
text(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# Red line corresponds to using R^2 cutoff
abline(h=0.8,col="red")
# Mean connectivity as fxn of solf-thresholding power
plot(sft$fitIndices[,1],sft$fitIndices[,5],xlab="Soft thresholding (power)",
     ylab="Mean connectivity",type="n",main="Mean connectivity")
text(sft$fitIndices[,1],sft$fitIndices[,5],labels=powers,cex=cex1,col="red")

# Generate adjacency & TOM (topological overlap) similarity matrices
softPower=7
#adj=WGCNA::adjacency(datExpr,type="unsigned",power=softPower)
TOM=WGCNA::TOMsimilarityFromExpr(datExpr,networkType="unsigned",
                                 TOMType="unsigned",power=softPower)
# row & col names set to genes
colnames(TOM) <- rownames(TOM) <- gene.names
dissTOM=1-TOM

# Module detection
geneTree <- hclust(as.dist(dissTOM),method="average")
plot(geneTree,xlab="",cex=0.3)

# Module identification
dynamicMods <- dynamicTreeCut::cutreeDynamic(dendro=geneTree,method="tree",
                                             minClusterSize=minModuleSize)
table(dynamicMods)    # look at module sizes

# assign colors to modules
dynamicColors <- WGCNA::labels2colors(dynamicMods)
table(dynamicColors)
# Plot module assignment under dendrogram
WGCNA::plotDendroAndColors(geneTree, dynamicColors,"Dynamic tree cut",
                           dendroLabels=FALSE,hang=0.03,addGuide=TRUE,
                           guideHang=0.05,main="Gene dendrogram & module colors")

# Discard unsaaigned gnes
restGenes <- (dynamicColors!="grey")
diss1 <- 1 - WGCNA::TOMsimilarityFromExpr(datExpr[,restGenes],power=softPower)

colnames(diss1) <- rownames(diss1) <- gene.names[restGenes]
hier1=hclust(as.dist(diss1),method="average")
WGCNA::plotDendroAndColors(hier1,dynamicColors[restGenes],"Dynamic tree cut",
                           dendroLabels=FALSE,hang=0.03,addGuide=TRUE,guideHang=0.05,
                           main="Gene dendrogram & module colors")
diag(diss1) <- NA   # set diagonal to NA

# Visualize TOMplot
WGCNA::sizeGrWindow(7,7)
WGCNA::TOMplot(diss1,hier1,as.character(dynamicColors[restGenes]))

# Extract modules
module_colors <- setdiff(unique(dynamicColors),"grey")
for(color in module_colors){
  module <- gene.names[which(dynamicColors==color)]
  write.table(module,paste0("module_",color,".txt"),
              sep="\t",row.names=FALSE,col.names=FALSE,
              quote=FALSE)
}

# gexp patterns of clustered genes
module.order <- unlist(tapply(1:ncol(datExpr),as.factor(dynamicColors),I))
m <- t(t(datExpr[,module.order])/apply(datExpr[,module.order],2,max))
heatmap(t(m),zlim=c(0,1),col=gray.colors(100),Rowv=NA,Colv=NA,labRow=NA,
        scale="none",RowSideColors=dynamicColors[module.order])

# Quantify module similarity by eigengene correlation
MEList <- WGCNA::moduleEigengenes(datExpr,colors=dynamicColors)
MEs <- MEList$eigengenes
WGCNA::plotEigengeneNetworks(MEs,"",marDendro=c(0,4,1,2),
                             marHeatmap=c(3,4,1,2))
