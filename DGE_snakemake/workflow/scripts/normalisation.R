#!/opt/conda/envs/R_env/bin/Rscript --vanilla
# Sanja Priselac
# 02/03/2021

########################
library(edgeR)
library(limma)
library(ggplot2)
library(rmutil)
########################

load(snakemake@input[[1]])
load(snakemake@input[[2]])

############################################################################################
# Normalisation
############################################################################################

y <- DGEList(txi$counts)

############################################################################################
# Quantile normalisation function 
############################################################################################

quartileNormalisation <- function(y) {
  ## y is DGEList object
  # save the indices of the sort
  orders <- apply(y$counts, 2, order, decreasing = TRUE)
  x.sort <- apply(y$counts, 2, sort, decreasing = TRUE)
  r.means <- rowMeans(x.sort)
  for (i in 1:ncol(x.sort)) {
    y$counts[orders[, i], i] <- r.means
  }
  return(y)
} 

############################################################################################
# Laplace Quantile normalisation function 
############################################################################################

quartileNormalisation.dp <- function(y, b, seed=100) {
  ## y is DGEList object
  # save the indices of the sort
  orders <- apply(y$counts, 2, order, decreasing = TRUE)
  x.sort <- apply(y$counts, 2, sort, decreasing = TRUE)
  r.means <- rowMeans(x.sort)
  
  ## laplace distributed random variables
  set.seed(seed)
  lap <- rlaplace(nrow(y$counts), m=0, s=b)
  ## add the noise and sort the vector, L3 
  r.means <- sort(r.means + lap, decreasing = TRUE)
  negatives <- which(r.means < 0)
  if (length(negatives) > 0) {
    warning(cat(length(negatives), 'r.means values negative and set to 0'))
    r.means[negatives] <- 0
  }
  for (i in 1:ncol(x.sort)) {
    y$counts[orders[, i], i] <- r.means
  }
  return(y)
} 


############################################################################################


y.qn <- quartileNormalisation(y)
unfilteredExpr <- cpm(y.qn, log=T)
#plotDensities(unfilteredExpr, col=myPalette[1:16], legend=FALSE)

cutoff <- 0 
#summary(unfilteredExpr > cutoff)
numSamplesWithExpression <- apply(unfilteredExpr, 1, function(x){ return(sum(x > cutoff)) })

## Most genes are either expressed in all samples, or in no sample; retain all genes expressed in at least half of the samples:
selectedGenes <- names(which(numSamplesWithExpression >= 8))
#length(selectedGenes) #9021 

# rebuild a new DGE object using only selected genes, and renormalize it
y <- DGEList(txi$counts[selectedGenes, ]) 

y_qn <- quartileNormalisation(y)
## calculate laplacian 
y_05 <- quartileNormalisation.dp(y, 0.5)
y_1 <- quartileNormalisation.dp(y, 1)
y_25 <- quartileNormalisation.dp(y, 2.5)
y_5 <- quartileNormalisation.dp(y, 5)
y_10 <- quartileNormalisation.dp(y, 10)

filteredExpr <- cpm(y_qn, log=T)
filteredExpr_05 <- cpm(y_05, log=T)
filteredExpr_1 <- cpm(y_1, log=T)
filteredExpr_25 <- cpm(y_25, log=T)
filteredExpr_5 <- cpm(y_5, log=T)
filteredExpr_10 <- cpm(y_10, log=T)

save(y_qn, file=snakemake@output[[1]])
save(y_05, file=snakemake@output[[2]])
save(y_1, file=snakemake@output[[3]])
save(y_25, file=snakemake@output[[4]])
save(y_5, file=snakemake@output[[5]])
save(y_10, file=snakemake@output[[6]])

save(filteredExpr, file=snakemake@output[[7]])
save(filteredExpr_05, file=snakemake@output[[8]])
save(filteredExpr_1, file=snakemake@output[[9]])
save(filteredExpr_25, file=snakemake@output[[10]])
save(filteredExpr_5, file=snakemake@output[[11]])
save(filteredExpr_10, file=snakemake@output[[12]])

############################################################################################
# Data clustering, only on quartile normalized data
############################################################################################

pca <- prcomp(t(filteredExpr), scale = T)
#plot(pca)
summary(pca)

loadings <- pca$rotation
scores <- pca$x

data <- data.frame(scores, as.factor(samples$resistance), as.factor(samples$treatment))
plotPCs <- function(i, j) {
  x_lab <- paste("PC", i,": ", round(summary(pca)$importance[2,i],3)*100, "% variance", sep="")
  y_lab <- paste0("PC", j,": ", round(summary(pca)$importance[2,j],3)*100, "% variance")
  
  p <- ggplot(data, aes(x=data[,i], y=data[,j], color=data[,17], shape=data[,18])) +
    geom_point() +
    xlab(x_lab) +
    ylab(y_lab) +
    labs(color="Resistance", shape="Treatment")
  return(p)
}

pc12 <- plotPCs(1, 2) #PC1 correlates with treatement 
ggsave("PCA12.pdf", plot=pc12, path = "/tmp/repo/DGE_snakemake/results/plots")

pc23 <- plotPCs(2, 3)
ggsave("PCA23.pdf", plot=pc23, path = "/tmp/repo/DGE_snakemake/results/plots")

pc34 <- plotPCs(3, 4) # PC4 correlates with the resistance 
ggsave("PCA34.pdf", plot=pc34, path = "/tmp/repo/DGE_snakemake/results/plots")


