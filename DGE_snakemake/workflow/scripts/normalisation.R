#!/opt/conda/envs/R_env/bin/Rscript --vanilla
# Sanja Priselac

########################
library(edgeR)
library(limma)
library(ggplot2)
########################

load(snakemake@input[[1]])
load(snakemake@input[[2]])

############################################################################################
# Normalisation
############################################################################################

y <- DGEList(txi$counts)
y <- calcNormFactors(y, method="upperquartile")
y$samples

#myPalette <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"))
## Plot CPM distribution (counts per million)
#plotDensities(cpm(y), col=myPalette[1:16], legend=FALSE) #really skewed
unfilteredExpr <- cpm(y, log=T)
#plotDensities(unfilteredExpr, col=myPalette[1:16], legend=FALSE)

cutoff <- 0 
summary(unfilteredExpr > cutoff)
numSamplesWithExpression <- apply(unfilteredExpr, 1, function(x){ return(sum(x > cutoff)) })
#hist(numSamplesWithExpression)

## Most genes are either expressed in all samples, or in no sample; retain all genes expressed in at least half of the samples:
selectedGenes <- names(which(numSamplesWithExpression >= 8))
length(selectedGenes) #8990

# rebuild a new DGE object using only selected genes, and renormalize it
y <- DGEList(txi$counts[selectedGenes, ]) 
y <- calcNormFactors(y,  method="upperquartile")
filteredExpr <- cpm(y, log=T)
#plotDensities(filteredExpr, col=myPalette[1:16], legend=FALSE) #"topright"

#plotDensities(filteredExpr, group=samples$resistance, col=myPalette[1:2], legend='topright')
#plotDensities(filteredExpr, group=samples$treatment, col=myPalette[3:4], legend='topright')

############################################################################################
# Data clustering
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


save(filteredExpr, file=snakemake@output[[1]])
save(y, file=snakemake@output[[2]])
