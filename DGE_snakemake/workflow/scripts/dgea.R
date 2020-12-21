#!/opt/conda/envs/R_env/bin/Rscript --vanilla
# Sanja Priselac

##################
library(ggplot2)
library(limma)
##################

load(snakemake@input[[1]])
load(snakemake@input[[2]])
load(snakemake@input[[3]])

condition <- factor(paste(samples$treatment, samples$resistance, sep="."))
design <- model.matrix(~ 0 + condition) 
colnames(design) <- gsub("condition", "", colnames(design))
## Apply the limma-voom method:
v <- voomWithQualityWeights(y, design, plot=T)

fit <- lmFit(v)

cont.matrix <- makeContrasts(
  treatmentInR  = Challenged.Resistant - Unchallenged.Resistant,
  treatmentInS  = Challenged.Susceptible - Unchallenged.Susceptible,
  treatment     = (Challenged.Resistant + Challenged.Susceptible) 
  - (Unchallenged.Resistant + Unchallenged.Susceptible),
  resistanceInC = Challenged.Resistant - Challenged.Susceptible,
  resistanceInU = Unchallenged.Resistant - Unchallenged.Susceptible,
  resistance    = (Challenged.Resistant + Unchallenged.Resistant) 
  - (Challenged.Susceptible + Unchallenged.Susceptible),
  levels=design)

#cont.matrix

####################################

fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
results <- decideTests(fit2, p.value=0.1) 
summary(results)

#volcanoplot(fit2, coef = 1, highlight = 0)
# fit2 for vulcano plot
# coeff: treatment in S, treatment in ...

volcano_plot <- function (i) {
  #xx <- as.matrix(fit2$coefficients)[,1] #same
  xx <- as.matrix(fit2@.Data[[1]][, i])
  name <- names(fit2$contrasts[1,])[i]
  yy <- as.matrix(fit2$p.value)[,i]
  yy <- -log10(yy)
  result <- results@.Data[,i]
  data <- data.frame(xx, yy, result)
  p <- ggplot(data, aes(x=xx, y=yy, col=result)) + 
    geom_point(size = 0.7) +
    theme_minimal() +
    xlab(paste("Coefficients for", name)) +
    ylab("-log10(p.value)") +
    ggtitle(paste0("Volcano plot for ", name))
  return(p)  
}

p1 <- volcano_plot(1)
ggsave(paste0("volcanoPlot_",names(fit2$contrasts[1,])[1], '.pdf'), plot=p1, path="/tmp/repo/DGE_snakemake/results/plots")
p2 <- volcano_plot(2)
ggsave(paste0("volcanoPlot_",names(fit2$contrasts[1,])[2], '.pdf'), plot=p2, path="/tmp/repo/DGE_snakemake/results/plots")
p3 <- volcano_plot(3)
ggsave(paste0("volcanoPlot_",names(fit2$contrasts[1,])[3], '.pdf'), plot=p3, path="/tmp/repo/DGE_snakemake/results/plots")
p4 <- volcano_plot(4)
ggsave(paste0("volcanoPlot_",names(fit2$contrasts[1,])[4], '.pdf'), plot=p4, path="/tmp/repo/DGE_snakemake/results/plots")
p5 <- volcano_plot(5)
ggsave(paste0("volcanoPlot_",names(fit2$contrasts[1,])[5], '.pdf'), plot=p5, path="/tmp/repo/DGE_snakemake/results/plots")
p6 <- volcano_plot(6)
ggsave(paste0("volcanoPlot_",names(fit2$contrasts[1,])[6], '.pdf'), plot=p6, path="/tmp/repo/DGE_snakemake/results/plots")

#vennDiagram(results[,c(1,2)])
#vennDiagram(results[,c(2,5)])
#vennDiagram(results[,c(2,6)])

###################################
# Extraction of differentially expressed genes
###################################

## Treatment:
treatmentGenes <- topTable(fit2, coef=3, p.value=0.1, number=Inf, sort.by="P")

## visualize the top 100 genes
#selectedExpression <- filteredExpr[rownames(treatmentGenes)[1:100],]
#heatmap.2(selectedExpression, scale="none", col = colors, margins = c(14, 6), trace='none', denscol="white", 
#          ColSideColors=myPalette[3:4][as.integer(as.factor(samples$treatment))])

## Resistance:
resistanceGenes <- topTable(fit2, coef=6, p.value=0.1, number=Inf, sort.by="P")

## visualize all DE genes
#selectedExpression <- filteredExpr[rownames(resistanceGenes),]
#heatmap.2(selectedExpression, scale="none", col = colors, margins = c(14, 6), trace='none', denscol="white", 
#          ColSideColors=myPalette[1:2][as.integer(as.factor(samples$resistance))])

save(treatmentGenes, file=snakemake@output[[1]])
save(resistanceGenes, file=snakemake@output[[2]])
