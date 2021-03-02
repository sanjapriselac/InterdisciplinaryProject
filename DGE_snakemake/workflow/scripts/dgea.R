#!/opt/conda/envs/R_env/bin/Rscript --vanilla
# Sanja Priselac
# 02/03/2021

##################
library(ggplot2)
library(limma)
library(VennDiagram)
##################

for (i in 1:13) {
  load(snakemake@input[[i]])
}

condition <- factor(paste(samples$treatment, samples$resistance, sep="."))
design <- model.matrix(~ 0 + condition) 
colnames(design) <- gsub("condition", "", colnames(design))

## Apply the limma-voom method:
v <- voomWithQualityWeights(y_qn, design, plot=T)

v.05 <- voomWithQualityWeights(y_05, design, plot=T)
v.1 <- voomWithQualityWeights(y_1, design, plot=T)
v.25 <- voomWithQualityWeights(y_25, design, plot=T)
v.5 <- voomWithQualityWeights(y_5, design, plot=T)
v.10 <- voomWithQualityWeights(y_10, design, plot=T)


fit <- lmFit(v)

fit.05 <- lmFit(v.05)
fit.1 <- lmFit(v.1)
fit.25 <- lmFit(v.25)
fit.5 <- lmFit(v.5)
fit.10 <- lmFit(v.10)

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

fit2.05 <- contrasts.fit(fit.05, cont.matrix)
fit2.05 <- eBayes(fit2.05)
fit2.1 <- contrasts.fit(fit.1, cont.matrix)
fit2.1 <- eBayes(fit2.1)
fit2.25 <- contrasts.fit(fit.25, cont.matrix)
fit2.25 <- eBayes(fit2.25)
fit2.5 <- contrasts.fit(fit.5, cont.matrix)
fit2.5 <- eBayes(fit2.5)
fit2.10 <- contrasts.fit(fit.10, cont.matrix)
fit2.10 <- eBayes(fit2.10)

results <- decideTests(fit2, p.value=0.1) 

results.05 <- decideTests(fit2.05, p.value=0.05) 
results.1 <- decideTests(fit2.1, p.value=0.05) 
results.25 <- decideTests(fit2.25, p.value=0.05) 
results.5 <- decideTests(fit2.5, p.value=0.05) 
results.10 <- decideTests(fit2.10, p.value=0.05) 

genes.up <- list(qn = which(results@.Data[, 1]==1), dp0.5 = which(results.05@.Data[, 1]==1), dp1 = which(results.1@.Data[, 1]==1), 
                 dp2.5 = which(results.25@.Data[, 1]==1), dp5 = which(results.5@.Data[, 1]==1), dp10 = which(results.10@.Data[, 1]==1))

genes.down <- list(qn = which(results@.Data[, 1]==-1), dp0.5 = which(results.05@.Data[, 1]==-1), dp1 = which(results.1@.Data[, 1]==-1), 
                   dp2.5 = which(results.25@.Data[, 1]==-1), dp5 = which(results.5@.Data[, 1]==-1), dp10 = which(results.10@.Data[, 1]==-1))

## Venn diagrams
venn.diagram(genes.up[1:4], filename = "/tmp/repo/DGE_snakemake/results/plots/venn.up1.png", 
             fill = c("#999999", "#E69F00", "#56B4E9", "#009E73"), imagetype="png")
venn.diagram(genes.up[c(1, 4:6)], filename = "/tmp/repo/DGE_snakemake/results/plots/venn.up2.png", 
             fill = c("#999999", "#E69F00", "#56B4E9", "#009E73"), imagetype="png")
venn.diagram(genes.down[1:4], filename = "/tmp/repo/DGE_snakemake/results/plots/venn.down1.png", 
             fill = c("#999999", "#E69F00", "#56B4E9", "#009E73"), imagetype="png")
venn.diagram(genes.down[c(1, 4:6)], filename = "/tmp/repo/DGE_snakemake/results/plots/venn.down2.png", 
             fill = c("#999999", "#E69F00", "#56B4E9", "#009E73"), imagetype="png")

###################################################################################################################
## Volacno Plots
###################################################################################################################

# i corresponds to the coeff, meaning the contrasts
volcano_plot <- function (i) {
  table <- topTable(fit2, coef=i, number=Inf, sort.by="P")
  gene_names <- attributes(table)$row.names   
  xx <- table$logFC
  name <- names(fit2$contrasts[1,])[i]
  yy <- table$adj.P.Val
  yy <- -log10(yy)
  result <- results@.Data[gene_names,i]
  data <- data.frame(xx, yy, result)
  p <- ggplot(data, aes(x=xx, y=yy, col=result)) + 
    geom_point(size = 0.7) +
    theme_minimal() +
    xlab(paste("Log folds for", name)) +
    ylab("-log10(adj.p.value)") +
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

###################################
# Extraction of differentially expressed genes
###################################

fits <- list(fit2, fit2.05, fit2.1, fit2.25, fit2.5, fit2.10)
for (i in 1:length(fits)) {
  ## Treatment:
  treatmentInRGenes <- topTable(fits[[i]], coef=1, number=Inf, sort.by="P")
  treatmentInSGenes <- topTable(fits[[i]], coef=2, number=Inf, sort.by="P")
  treatmentGenes <- topTable(fits[[i]], coef=3, number=Inf, sort.by="P") 
  
  ## Resistance:
  resistanceInCGenes <- topTable(fits[[i]], coef=4, number=Inf, sort.by="P")
  resistanceInUGenes <- topTable(fits[[i]], coef=5, number=Inf, sort.by="P")
  resistanceGenes <- topTable(fits[[i]], coef=6, number=Inf, sort.by="P")
  
  write.csv(treatmentInRGenes, file=snakemake@output[[i]])
  write.csv(treatmentInSGenes, file=snakemake@output[[6+i]])
  write.csv(treatmentGenes, file=snakemake@output[[12+i]])
  write.csv(resistanceInCGenes, file=snakemake@output[[18+i]])
  write.csv(resistanceInUGenes, file=snakemake@output[[24+i]])
  write.csv(resistanceGenes, file=snakemake@output[[30+i]])
}

