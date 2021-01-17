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

## Treatment:
treatmentInRGenes <- topTable(fit2, coef=1, number=Inf, sort.by="P")
treatmentInSGenes <- topTable(fit2, coef=2, number=Inf, sort.by="P")
treatmentGenes <- topTable(fit2, coef=3, number=Inf, sort.by="P") 

## Resistance:
resistanceInCGenes <- topTable(fit2, coef=4, number=Inf, sort.by="P")
resistanceInUGenes <- topTable(fit2, coef=5, number=Inf, sort.by="P")
resistanceGenes <- topTable(fit2, coef=6, number=Inf, sort.by="P")

write.csv(treatmentInRGenes, file=snakemake@output[[1]])
write.csv(treatmentInSGenes, file=snakemake@output[[2]])
write.csv(treatmentGenes, file=snakemake@output[[3]])
write.csv(resistanceInCGenes, file=snakemake@output[[4]])
write.csv(resistanceInUGenes, file=snakemake@output[[5]])
write.csv(resistanceGenes, file=snakemake@output[[6]])