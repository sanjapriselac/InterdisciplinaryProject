#!/opt/conda/envs/R_env/bin/Rscript --vanilla
# Sanja Priselac
# 17/03/2021

##############################


###########################################################################################################
ttab.qn <- read.csv(snakemake@input[[1]])
ttab.dp0.05 <- read.csv(snakemake@input[[2]])
ttab.dp0.1 <- read.csv(snakemake@input[[3]])
ttab.dp1 <- read.csv(snakemake@input[[4]])
ttab.dp5 <- read.csv(snakemake@input[[5]])
ttab.dp10 <- read.csv(snakemake@input[[6]])

rtab.qn <- read.csv(snakemake@input[[7]])
rtab.dp0.05 <- read.csv(snakemake@input[[8]])
rtab.dp0.1 <- read.csv(snakemake@input[[9]])
rtab.dp1 <- read.csv(snakemake@input[[10]])
rtab.dp5 <- read.csv(snakemake@input[[11]])
rtab.dp10 <- read.csv(snakemake@input[[12]])
###########################################################################################################

## tab main is the table where the first n p values are taken
test.orders <- function(tab.main, tab.sub, n=30) {
  tab.main <- (tab.main[order(tab.main$adj.P.Val), ])[1:n, c('X', 'adj.P.Val')] 
  tab.test <- merge(tab.main, tab.sub[,  c('X', 'adj.P.Val')], by='X')
  t <- friedman.test(as.matrix(tab.test[, -1]))
  return(t$p.value)
}

results <- data.frame(matrix(ncol=5, nrow = 2))
colnames(results) <- c('eps0.05', 'eps0.1', 'eps1', 'eps5', 'eps10')
rownames(results) <- c('p-value treat', 'p-value res')

treatment <- list(ttab.qn, ttab.dp0.05, ttab.dp0.1, ttab.dp1, ttab.dp5, ttab.dp10)
resistance <-  list(rtab.qn, rtab.dp0.05, rtab.dp0.1, rtab.dp1, rtab.dp5, rtab.dp10)

for (i in 2:length(treatment)) {
  results[1, i-1] <- test.orders(treatment[[1]], treatment[[i]])
  results[2, i-1] <- test.orders(resistance[[1]], resistance[[i]])
}

write.csv(results, file=snakemake@output[[1]])

