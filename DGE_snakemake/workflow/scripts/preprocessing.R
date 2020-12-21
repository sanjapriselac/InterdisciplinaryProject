#!/opt/conda/envs/R_env/bin/Rscript --vanilla
# Sanja Priselac, 18.12.2020

#######################
library(GEOquery)
library(tximport)
library(biomaRt)
library(rhdf5)
#library(optparse)
#######################

#my_options = list(
#  make_option(c("-i", "--inputfile"), default='/tmp/repo/DGE_snakemake/results/kallisto'),
#  make_option(c("-o", "--outputfile"), default='/tmp/repo/DGE_snakemake/results/R')
#)

#args <- parse_args(OptionParser(option_list=my_options))

dataFolder <- '/tmp/repo/DGE_snakemake/data/additional/'
gse <- getGEO(filename=file.path(dataFolder, "GSE59411_series_matrix.txt"))
samples <- pData(phenoData(gse))[,c(1,2, 8:12, 45)]

samples$library_accession <- unlist(
  lapply(strsplit(as.character(samples$supplementary_file_1), split="/"), tail, n=1))

ENA_eunID <- read.csv("/tmp/repo/DGE_snakemake/data/additional/samples_info.csv")

samples <- merge(samples, ENA_eunID, by.x="title", by.y="Sample_Title")
samples <- samples[,c(2, 4:8, 3, 1, 10)]

## Rename the columns
names(samples)[3] <- "treatment"
names(samples)[4] <- "dgrp_line"
names(samples)[5] <- "resistance"

## Simplify the levels of experimental factors
samples$treatment <- gsub("treatment: ", "", samples$treatment)
samples$dgrp_line <- gsub("dgrp line: ", "", samples$dgrp_line)
samples$resistance <- gsub("resistance: ", "", samples$resistance)

## Path to Kallisto result files
## promijeni
files <- file.path(snakemake@input[[1]], samples$run_accession, "abundance.tsv")
#  paste('k', samples$run_accession, '.fastq.gz', sep='')
names(files) <- samples$title

#import the transciptID and geneID, Choose D. melanogaster dataset in Ensembl release 84
mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "dmelanogaster_gene_ensembl", host="mar2016.archive.ensembl.org")
## retrieve the mapping between transcript and gene IDs
transcript2gene <- getBM(attributes=c("ensembl_transcript_id", "ensembl_gene_id"), mart=mart) 
transcript2gene <- transcript2gene[order(transcript2gene$ensembl_gene_id), ]

# import Kallisto results
txi <- tximport(files, type = "kallisto", tx2gene = transcript2gene, countsFromAbundance="scaledTPM")

save(samples, file=snakemake@output[[1]])
save(txi, file=snakemake@output[[2]])

