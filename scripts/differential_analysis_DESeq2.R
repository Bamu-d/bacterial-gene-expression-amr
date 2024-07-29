#set working directory and load libraries

library(tximport)
library(ggplot2)
library(rhdf5) # enables the reading of abundace.h5 files
library(DESeq2)
library(genefilter)

#read and parse  data
samples = read.table("samples.tsv", sep = "\t", header = TRUE)
head(samples, 5)

#define path to  seq data and bring them to a common directory
files <- file.path("RNAseq_Quant", samples$sample, "abundance.h5")
head(files, 3)

#assign names to  abundances
names(files) <- samples$sample
head(files, 5)

#read transcript level information and construct counts table.
txi.kallisto <- tximport(files, type = "kallisto", txOut = TRUE)
countdata <- round(txi.kallisto$counts)

#  remove genes with very low counts prior to the DESeq analysis
filter    <- which( rowSums(countdata >= 10) >= 2 )
countdata <- countdata[filter,]

#generate  metadata which contains information required for  deseqdata object
metadata <- read.table("samples.tsv", sep = "\t", header = TRUE)
head(metadata, 5)
rownames(metadata) <- metadata$sample

#create  table for DESeq analysis downstream
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=metadata, design= ~batch + strain)

#run DESeq pipeline
fit <- DESeq(dds)

#visualise significant results
results <- results(fit)

#filter results in terms of log2FC as well as adj p values 
summary(results)
resultsNames(fit)

# write out the results (contrasts) between strains and our reference strain (EK12)

strains <- samples$strain
strains <- unique(strains)
strains <- strains[-1]

for(istrain in strains){
  istraintable <- results(fit, contrast=c("strain", istrain, "NT12001"))
  write.csv(as.data.frame(istraintable),paste("DESeqOutput/",istrain,".csv"))
}