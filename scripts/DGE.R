#Differential gene expression using DESeq2 
#set working directory

#load libraries
library(tximport)
library(ggplot2)
library(rhdf5) # enables the reading of abundace.h5 files
library(DESeq2)
library(genefilter)
#quantification of genes done using kallisto

#define your directory by using the path that contains your data and infos
dir<-'/Users/fufordamaris/Desktop/My Data'

#read and parse your data
samples = read.table("/Users/fufordamaris/Desktop/My Data/samples.tsv", sep = "\t", header = TRUE)
head(samples, 5)

#define path to your seq data and bring them to a common directory
files <- file.path(dir, "RNAseq_Quant", samples$sample, "abundance.h5")
head(files, 3)

#assign names to your abundances
names(files) <- samples$sample
head(files, 5)

#read transcript level information and construct counts table.
txi.kallisto <- tximport(files, type = "kallisto", txOut = TRUE)
countdata <- round(txi.kallisto$counts)

# ! added by Chris !
# I suggest to remove genes with very low counts prior to the DESeq analysis
filter    <- which( rowSums(countdata >= 10) >= 2 )
countdata <- countdata[filter,]

#generate your metadata which contains information required for your deseqdata object
metadata <- read.table("/Users/fufordamaris/Desktop/My Data/samples.tsv", sep = "\t", header = TRUE)
head(metadata, 5)
rownames(metadata) <- metadata$sample

#create your table for DESeq analysis downstream
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=metadata, design= ~batch + strain)

# ! adjusted by Chris !
# you have to assign the result of the DESeq call to a variable (named fit here)
#run DESeq pipeline
fit <- DESeq(dds)
#warning message: 33 rows did not converge in beta, labelled in mcols(object)$betaConv. Use larger maxit argument with nbinomWaldTest
# Chris: this is just a warning, no error

# in our deseq object, dds 33 rows did not converge in beta, labelled in mcols(object)$betaConv.

#visualising significant results
results <- results(fit)
#filter results in terms of log2FC as well as adj p values 
summary(results)
resultsNames(fit)

#Comparing results
NT12002 <- results(fit, contrast=c("strain", "NT12002", "NT12001")) #comparing expression between two strains, the position of the last strain in our baseline for the comparison
write.csv(as.data.frame(NT12002), "/Users/fufordamaris/Desktop/My Data/DESeqOutput/NT12002.csv") # save our results to a csv file for further processing 

#plotting MA plots
plotMA(results)
     NT12004 <- results(fit, contrast=c("strain", "NT12002", "NT12001"))
plotMA(NT12002, main='NT12002')


#Generating a loop to write out the results (contrasts) between strains and our reference strain (EK12)

#read table with strain info, ensuring to eliminating the row containing strain id of ref strain and this will create error message when generating our contrasts downstream
strains <- samples$strain
strains <- unique(strains)
strains <- strains[-1]

for(istrain in strains){
  istraintable <- results(fit, contrast=c("strain", istrain, "NT12001"))
  write.csv(as.data.frame(istraintable),paste("/Users/fufordamaris/Desktop/My Data/DESeqOutput/",istrain,".csv"))
}

#Visualisation of results
#MA plots
for(istrain in strains){
  png(file=paste("/Users/fufordamaris/Desktop/My Data/DESeqOutput/Chart_plots/",istrain,".png", sep = ""))
  istraintable <- results(fit, contrast=c("strain", istrain, "NT12001"))
  plotMA(istraintable, main=paste(istrain,"vs NT2001"))
  dev.off()
}

#filtering results to get most significant genes only 

for(istrain in strains){
  istraintable <- results(fit, alpha = 0.05, lfcThreshold = 1, contrast=c("strain", istrain, "NT12001"))
  write.csv(as.data.frame(istraintable),paste("/Users/fufordamaris/Desktop/My Data/DESeqOutput/deseq2_subset/",istrain,".csv"))
}

#plotting MA plots using loop 

for(istrain in strains){
  png(file=paste("/Users/fufordamaris/Desktop/My Data/DESeqOutput/deseq2_subset/MA plots/",istrain,".png", sep = ""))
  istraintable <- results(fit, alpha = 0.05, lfcThreshold = 1, contrast=c("strain", istrain, "NT12001"))
  plotMA(istraintable, main=paste(istrain,"vs NT2001"))
  dev.off()
}
 




#other normalization   

# correcting for batch effect using limma 
# install and load package library 
#if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
#BiocManager::install("limma")


#normalising counts using vst 
vsd <- vst(dds, blind=FALSE)

#taking care of batch effects using limma
library(limma)

#prep data and correct for the batch effects
mat <- assay(vsd)
vst_batch_corrected <- limma::removeBatchEffect(mat, vsd$batch) # removing batch effects
write.csv(vst_batch_corrected,"vst_batch_corrected.csv")  # write output to a dataframe 
write.csv(mat, "batchwithbatch.csv") # vst with batch effect


#visualistation variable genes across samples 
# using normalized, vst-transformed , and batch corrected data
# calculate variance  per gene across samples

variablegenes <- rowVars(vst_batch_corrected)
write.csv(variablegenes, 'genechanges.csv')
# sort, so that the most variable genes will be on top of the object
topgenes <- head(order(variablegenes, decreasing = TRUE),30)
heatmap.2(assay(vsd)[topgenes,], scale = "row", trace="none", dendogram="column",margins=c(8,18), main='Top 30 Variable genes')

#1. sample distance 
#method Euclidean distance between samples
#estimating sample distances : assess overall similarity or differences between samples
#transpose our matrix since dist function expects different samples to be rows for its argument and different dimensions(genes) to be columns
#method Euclidean distance between samples

sampleDists <- dist(t(assay(vsd))) # t transposes our data matrix so we have samples occupying rows and columns
sampleDists
#visualise our results 
sampleDistsMatrix <- as.matrix(sampleDists) 
rownames(sampleDistsMatrix) <- paste(vsd$strain)
colnames(sampleDistsMatrix) <- NULL
library("gplots")
library("RColorBrewer")
colors <- colorRampPalette(rev(brewer.pal(7,"RdYlBu")))(100)
fontsize <- 3
jpeg(filename = "sampledistance_HM.jpeg")
pheatmap(sampleDistsMatrix,clustering_distance_rows = sampleDists, clustering_distance_cols = sampleDists,
         col = colors, cellheight = 3, fontsize_row = fontsize)



#normalisation using rld 
rld <- rlog(dds)
#correcting for batch effects 
rld_batch_corrected <- limma::removeBatchEffect(assay(rld), rld$batch)
write.csv(rld_batch_corrected,"rld_batch_corrected.csv") 





