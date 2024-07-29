source("http://bibs.snu.ac.kr/software/LPEseq/LPEseq.R")
install.packages("local_folder/LPEseq_version.tar.gz", repos=NULL,
                 type="source")

#set working directory and load libraries. 

library(tximport)

#read and parse your data

samples = read.table("sample_info.txt", sep = "\t", header = TRUE)
head(samples, 5)

#define path to your seq data and bring them to a common directory
files <- file.path("pseudomonas_quant", samples$Identifier, "abundance.h5")
head(files, 3)

#assign names to your abundances
names(files) <- samples$Identifier
head(files, 5)

#read transcript level information and construct counts table.
txi.kallisto <- tximport(files, type = "kallisto", txOut = TRUE)
countdata <- round(txi.kallisto$counts)

#NORMALISATION 
#Helps remove the effect of sequencing depth
#we are making use of raw counts, if normalised counts were used, a simple log2 transformation will be sufficient. 
counts.norm <- LPEseq.normalise(countdata)

#TESTING DIFFERENTIAL EXPRESSION 

test.result <- LPEseq.test(counts.norm[, 'PA14'], counts.norm[, 'ESP055'])
head(test.result, 3)


#read table with strain info, ensuring to eliminating the row containing strain id of ref strain and this will create error message when generating our contrasts downstream
strains <- samples$Identifier

for (straintype in strains){
  strain.result <- LPEseq.test(counts.norm[, 'PA14'], counts.norm[, straintype])
  write.csv(strain.result, paste('DGE/',straintype,'.csv'))
}
