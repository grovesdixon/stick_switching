#initialize_counts.R
#prepare reads for DESeq and get gist of data

rm(list=ls())
source('my_functions.R')


rnaDir='tagSeq/data/'

#upload read counts features counts
countsFile=paste(rnaDir, 'feature_counts_out.tsv', sep='')
counts = read.table(countsFile, header = T, row.names='Geneid')
counts = counts[,6:ncol(counts)]
revisedNames = sub('_dupsRemoved.bam', '', colnames(counts))
colnames(counts) = revisedNames
counts = sort_counts(counts, subout='tagSeq.')

#look at counts matrix
head(counts)
dim(counts)
length(unique(colnames(counts))) #96 samples
nrow(counts) #28188 genes
tots = apply(counts, 2, sum)
print(paste("Mean read count per sample =", paste(round(mean(tots) / 1e6, 2), "million reads")))


#remove genes with low coverage
cut=3
cc=counts
means=apply(cc,1,mean) #get the average count for each gene using apply()
table(means>cut)       #look at how many genes have average count greater than the chosen 'cut'
counts=cc[means>cut,]  #reduce the counts matrix to only those genes
nrow(counts)           #show the number of genes with enough coverage to keep



#SET UP COLDATA
#read in
coldata = read_csv('metadata/sample_information_table.csv') #this is the corrected version

#check agreement (make sure the sample numbers in the coldata table match those in the counts column names)
samples = sub('tagSeq.', '', colnames(counts), fixed=TRUE)
snums = sapply(samples, function(x) strsplit(x, '_')[[1]][1]) %>% 
  as.numeric()
sum(snums==coldata$sampleNumber)==nrow(coldata) #do they match?


#------- GET RAW VARIANCE STABILIZED COUNTS ------------#
#use DESEq2 to get a matrix of variance stabilized counts
#this is similar to getting an FPKM matrix
#set up input matrix for DESeq
ddsHTSeq<-DESeqDataSetFromMatrix(counts,
	colData = coldata,
	design = formula(~1))

#run DESeq
dds = DESeq(ddsHTSeq)

#get DEseq results
res = results(dds)

#get variance stabilized counts
rld = vst(dds)
rld.df=assay(rld)
colnames(rld.df) = colnames(counts)



#=====================================================================================
#
#  Code chunk 2
# transpose the dataset you have samples as rows and genes as columns
#=====================================================================================
#transpose the rld matrix
datExpr0 = as.data.frame(t(rld.df));
dim(rld.df)
dim(datExpr0)

#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================

#check that the dataset doesn't have geneswith too many missing values
#these would likely represent lowly expressed genes and under sequenced samples
library(WGCNA)
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK



#=====================================================================================
#
#  Code chunk 4

#=====================================================================================
#removing genes that were flagged with too many missing values
#note how many genes we have right now
before = ncol(datExpr0)
print(before)


if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
     printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
     printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
rld.df=t(datExpr0)
rld=rld[rownames(rld.df),]
dim(rld.df)
dim(rld)
nrow(datExpr0)
after = ncol(datExpr0)
print(paste(before - after, "Genes With Too Many Missing Values Were Removed"))

#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================

#build sample heatmaps 
library(pheatmap)
phm=pheatmap(cor(rld.df))

#build pca
NTOP = round(nrow(rld.df) / 5, 0)
pca_df = build_pca(rld.df, coldata,
                   ntop = NTOP,
                   pcs = 10)

#plot pca
addx=3
addy=2
pca_df %>% 
  ggplot(aes(x=PC1, y=PC2, color=colony)) +
  geom_point(size=5)


#now cluster samples based on gene expression to identify outliers
sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,5,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)


#save the outlier names so you can optionally remove them in other analyses
# save(outlierNames, file = 'datasets/outliers.Rdata')
dim(counts)
dim(coldata)

#save the variance stabilized counts
rldFile=paste(rnaDir, 'rld.Rdata', sep='')
save(rld, coldata, file=rldFile)

#save the filtered raw counts for DESeq2 analysis
deseqInfile=paste(rnaDir, 'deseqInput.Rdata', sep='')

#standardize the colnames and rownames for counts and coldata
snum = sapply(colnames(counts), function(x) strsplit(x, '_')[[1]][1])
snum = as.numeric(sub('tagSeq.', '', snum, fixed=TRUE))
colnames(counts) = snum
coldata = data.frame(coldata)
rownames(coldata) = snum
save(counts, coldata, file=deseqInfile)

