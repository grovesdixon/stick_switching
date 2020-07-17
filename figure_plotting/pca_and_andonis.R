#pca_and_adonis.R
#plot PCAs and adonis partitioning of variance piecharts 

rm(list=ls())
library(vegan)
library(gtools)
library(DESeq2)
library(limma)
source('my_functions.R')


# mdRAD ---------------------------------------------------------------

#load fpkm data for gbm
ll=load('mdRAD/data/mdRAD_geneFPKMs.Rdata')
ll


#build pca
NTOP = round(nrow(geneFpkm) / 5, 0)
pca_df = build_pca(geneFpkm, coldata,
          ntop = NTOP,
          pcs = 10)

#plot for colony
pca_df %>% 
  ggplot(aes(x=PC1, y=PC2, color=colony)) +
  geom_point(size=5)



# PCA FOR TAG-SEQ ---------------------------------------------------------


#load
ll = load('tagSeq/data/rld.Rdata')
ll

#build pca
rld.df = data.frame(assay(rld))
NTOP = round(nrow(rld.df) / 5, 0)
tag_pca_df = build_pca(rld.df, coldata,
                       ntop = NTOP,
                       pcs = 2)

tag_plt = tag_pca_df %>% 
  ggplot(aes(x=PC1, y=PC2, color=colony)) +
  geom_point(size=4)
tag_plt


