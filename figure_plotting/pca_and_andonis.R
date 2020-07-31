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



# CONTROL FOR GENOTYPE  ---------------------------------------------------

geno_cont <- removeBatchEffect(rld.df, batch=coldata$colony)

#get simplified to early (initial and 3 days) and late (14 days and final set)
coldata = coldata %>% 
  mutate(long_term_time = if_else(timepointDescription %in% c("initial", "3 days of heat"),
                                  'early',
                                  'late'))

#rebuild pca
tag_pca_df = build_pca(geno_cont, coldata,
                       ntop = NTOP,
                       pcs = 2)

#plot with timepoint types color coded
tag_pca_df %>% 
  mutate(timepointTemp = if_else(is.na(timepointTemp),
                                 'unknown',
                                 timepointTemp)) %>% 
  ggplot(aes(x=PC1, y=PC2, color=long_term_time, shape=timepointDescription)) +
  geom_point(size=4)



#CONTROL FOR GENOTYPE AND LONG TERM TIMEPOINT
geno_ltime_cont <- removeBatchEffect(rld.df, batch=coldata$colony, batch2 = coldata$timepointDescription)

#rebuild pca
tag_pca_df = build_pca(geno_ltime_cont, coldata,
                       ntop = NTOP,
                       pcs = 2)

#plot with timepointTemp color coded
tag_pca_df %>% 
  mutate(timepointTemp = if_else(is.na(timepointTemp),
                                 'unknown',
                                 timepointTemp)) %>% 
  ggplot(aes(x=PC1, y=PC2, color=timepointTemp)) +
  geom_point(size=4)

#this doesn't clarify things any :(


# LOOK WITHIN EARLY AND LATE TIMES ----------------------------------------

#modify the column names of rld.df to easier subsetting
colnames(rld.df) = 1:19

#subset coldata by the early and late sample times
early = coldata %>% 
  filter(long_term_time == 'early')
late = coldata %>% 
  filter(long_term_time == 'late')


#run for early
pca_for_sub = function(rld.df, coldata_sub){
  sub_rld = rld.df[, coldata_sub$sampleNumber] #subset by long term timepoint
  sub_pca = build_pca(sub_rld, coldata_sub,
                      ntop = NTOP,
                      pcs = 3) #build pca for subset
  return(sub_pca)
}
e_pca = pca_for_sub(geno_cont, early)

#plot early 
e_pca %>% 
  ggplot(aes(x=PC1, y=PC2, color=timepointTemp)) +
  geom_point(size=4)

#plot early without unknowns
e_pca %>% 
  filter(!is.na(timepointTemp)) %>% 
  ggplot(aes(x=PC1, y=PC2, color=timepointTemp)) +
  geom_point(size=4)


#REPEAT FOR LATE
#subset
l_pca = pca_for_sub(geno_cont, late)
#plot
l_pca %>% 
  ggplot(aes(x=PC1, y=PC2, color=timepointTemp)) +
  geom_point(size=4)
#plot without NAs
l_pca %>% 
  filter(!is.na(timepointTemp)) %>% 
  ggplot(aes(x=PC1, y=PC2, color=timepointTemp)) +
  geom_point(size=4)

#mostly right?


