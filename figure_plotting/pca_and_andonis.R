#pca_and_adonis.R
#plot PCAs and adonis partitioning of variance piecharts 

rm(list=ls())
library(vegan)
library(gtools)
library(DESeq2)
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

rownames(pca_df)==colnames(geneFpkm)
sample = sapply(colnames(geneFpkm), function(x) strsplit(x, '_')[[1]][1])
snum = sub('mdRAD.', '', sample, fixed = TRUE)
data.frame(colnames = colnames(counts),
           sample = sample,
           snum = snum)
snum == coldata$sampleNumber


#plot for colony
pca_df %>% 
  ggplot(aes(x=PC1, y=PC2, color=colony)) +
  geom_point(size=5)

#check mixups
pca_df %>% 
  filter(PC1 < -500) %>% 
  pull(colony) %>% 
  table()

#check mixups
n1_mr_mislab = pca_df %>% 
  rownames_to_column('sampleID') %>% 
  filter(PC1 < -500,
         colony != 'N1') 


#check mixups
l1_mr_mislab = pca_df %>% 
  filter(PC1 > 0,
         PC2 > 0,
         colony != 'L1')

#check mixups
n4_mr_mislab = pca_df %>% 
  rownames_to_column('sampleID') %>% 
  filter(PC1 > 0,
         PC2 < 0,
         colony != 'N4')



miss = rbind(n1_mr_mislab, n4_mr_mislab)


miss %>% 
  write_csv(path='~/Desktop/mdRAD_missmatches.csv')





#STUFF TO REVISE FOR THIS PROJECT BELOW






#SUBSET FOR THE TWO ENZYMES
fColdata = mrColdata %>% 
  filter(enzyme=='F')
fFpkm = geneFpkm[,fColdata$Run]
mColdata = mrColdata %>% 
  filter(enzyme=='M')
mFpkm = geneFpkm[,mColdata$Run]

#PLOT OVERALL RESULTS
#overall pca
mr.pca = plot_vsd_pca(geneFpkm, mrColdata, shape.col = 'enzyme') + 
  theme(legend.position = 'right') +
  labs(shape='enzyme',
       color='colony')

#overal adonis

ad=adonis(t(geneFpkm)~genotype+tissue+enzyme)
labs=c("genotype","tissue", "enzyme", "residuals")
aovTab = ad$aov.tab
r2s = aovTab$R2[1:4]
ps = aovTab$`Pr(>F)`[1:3]
stars = stars.pval(ps)
labs[1] = paste(labs[1], stars[1], sep='')
labs[2] = paste(labs[2], stars[2], sep='')
labs[3] = paste(labs[3], stars[3], sep='')
cols=append(gg_color_hue(3), 'grey')
# pie(ad$aov.tab$R2[1:3],labels=labs,col=cols,main="WGBS")
pidat = data.frame(r2 = r2s/sum(r2s)*100,
                   p = append(ps, NA),
                   lab=factor(labs, levels=labs))
#plot with ggplot
mr.pie = pidat %>% 
  ggplot() +
  geom_bar(aes(fill=lab, y=r2, x=''), stat='identity', color='black') +
  coord_polar("y", start=.5) +
  scale_fill_manual(values=cols) +
  theme(axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())

#--- PLOT INDIVIDUAL ENZYMES ---#

#FSPE1
#pca
f.pca = plot_vsd_pca(fFpkm, fColdata) + 
  theme(legend.position = 'right') +
  labs(shape='tissue',
       color='colony')

#adonis
f.pie = plot_adonis_from_vsd(fFpkm, fColdata)

#MSPJ1
#pca
m.pca = plot_vsd_pca(mFpkm, mColdata) + 
  theme(legend.position = 'right') +
  labs(shape='tissue',
       color='colony')

#adonis
m.pie = plot_adonis_from_vsd(mFpkm, mColdata)



# ASSEMBLE ALL RESULTS ----------------------------------------------------

#assemble pcas
pcaList = list(pm.pca + 
                 labs(title='WGBS') +
                 theme(legend.position = 'none'),
               mbd.pca + 
                 labs(title='MBD-seq') +
                 theme(legend.position = 'none'),
               f.pca + 
                 labs(title='FspE1') +
                 theme(legend.position = 'none'),
               m.pca +
                 labs(title='MspJ1'))
pcaList = lapply(pcaList, function(x) return(x+theme(plot.title = element_text(hjust=0.5))))


#assemble pie charts
pieList = list(pm.pie + labs(title='WGBS'),
               mbd.pie + labs(title='MBD-seq'),
               f.pie + labs(title='FspE1'),
               m.pie + labs(title='MspJ1'))
pieList = lapply(pieList, function(x) x+theme(plot.title=element_text(hjust=0.5)))


plot_grid(plotlist = pcaList, nrow=1, rel_widths = c(1,1,1,1.1))
plot_grid(plotlist = pieList, nrow=2)


#plot full mdRAD
plot_grid(mr.pca,
          mr.pie)

