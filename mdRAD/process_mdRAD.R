#process_mdRAD.R

rm(list=ls())
source('my_functions.R')

# READ IN BEDTOOLS COUNTS -------------------------------------------------

bedtoolsFileList=c('mdRAD/data/geneBoundaries.bed.tsv')
names=c('gene')

datList = lapply(bedtoolsFileList, function(x) read_bedtools(x))
posList = lapply(datList, function(x) return(x[['pos']]))
lengthList = lapply(posList, function(x) return(x$end-x$start))
countsList0 = lapply(datList, function(x) return(x[['counts']]))
countsList = lapply(countsList0, function(x) sort_counts(x))
names(countsList)=names
names(posList)=names
names(lengthList)=names

#TEMP FIX FOR NAMES 
#this is here now to handle the plate rotation, won't be necessary if matrix comes from renamed fastq files
#see troubleshoot_ids.R
ll=load('troubleshoot_ids/name_rotations.Rdata')
for (n in names(countsList)){
  old_dat = countsList[[n]]
  colnames(old_dat) = rs_dat$rot_sampleID
  new_dat = old_dat[,as.character(rs_dat$sampleID)]
  countsList[[n]] = new_dat
}

# SET UP COLDATA ----------------------------------------------------------

#revise names in counts table
counts=countsList[[1]]
sample = sapply(colnames(counts), function(x) strsplit(x, '_')[[1]][1])
snum = sub('mdRAD.', '', sample, fixed = TRUE)
data.frame(colnames = colnames(counts),
           sample = sample,
           snum = snum)

#load sample data
coldata = read_csv('metadata/sample_information_table.csv') #this is the correct one

#check alignment
sum(coldata$sampleNumber==as.numeric(snum))==nrow(coldata)


# GET FPKM -----------------------------------------------------

#run for entire set, and each enzyme individually
fpkmList = map(names, function(x) get_fpkm(x))
names(fpkmList)=names


#output the GBM fpkms for pca_and_adonis.R
geneFpkm = fpkmList[['gene']]
save(geneFpkm, coldata, file='mdRAD/data/mdRAD_geneFPKMs.Rdata')


#get mean fkpm accross each dataset
lvlList = lapply(fpkmList, function(x) get_means(x))
names(lvlList)=names


#re-append positions
lapply(lvlList, head)
for (n in names){
  print(n)
  lvlList[[n]] = cbind(posList[[n]], lvlList[[n]])
  # lvlList[[n]]$tag<-NULL
}
lapply(lvlList, head)

save(lvlList, file='mdRAD/data/methLevelList.Rdata')

#CHECK DISTRIBUTIONS
pltList = list()
for (n in names){
  ldf = lvlList[[n]]
  plt=ldf %>% 
    ggplot(aes(x=mn)) +
    geom_histogram() +
    labs(subtitle=n)
  pltList[[n]]=plt
}
plot_grid(plotlist=pltList, nrow=1)

#CHECK CORRELATION WITH BENCHMARKING
ll=load('mdRAD/data/benchmarking_gbm_lvls.Rdata')
head(gbm.dat)
new_gbm = lvlList['gene'] %>% 
  data.frame() %>% 
  as_tibble()
colnames(new_gbm) = c('chr', 'start', 'end', 'name', 'tag', 'new_fpkm')
sum(new_gbm$name %in% gbm.dat$name)
mdat = new_gbm %>% 
  left_join(gbm.dat) %>% 
  dplyr::select(name, new_fpkm, l.fracMeth, mrB, mrF, mrM, mbd.score)

#build scatterplots
lmdat = mdat %>% 
  pivot_longer(l.fracMeth:mbd.score, 
               names_to = 'assay',
               values_to = 'gbm_measure')
lmdat %>% 
  ggplot(aes(x=gbm_measure, y=new_fpkm)) +
  geom_point() +
  geom_smooth(method='lm', se=FALSE) +
  facet_grid(~assay, scales='free')







#BELOW STILL NEEDS TO BE WORKED OUT







# GET RESPONSES -----------------------------------------------------------


#for tissue
CONTRAST=c('tissue', 't', 's')
t.resList = lapply(countsList, function(x) get_response(x, coldata, CONTRAST))
names(t.resList)=names


#for genotype
CONTRAST=c('genotype', 'L5', 'N12')
g.resList = lapply(countsList, function(x) get_response(x, coldata, CONTRAST))
names(g.resList)=names


#CHECK VOLCANOS

#for tissue
tpltList=list()
for (n in names){
  print(n)
  rdf = data.frame(t.resList[[n]])
  plt=plot_volcano_general(rdf)
  tpltList[[n]]=plt
}
plot_grid(plotlist=tpltList, nrow=3)

#for genotypes
gpltList=list()
for (n in names){
  print(n)
  rdf = data.frame(g.resList[[n]])
  plt=plot_volcano_general(rdf)
  gpltList[[n]]=plt
}
plot_grid(plotlist=gpltList, nrow=3)



#RE-APPEND POSITIONS

#cehck names match up
lapply(t.resList, head)
lapply(t.resList, nrow)
lapply(posList, nrow)

for (n in names){
  pos=posList[[n]]
  tag = paste(pos$name, paste(paste(pos$chr, pos$start, sep='_'), pos$end, sep='_'), sep='_')
  match=sum(rownames(t.resList[[n]])==tag)==nrow(t.resList[[n]])
  print(paste(n,match,sep=' = '))
}


#RE-APPEND POSITIONS
#re-append
mr.t.list = reappend_positions_to_res(t.resList, names)
mr.g.list = reappend_positions_to_res(g.resList, names)
lapply(mr.t.list, head)
lapply(mr.g.list, head)


# REPEAT RESPONSE ANALYSIS FOR 8-SAMPLE SET -------------------------------

#SET UP NEW COLDATA
es.counts=es.countsList[[1]]
es.sample = sapply(colnames(es.counts), function(x) strsplit(x, '_')[[1]][1])
es.enzyme = sub('mr',
             '',
             sapply(es.sample, function(x) strsplit(x, '.', fixed=TRUE)[[1]][1]))
es.genotype = sapply(es.sample, function(x) strsplit(x, '.', fixed=TRUE)[[1]][2])
es.tissue = substr(sapply(es.sample, function(x) strsplit(x, '.', fixed=TRUE)[[1]][3]),
                1,
                1)


es.coldata = data.frame(Run=es.sample,
                     genotype=es.genotype,
                     tissue=es.tissue,
                     enzyme=es.enzyme)
es.coldata #coldata for the 8-sample reduced set
nrow(es.coldata)

#GET RESPONSES

#for tissue
CONTRAST=c('tissue', 't', 's')
es.t.resList = lapply(es.countsList, function(x) get_response_singleEnzyme(x, es.coldata, CONTRAST))
names(es.t.resList)=es.names


#for genotype
CONTRAST=c('genotype', 'L5', 'N12')
es.g.resList = lapply(es.countsList, function(x) get_response_singleEnzyme(x, es.coldata, CONTRAST))
names(es.g.resList)=es.names

#RE-APPEND POSITIONS

#cehck names match up
lapply(es.t.resList, head)
lapply(es.t.resList, nrow)
lapply(es.posList, nrow)

for (n in es.names){
  pos=es.posList[[n]]
  tag = paste(pos$name, paste(paste(pos$chr, pos$start, sep='_'), pos$end, sep='_'), sep='_')
  match=sum(rownames(es.t.resList[[n]])==tag)==nrow(es.t.resList[[n]])
  print(paste(n,match,sep=' = '))
}


#RE-APPEND POSITIONS
#re-append
es.mr.t.list = reappend_positions_to_res(es.t.resList, es.names)
es.mr.g.list = reappend_positions_to_res(es.g.resList, es.names)
lapply(es.mr.t.list, head)
lapply(es.mr.g.list, head)

#check volcano
gdat = data.frame(es.mr.g.list[['gene']])
x=plot_volcano_general(gdat)
x


#SAVE THE 8-SAMPLE REDUCED RESULTS
save(es.mr.t.list, file='methylRAD/datasets/eightSample_tissue_responses.Rdata')
save(es.mr.g.list, file='methylRAD/datasets/eightSample_genotype_responses.Rdata')



# ADD RESPONSE RESULTS FOR SMALL WINDOWS ----------------------------------
#PROBABLY NOT GOING TO USE THESE, SINCE THEY ARE SO CUMBERSOME BUT KEEPING FOR REFERENCE
filePath='methylRAD/datasets/mr_500bp_window_genotype_allResponse.tsv'
t500 = upload_response_to_add('methylRAD/datasets/mr_500bp_window_tissue_allResponse.tsv')
g500 = upload_response_to_add('methylRAD/datasets/mr_500bp_window_genotype_allResponse.tsv')
mr.t.list[['500bp']] = t500
mr.g.list[['500bp']] = g500



#SAVE
save(mr.t.list, file='methylRAD/datasets/tissue_responses.Rdata')
save(mr.g.list, file='methylRAD/datasets/genotype_responses.Rdata')



# GET VSD RESULTS ---------------------------------------------------------

#Didn't end up needing this but keeping for refernce

# #subset for just gene and 1kb
# subCounts = list(countsList[['gene']],
#                  countsList[['1kb']])
# 
# #run vst on each set of counts
# vsdList0 = lapply(subCounts, function(x) get_vsd(x))
# names(vsdList0)=c('gene', '1kb')
# vsdList=vsdList0
# 
# #save
# save(vsdList, coldata, file='methylRAD/datasets/vst_results.Rdata')


