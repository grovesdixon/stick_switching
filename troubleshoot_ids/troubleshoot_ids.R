#toubleshoot_ids.R


rm(list=ls())
library(vegan)
library(gtools)
library(DESeq2)
source('my_functions.R')

# COMPARE MISLABELED WITH GSAF UPLOAD -------------------------------------

#upload the barcodes submitted to gsaf
library(readxl)
gsaf = read_excel('metadata/GSAFSampleUpload.xlsx') %>% 
  filter(!grepl('^JP', SampleName)) %>% 
  mutate(num = sub('tagSeq.', '', SampleName, fixed=TRUE),
         num = sub('mdRAD.', '', num, fixed=TRUE),
         num = as.numeric(num)) %>% 
  dplyr::rename(i7 = `Barcode i7 ("standard")`,
                i5 = `Barcode i5 (Dual index)`) %>% 
  dplyr::select(-`Sample Description`)
sum(table(gsaf$num)==2) #two reps of each sample number (one for tagseq one for mdRAD)


#double-check that numbers are always matched
tag = gsaf %>% 
  filter(grepl('^tagSeq', SampleName))
md = gsaf %>% 
  filter(grepl('^mdRAD', SampleName))
match = tag %>% 
  full_join(md, by = 'num')
sum(match$i7.x==match$i7.y)
#good, they did


#upload my barcode sequence-number table
barcodes = read_csv('metadata/oligo_barcode_sequences_rc.csv')

#slice out i7 sequences
i7bcs = barcodes %>% 
  filter(grepl('^ILL', oligo)) %>% 
  set_names(c('oligo', 'for', 'i7')) %>% 
  dplyr::select(oligo, i7)

#check presence in gsaf upload df
i7_seqs = unique(i7bcs$i7)
length(i7_seqs)
length(unique(gsaf$i7))
sum(i7_seqs %in% gsaf$i7)
#these were sequenced with Novaseq, so i7 barcodes were given as reverse compliment
#we used 36 of the 72 barcodes that we have in this study
#all samples have a fastq from GSAF, so we know all barcodes made it in


#splice out i5 sequences ()
i5bcs = barcodes %>% 
  filter(grepl('^TruSeq', oligo)) %>% 
  set_names(c('oligo', 'i5', 'rc')) %>% 
  dplyr::select(oligo, i5)

#check presence in gsaf upload df
i5_seqs = unique(i5bcs$i5)
length(i5_seqs)
length(unique(gsaf$i5))
sum(i5_seqs %in% gsaf$i5)
#i5 barcodes were given in forward direction
#we used all 6 of our i5 (TrueSeqUN) barcodes for this study
#all samples have a fastq from GSAF, so we know all barcodes made it in


#ADD IN THE BARCODE IDS TO GSAF
gsaf2 = gsaf %>% 
  left_join(i7bcs, by = 'i7') %>% 
  dplyr::rename(i7oligo = oligo) %>% 
  left_join(i5bcs, by = 'i5') %>% 
  dplyr::rename(i5oligo = oligo) %>% 
  arrange(num) %>% 
  dplyr::rename(sampleNumber=num)


gsaf_tag = gsaf2 %>% 
  filter(grepl('tagSeq', SampleName))

gsaf_md = gsaf2 %>% 
  filter(grepl('mdRAD', SampleName))



# TEST ROTATION HYPOTHESIS ------------------------------------------------

#upload the tagseq pca dataframe
ll=load('metadata/raw_tagseq_pca_df.Rdata')
ll

#show the mixups
pca_df %>% 
  ggplot(aes(x=PC1, y=PC2, color=colony)) +
  geom_point(size=5)

#idea here is to search for a way to rotate the 96-well plate that fixes things

#read in the plate with genotype labels
gplate = read_excel('troubleshoot_ids/genotype_plate.xlsx') %>% 
  column_to_rownames('row')
#read in teh plate with sample numbers
splate = read_excel('troubleshoot_ids/sample_number_plate.xlsx') %>% 
  column_to_rownames('row')

#these can be built from the data as well
column_to_plate = function(column_string){
  m=coldata %>% 
    dplyr::select(c('sampleNumber', column_string)) %>% 
    arrange(sampleNumber) %>% 
    pull(column_string) %>% 
    matrix(nrow=8, ncol=12)
  colnames(m) = paste('c', 1:12)
  rownames(m) = LETTERS[1:8]
  return(m)
}

#check
gplate_from_dat=column_to_plate('colony')
gplate_from_dat == gplate

#build other types of plates
treat_plate = column_to_plate('treatment')
time_plate = column_to_plate('timepointDescription')


#rotate the splate 180 degrees
rotate_plate = function(plate){
  rot_plate = matrix(nrow=nrow(plate),
                     ncol=ncol(plate))
  colnames(rot_plate) = colnames(plate)
  rownames(rot_plate) = rownames(plate)
  rev_plate = rev(plate)
  for (i in 1:ncol(rot_plate)){
    rot_plate[,i] = rev(rev_plate[,i])
  }
  return(rot_plate)
}
r_splate = rotate_plate(splate)
r_gplate = rotate_plate(gplate)


#pair up plates

plate_list = list(splate, gplate)
colnames = c('plateSampleNumber', 'plateColony')
combine_plates = function(plate_list,
                          colnames){
  coldat_list = list()
  for (i in 1:length(colnames)){
    cname = colnames[i]
    plate = as.matrix(plate_list[[i]])
    coldat = c(plate)
    coldat_list[[cname]] = coldat
  }
  data.frame(coldat_list) %>% 
    as_tibble()
}

ori_plate_dat = combine_plates(list(splate, gplate, treat_plate, time_plate),
                           c('sampleNumber', 'plateColony', 'treat', 'time'))
rot_plate_dat = combine_plates(list(r_splate, gplate),
                               c('sampleNumber', 'rot_plateColony'))
plate_dat = ori_plate_dat %>% 
  left_join(rot_plate_dat, by = 'sampleNumber')
sum(plate_dat$plateColony != plate_dat$rot_plateColony)

plate_dat %>% 
  filter(plateColony != rot_plateColony) %>% 
  data.frame()


mod_pca = pca_df %>% 
  left_join(plate_dat, by='sampleNumber') %>% 
  as_tibble()

ori = mod_pca %>% 
  ggplot(aes(x=PC1, y=PC2, color=colony)) +
  geom_point(size=4)

plate = mod_pca %>% 
  ggplot(aes(x=PC1, y=PC2, color=plateColony)) +
  geom_point(size=4)

rot = mod_pca %>% 
  ggplot(aes(x=PC1, y=PC2, color=rot_plateColony)) +
  geom_point(size=4)

plot_grid(ori, plate, rot, nrow=3)
#so you can rotate the entire plate without messing up very many L1s

#now create an 'isolated' rotationto fix the cluster that's wrong in tagseq (see barcoding_plates_GD_5-1-20)
splate
to_rotate = splate[,c(7:9)]
left = splate[,c(1:6)]
right = splate[,c(10:12)]
rotated0 = rotate_plate(to_rotate)
rotated = cbind(left, rotated0, right)
iso_rot = combine_plates(list(rotated, gplate, treat_plate, time_plate),
                         c('sampleNumber', 'rot_colony', 'rot_treat', 'rot_time'))

mod_pca2 = pca_df %>% 
  left_join(iso_rot, by='sampleNumber') %>% 
  as_tibble()

#plot original
mod_pca2 %>% 
  ggplot(aes(x=PC1, y=PC2, color=colony)) +
  geom_point(size=4)

mod_pca2 %>% 
  ggplot(aes(x=PC1, y=PC2, color=rot_plateColony)) +
  geom_point(size=4)


mod_coldata = coldata %>% 
  left_join(iso_rot, by = 'sampleNumber') %>% 
  rownames_to_column('sampleID') %>% 
  as_tibble() %>% 
  unite('ROT',rot_colony, rot_time, rot_treat, sep='_', remove = FALSE) %>% 
  unite('ORI', colony, timepointDescription, treatment, sep='_', remove = FALSE) %>% 
  mutate(changed = ROT != ORI) %>% 
  dplyr::select(-ROT, -ORI)

sum(mod_coldata$changed)


# EXAMINE EFFECTS OF HYPOTHETICAL ROTATION --------------------------------

#build a mod_coldata file that includes feature columsn resulting from rotations as above
mod_coldata

#upload the rld results from tagseq
ll=load('tagSeq/data/rld.Rdata')
ll
#modify
rld.df = data.frame(assay(rld))
colnames(rld.df) = sapply(colnames(rld.df), function(x) strsplit(x, '_')[[1]][1])

#build pca
NTOP = round(nrow(rld.df) / 5, 0)
pca_df = build_pca(rld.df, mod_coldata,
                   ntop = NTOP,
                   pcs = 10)

#check re-assignment
plot_pca_comparison(rld.df, mod_coldata, COLOR='colony', SHAPE='timepointDescription', ROT_COLOR='rot_colony', ROT_SHAPE='rot_time')


#check colony re-assignments
pca_df %>% 
  ggplot(aes(x=PC1, y=PC2, color=colony)) +
  geom_point(size=5)
pca_df %>% 
  ggplot(aes(x=PC1, y=PC2, color=rot_colony)) +
  geom_point(size=5)


#build PCAs for individual colonies after assignments
############### L1
#subset for L1
geno_string = 'L1'
gcoldata = mod_coldata %>% 
  filter(grepl(geno_string, colony))
gids = paste('tagSeq.', gcoldata$sampleNumber, sep='')
grld = rld.df[,gids]
dim(gcoldata)
dim(grld)

#build pca
NTOP = round(nrow(grld) / 10, 0)
pca_df = build_pca(grld, gcoldata,
                   ntop = NTOP,
                   pcs = 10)
#check colony time reassignments
plot_pca_comparison(grld, gcoldata, COLOR='timepointDescription', SHAPE='treatment', ROT_COLOR='rot_time', ROT_SHAPE='rot_treat')

#SUBSET FOR 14 DAY TIMEPOINT
gt_coldata = mod_coldata %>% 
  filter(grepl(geno_string, colony),
         timepointDescription=='14 days of heat')
gids = paste('tagSeq.', gt_coldata$sampleNumber, sep='')
gt_rld = rld.df[,gids]
plot_pca_comparison(gt_rld, gt_coldata, COLOR='timepointTemp', SHAPE='treatment', ROT_COLOR='timepointTemp', ROT_SHAPE='rot_treat')


############## CONTROL FOR GENOTYPE, THEN CHECK FOR DOMINANT VARIATION
#control for genotype
library(limma)
cont <- removeBatchEffect(rld.df, batch=mod_coldata$rot_colony)



#plot
plot_pca_comparison(cont, mod_coldata)


#get rid of the initial timepoint, and repeat
noi_coldata = mod_coldata %>% 
  filter(timepointDescription != 'initial')
gids = paste('tagSeq.', noi_coldata$sampleNumber, sep='')
noi_cont = cont[,gids]
plot_pca_comparison(noi_cont, noi_coldata)


#get rid of the all but middle timepoints
noi_coldata = mod_coldata %>% 
  filter(!timepointDescription %in% c('initial', 'final meth experiment fix'))
gids = paste('tagSeq.', noi_coldata$sampleNumber, sep='')
noi_cont = cont[,gids]
plot_pca_comparison(noi_cont, noi_coldata)

sum(mod_coldata$timepointDescription==mod_coldata$rot_time)


########## CONTROL FOR GENOTYPE AND TIME THEN LOOK AT RE-ASSIGNMENTS
library(limma)
dcont <- removeBatchEffect(rld.df,
                          batch=mod_coldata$rot_colony,
                          batch2=mod_coldata$rot_time)

#plot
plot_pca_comparison(dcont, mod_coldata,COLOR='timepointTemp', SHAPE='treatment', ROT_COLOR='rot_treat', ROT_SHAPE='rot_treat')



# subset for 14 hours -----------------------------------------------------

ts_coldata = mod_coldata %>% 
  filter(timepointDescription=="14 days of heat")
gids = paste('tagSeq.', ts_coldata$sampleNumber, sep='')
ts_rld = rld.df[,gids]
ts_cont = removeBatchEffect(ts_rld,
                            batch=ts_coldata$rot_colony)
plot_pca_comparison(ts_cont, ts_coldata,COLOR='timepointTemp', SHAPE='treatment', ROT_COLOR='rot_treat', ROT_SHAPE='rot_treat')



# TAGSEQ ------------------------------------------------------------------
#gather mislabeled samples for tagseq

ll=load('metadata/raw_tagseq_pca_df.Rdata')
ll


pca_df %>% 
  ggplot(aes(x=PC1, y=PC2, color=colony)) +
  geom_point(size=5)


#identify mismatch
n4_ts_mislab = pca_df %>% 
  filter(PC2 > 0,
         PC1 < 0,
         colony !='N4')  %>% 
  mutate(misAssay = 'ts',
         misGeno = 'N4')

#identify mismatch
n1_ts_mislab = pca_df %>% 
  filter(PC2 < 0,
         colony !='L1',
         colony !='N1') %>% 
  mutate(misAssay = 'ts',
         misGeno = 'N1')

#merge tagseq misses with sample numbers to check for pairings
ts_misses = rbind(n4_ts_mislab  %>% 
                    dplyr::select(sampleNumber, colony, age, treatment, misGeno) %>% 
                    dplyr::rename(colonyLabel = colony,
                                  colonyPCA = misGeno) %>% 
                    left_join(gsaf_tag, by = 'sampleNumber'),
                  n1_ts_mislab  %>% 
                    dplyr::select(sampleNumber, colony, age, treatment, misGeno) %>% 
                    dplyr::rename(colonyLabel = colony,
                                  colonyPCA = misGeno) %>% 
                    left_join(gsaf_tag, by = 'sampleNumber'))

#check treatment counts
ts_misses %>% 
  group_by(colonyPCA, treatment) %>% 
  summarize(N=n())
#so they line up

table(ts_misses$i7oligo)
table(ts_misses$treatment)
table(ts_misses$colonyPCA)

#build a pairing table
n1pair = n1_ts_mislab %>% 
  dplyr::select(sampleNumber, colonyLabel, age, treatment)
n4pair = n4_ts_mislab %>% 
  dplyr::select(sampleNumber, colonyLabel, age, treatment)
cbind(n1pair, n4pair)

# BUILD A GENOTYPE 96 WELL PLATE ------------------------------------------

num_col = gsaf_tag %>% 
  left_join(coldata, by = 'sampleNumber') %>% 
  arrange(sampleNumber) %>% 
  dplyr::select(sampleNumber,
                colony)

nums = num_col$sampleNumber
genos = num_col$colony

matrix(nums, nrow = 8, ncol=12)
gm = matrix(genos, nrow=8, ncol=12)
colnames(gm) = 1:12
rownames(gm) = LETTERS[1:8]
write.table(gm, file='troubleshoot_ids/genotype96well.tsv', quote=FALSE)

# mdRAD ---------------------------------------------------------------
#gather mislabeled samples for mdRAD

#load fpkm data for gbm
ll=load('mdRAD/data/mdRAD_geneFPKMs.Rdata')
ll

#build pca
NTOP = round(nrow(geneFpkm) / 10, 0)
pca_df = build_pca(geneFpkm, coldata,
                   ntop = NTOP,
                   pcs = 10)

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
  filter(PC1 < -500,
         colony != 'N1') %>% 
  mutate(misAssay = 'mr',
         misGeno = 'N1')


#check mixups
l1_mr_mislab = pca_df %>% 
  filter(PC1 > 0,
         PC2 > 0,
         colony != 'L1') %>% 
  mutate(misAssay = 'mr',
         misGeno = 'L1')

#check mixups
n4_mr_mislab = pca_df %>% 
  filter(PC1 > 0,
         PC2 < 0,
         colony != 'N4') %>% 
  mutate(misAssay = 'mr',
         misGeno = 'N4')



# ASSEMBLE ALL MISLABELES -------------------------------------------------

miss = list(l1_mr_mislab,
                 n1_mr_mislab,
                 n1_ts_mislab,
                 n4_mr_mislab,
                 n4_ts_mislab) %>% 
  purrr::reduce(rbind) %>% 
  as_tibble() %>% 
  dplyr::select(sampleNumber, colony, age, treatment, timepointTemp, notebookPg, misAssay, misGeno)

miss




#SLICE OUT THE MISLABELED SAMPLE NUMBERS
cdat = coldata

miss_nums = unique(miss$sampleNumber)
length(miss_nums)
mgsaf = gsaf2 %>% 
  filter(sampleNumber %in% miss$sampleNumber) %>% 
  left_join(cdat, by = 'sampleNumber') %>% 
  dplyr::select(SampleName:treatment) %>% 
  separate(SampleName, into = c('assay', 'snum'))
table(mgsaf$i7oligo)


# TRY TO RESOLVE FROM TAGSEQ ----------------------------------------------

#could not find any mistakes in the spreadsheets
#now try to use tagseq to figure out how mismaches occured

ll=load('metadata/raw_tagseq_pca_df.Rdata')
ll

#genotype as called before
pca_df %>% 
  ggplot(aes(x=PC1, y=PC2, color=colony)) +
  geom_point(size=4)

#modify coldata for troubleshooting
early = c('initial', '3 days of heat')
tag_miss = miss %>% 
  filter(misAssay=='ts') %>% 
  dplyr::select(sampleNumber, misGeno) %>% 
  rename(genoByData = misGeno)


coldata = coldata %>% 
  left_join(tag_miss, by='sampleNumber') %>% 
  mutate(broad_time = if_else(timepointDescription %in% early,
                              'first3days',
                              'last3days'),
         genoByData = if_else(is.na(genoByData),
                              colony,
                              genoByData))



# PLOT WITH GENOTYPES FROM PCA --------------------------------------------

NTOP = round(nrow(rld.df) / 10, 0)
pca_df = build_pca(rld.df, coldata,
                   ntop = NTOP,
                   pcs = 10)

#original genotype calls
pca_df %>% 
  ggplot(aes(x=PC1, y=PC2, color=colony)) +
  geom_point(size=4)

#fixed genotype calls
pca_df %>% 
  ggplot(aes(x=PC1, y=PC2, color=genoByData)) +
  geom_point(size=4)


#sample timepoint 
pca_df %>% 
  ggplot(aes(x=PC1, y=PC2, color=timepointDescription)) +
  geom_point(size=4)


# SET THINGS UP BY SAMPLE NUMBER FOR SIMPLICITY ---------------------------

#rld
rld_bynum = rld.df
rev1 = sapply(colnames(rld.df), function(x) strsplit(x, '_')[[1]][1])
num_str = sub('tagSeq.', '', rev1, fixed=TRUE)
colnames(rld_bynum) = num_str

#coldata
coldata$num_str = as.character(coldata$sampleNumber)



# USE L1 TO SEE WHAT VARIATION SHOULD LOOK LIKE ---------------------------



geno_coldata = coldata %>% 
  filter(colony=='L1')
geno_rld = rld_bynum[,geno_coldata$num_str]

geno_pca = build_pca(geno_rld, geno_coldata, ntop=NTOP)

#with all timepoints
geno_pca %>% 
  ggplot(aes(x=PC1, y=PC2, color = timepointDescription, shape=treatment)) +
  geom_point(size=4) +
  labs(subtitle='L1 colony')

#with broad timepoints
geno_pca %>% 
  ggplot(aes(x=PC1, y=PC2, color = broad_time, shape=age)) +
  geom_point(size=4)

#CHECK WITHIN TIMEPOINTS
NTOP=10000
unique(geno_coldata$timepointDescription)

#check 3 days
tp="3 days of heat"
sub_coldata = geno_coldata %>% 
  filter(timepointDescription==tp)
sub_rld = rld_bynum[,sub_coldata$num_str]
sub_pca = build_pca(sub_rld, sub_coldata, ntop=NTOP)
sub_pca %>% 
  ggplot(aes(x=PC1, y=PC2, color=timepointTemp, shape=treatment)) +
  geom_point(size=4)


tp="14 days of heat"
sub_coldata = geno_coldata %>% 
  filter(timepointDescription==tp)
sub_rld = rld_bynum[,sub_coldata$num_str]
sub_pca = build_pca(sub_rld, sub_coldata, ntop=NTOP)
sub_pca %>% 
  ggplot(aes(x=PC1, y=PC2, color=timepointTemp, shape=treatment)) +
  geom_point(size=4)


# CHECK WITH N1 as called -----------------------------------------------------------

geno_coldata = coldata %>% 
  filter(colony=='N1')
geno_rld = rld.df[,geno_coldata$sampleNumber]

geno_pca = build_pca(geno_rld, geno_coldata, ntop=NTOP)

#with broad timepoints
geno_pca %>% 
  ggplot(aes(x=PC1, y=PC2, color = broad_time, shape=genoByData)) +
  geom_point(size=4)


# CHECK WITH N1 by data -----------------------------------------------------------

geno_coldata = coldata %>% 
  filter(genoByData=='N1')
geno_rld = rld.df[,geno_coldata$sampleNumber]

geno_pca = build_pca(geno_rld, geno_coldata, ntop=NTOP)


#with all timepoints
geno_pca %>% 
  ggplot(aes(x=PC1, y=PC2, color = timepointDescription, shape=treatment)) +
  geom_point(size=4)

#with broad timepoints
geno_pca %>% 
  ggplot(aes(x=PC1, y=PC2, color = broad_time, shape=colony)) +
  geom_point(size=4)




# CHECK WITH N4 -----------------------------------------------------------

geno_coldata = coldata %>% 
  filter(colony=='N4')
geno_rld = rld.df[,geno_coldata$sampleNumber]

geno_pca = build_pca(geno_rld, geno_coldata, ntop=NTOP)

#with broad timepoints
geno_pca %>% 
  ggplot(aes(x=PC1, y=PC2, color = broad_time, shape=genoByData)) +
  geom_point(size=4)


# CHECK WITH N4 FROM DATA -------------------------------------------------

geno_coldata = coldata %>% 
  filter(genoByData=='N4')
geno_rld = rld.df[,geno_coldata$sampleNumber]

geno_pca = build_pca(geno_rld, geno_coldata, ntop=NTOP)


#with broad timepoints
geno_pca %>% 
  ggplot(aes(x=PC1, y=PC2, color = broad_time, shape=genoByData)) +
  geom_point(size=4)

#with broad timepoints
geno_pca %>% 
  ggplot(aes(x=PC1, y=PC2, color = broad_time, shape=genoByData)) +
  geom_point(size=4)








#build pca
NTOP = round(nrow(geneFpkm) / 5, 0)
pca_df = build_pca(geneFpkm, coldata,
                   ntop = NTOP,
                   pcs = 10)

#look into n1 based on the data
n1_bydat = pca_df %>% 
  rownames_to_column('sampleID') %>% 
  filter(PC1 < 0 & PC2 > 0)


geno_rld = rld.df[,n1_bydat$sampleID]
geno_coldata = coldata[n1_bydat$sampleNumber,]
dim(geno_rld)
dim(geno_coldata)
geno_pca = build_pca(geno_rld, geno_coldata,
                   ntop = NTOP,
                   pcs = 10) %>% 
  rownames_to_column('sampleID') %>% 
  as_tibble()

g = geno_pca %>% 
  ggplot(aes(x=PC1, y=PC2, color = colony)) +
  geom_point(size=4)

t = geno_pca %>% 
  ggplot(aes(x=PC1, y=PC2, color = treatment, shape=colony)) +
  geom_point(size=4)

geno_pca %>% 
  filter(colony=='N4') %>% 
  ggplot(aes(x=PC1, y=PC2, color = timepointTemp, shape=colony)) +
  geom_point(size=4)



