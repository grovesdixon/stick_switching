#assign_genotypes.R

library(vcfR)
library(adegenet)
library(readxl)
source('my_functions.R')

# load genotypes ----------------------------------------------------------

vcfInput = 'troubleshoot_ids/all_filt3_imputed.vcf'
N.CORES = 3
gll=vcfR2genlight(read.vcfR(vcfInput))
print(gll)


# check for missing data --------------------------------------------------

missing = glNA(gll)
has_missing = missing > 0
n_has_missing = sum(has_missing)
tot_loci = ncol(gll)
pct_missing = round((n_has_missing / tot_loci), 3)*100
pct_string = paste('(', pct_missing, '%)', sep='')
good = which(!has_missing)
if (n_has_missing > 0){
  print(paste('WARNING.', n_has_missing, 'of', tot_loci, pct_string, 'loci have at least one NA and will be removed before PCA'))
  gll = gll[,good]
  print('Summary of new genlight object:')
  print(gll)
}


#run PCA
pca=glPca(gll, nf=2, n.cores=N.CORES)
p=pca$scores


# merge with trait labels -------------------------------------------------

#here use the original sample spreadsheet for illustration
#the one in metadata has been corrected
coldata = read_csv('troubleshoot_ids/original_wrong_sample_information_table.csv')

# plot --------------------------------------------------------------------

geno_pca_df = pca$scores %>% 
  data.frame() %>% 
  rownames_to_column('sample') %>% 
  mutate(sampleNumber = sample_numbers_from_names(sample),
         assay = substr(sample, start=1, stop=5)) %>% 
  left_join(coldata) %>% 
  as_tibble()

#function to plot
plot_pca = function(pca_df, pc1, pc2){
  pca_df %>% 
    mutate(assay = substr(sample, start=1, stop=5)) %>% 
    ggplot(aes_string(x=pc1, y=pc2, color = 'colony', shape='assay')) +
    geom_point(size=4)
}

gplt = plot_pca(geno_pca_df, 'PC1', 'PC2')

#get stats on clusters
mod_df = geno_pca_df %>% 
  mutate(gi_colony = if_else(PC2 > 5,
                             'L1',
                             'unassigned'),
         gi_colony = if_else(PC2 < 5 & PC1 < 0,
                             'N4',
                             gi_colony),
         gi_colony = if_else(PC2 < 5 & PC1 > 0,
                             'N1',
                             gi_colony)) %>% 
  dplyr::rename(geno_PC1 = PC1,
         geno_PC2 = PC2)

#check inferred assignments
gi_plt = mod_df %>% 
  ggplot(aes(x=geno_PC1, y=geno_PC2, color=gi_colony)) +
  geom_point(size=4) +
  labs(color = 'inferred colony')


plot_grid(gplt, gi_plt)

#look at mdRAD mismatches
mod_df %>% 
  filter(assay=='mdRAD',
         colony != gi_colony) %>% 
  arrange(sampleNumber) %>% 
  tail


# look at mdRAD results ----------------------------------------------------

library(DESeq2)
#load fpkm data for gbm
ll=load('troubleshoot_ids/mdRAD_geneFPKMs.Rdata')
ll

#build pca
NTOP = round(nrow(geneFpkm) / 5, 0)
md_pca_df = build_pca(geneFpkm, coldata,
                   ntop = NTOP,
                   pcs = 2)

#plot for colony
ori_md_plt = md_pca_df %>% 
  ggplot(aes(x=PC1, y=PC2, color=colony)) +
  geom_point(size=5)
ori_md_plt

#view mismatches
md_N4_mm0 = md_pca_df %>% 
  filter(PC1 < 0,
         colony != 'N1')

md_N1_mm0 = md_pca_df %>% 
  filter(PC1 > 0,
         PC2 < 0,
         colony != 'N4')

md_L1_mm0 = md_pca_df %>% 
  filter(PC1 > 0,
         PC2 > 0,
         colony != 'L1')

md_pca_df %>% 
  filter(sampleNumber %in% c(1,96))




# hypothesize plate rotation ----------------------------------------------

#function to convert a column from the coldata into 96-well plate format based on how they were prepared
column_to_plate = function(coldata, column_string){
  m=coldata %>% 
    dplyr::select(c('sampleNumber', column_string)) %>% 
    arrange(sampleNumber) %>% 
    pull(column_string) %>% 
    matrix(nrow=8, ncol=12)
  colnames(m) = paste('c', 1:12, sep='')
  rownames(m) = LETTERS[1:8]
  m %>% 
    as_tibble()
}

#make genotype and samle number 96-well plate
gplate=column_to_plate(coldata, 'colony')
splate = column_to_plate(coldata, 'sampleNumber')

#double-check against the manually built one
gplate_manual = read_excel('troubleshoot_ids/genotype_plate.xlsx') %>% 
  column_to_rownames('row')
sum(gplate==gplate_manual)

#function to rotate a palte 180 degrees
rotate_plate = function(plate){
  rot_plate = matrix(nrow=nrow(plate),
                     ncol=ncol(plate))
  colnames(rot_plate) = colnames(plate)
  rownames(rot_plate) = rownames(plate)
  rev_plate = rev(data.frame(plate))
  for (i in 1:ncol(rot_plate)){
    rot_plate[,i] = rev(rev_plate[,i])
  }
  return(rot_plate)
}

#rotate the plates
r_gplate = rotate_plate(gplate)
r_splate = rotate_plate(splate)

#write them out for reference
r_gplate %>% 
  data.frame() %>% 
  write_csv(path='troubleshoot_ids/rotated_genotype_plate.csv')
r_splate %>% 
  data.frame() %>% 
  write_csv(path='troubleshoot_ids/rotated_sampleNumber_plate.csv')


#function to combine two 96-well plate back into dataframe
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

#combine the rotated genotypes with the original sample names
rot_plate_dat = combine_plates(list(splate, r_splate, r_gplate),
                               c('sampleNumber', 'rot_sampleNumber', 'rot_colony'))
rot_plate_dat

#merge back with trait data
md_pca_df = md_pca_df %>% 
  left_join(rot_plate_dat, by = 'sampleNumber')

#plot original and roated genotypes
rot_md_plt = md_pca_df %>% 
  ggplot(aes(x=PC1, y=PC2, color = rot_colony)) +
  geom_point(size=4)
plot_grid(ori_md_plt, rot_md_plt + labs(color='rotated colony'))

#assemble the genotype mismatches
md_N4_mm = md_pca_df %>% 
  filter(PC1 < 0,
         rot_colony != 'N4')

md_N1_mm = md_pca_df %>% 
  filter(PC1 > 0,
         PC2 < 0,
         rot_colony != 'N1')
md_mm = rbind(md_N4_mm,
              md_N1_mm)

#count the mismatches
table(md_mm$colony)
table(md_mm$rot_colony)


# look at the tagseq data ----------------------------------------------------

#load
ll = load('troubleshoot_ids/rld.Rdata')
ll

#build pca
rld.df = data.frame(assay(rld))
NTOP = round(nrow(rld.df) / 5, 0)
tag_pca_df = build_pca(rld.df, coldata,
                   ntop = NTOP,
                   pcs = 2)

ori_tag_plt = tag_pca_df %>% 
  ggplot(aes(x=PC1, y=PC2, color=colony)) +
  geom_point(size=4)


#assemble the genotype mismatches
tag_N4_mm = tag_pca_df %>% 
  filter(PC1 < 0,
         PC2 > 0,
         colony != 'N4')

tag_N1_mm = tag_pca_df %>% 
  filter(PC2 < 0,
         !colony %in% c('N1', 'L1'))
tag_mm = rbind(tag_N4_mm,
               tag_N1_mm)

#count the mismatches
table(tag_mm$colony)


# compare mdRAD and tagseq mismatches -------------------------------------

#the same sample number are mismached for genotype
md_nums = md_mm$rot_sampleNumber
tag_nums = tag_mm %>% 
  arrange(sampleNumber) %>% 
  pull(sampleNumber)
sum(md_nums %in% tag_nums) == length(md_nums)
order(md_nums)



# upload new correated sample traits --------------------------------------
#(this one matches the one in metadata/)
#noticed mistake when entering in data from the 
#field notebook for the exact samples identified 

fixed_coldata = read_csv('metadata/sample_information_table.csv')

#merge the mdRAD data with fixed traits using the ROATED sample numbers

#re-upload the mdRAD data
ll=load('troubleshoot_ids/mdRAD_geneFPKMs.Rdata')
ll

#organize ids and numbers from counts matrix for fixing
samples = colnames(geneFpkm)
numbers0 = sapply(samples, function(x) strsplit(x, '_')[[1]][1])
numbers = as.numeric(sub('mdRAD.', '', numbers0, fixed=TRUE))
sdat = data.frame(sampleID = samples,
                  sampleNumber = numbers)

#get the rotated numbers
rs_dat = md_pca_df %>% 
  dplyr::select(sampleNumber, rot_sampleNumber) %>% 
  left_join(sdat, by = 'sampleNumber') %>% 
  as_tibble()

#pull rotated sample IDs
rot_sampleID = rs_dat %>% 
  arrange(rot_sampleNumber) %>% 
  pull(sampleID) %>% 
  as.character()
rs_dat$rot_sampleID = rot_sampleID
rs_dat

#rotate the fpkm matrix 
rot_geneFpkm = geneFpkm[,rot_sampleID]
colnames(rot_geneFpkm) = rs_dat$sampleID

#rebuild the pca
NTOP = round(nrow(rot_geneFpkm) / 5, 0)
fixed_md_pca_df = build_pca(rot_geneFpkm, fixed_coldata,
                      ntop = NTOP,
                      pcs = 2)
fixed_md_plt = fixed_md_pca_df %>% 
  ggplot(aes(x=PC1, y = PC2, color=colony)) +
  geom_point(size=4)

plot_grid(ori_md_plt,
          fixed_md_plt)
ori_md_plt



# replot the tag-seq data with fixed coldata ------------------------------

#build pca
NTOP = round(nrow(rld.df) / 5, 0)
fixed_tag_pca_df = build_pca(rld.df, fixed_coldata,
                       ntop = NTOP,
                       pcs = 2)

fixed_tag_plt = fixed_tag_pca_df %>% 
  ggplot(aes(x=PC1, y=PC2, color=colony)) +
  geom_point(size=4)



# replot the genotypes with fixed info ------------------------------------

#remerge genotype data with fixed coldata
fixed_geno_pca_df =  pca$scores %>% 
  data.frame() %>% 
  rownames_to_column('sample') %>% 
  mutate(sampleNumber = sample_numbers_from_names(sample),
         assay = substr(sample, start=1, stop=5)) %>% 
  left_join(fixed_coldata) %>% 
  as_tibble()


#pull the genotype data from mdRAD
fixed_md_geno0 = fixed_geno_pca_df %>% 
  filter(assay=='mdRAD') %>% 
  dplyr::select(sampleNumber, colony, PC1, PC2, assay) %>% 
  left_join(rs_dat)

#set colony ids for rotated sample Numbers
fixed_md_colonies = fixed_md_geno0$colony
names(fixed_md_colonies) = as.character(fixed_md_geno0$rot_sampleNumber)

#build fixed df
fixed_md_geno = fixed_md_geno0 %>% 
  mutate(colony = fixed_md_colonies[as.character(sampleNumber)]) %>% 
  dplyr::select(sampleNumber, colony, PC1, PC2, assay)


#fix for tagseq
fixed_tag_geno = fixed_geno_pca_df %>% 
  filter(assay=='tagSe') %>% 
  dplyr::select(sampleNumber, colony, PC1, PC2, assay)

#bind them and plot
fixed_geno = rbind(fixed_md_geno, fixed_tag_geno)
fixed_gplt = fixed_geno %>% 
  ggplot(aes(x=PC1, y=PC2, color=colony, shape = assay)) +
  geom_point(size=4)



# plot to show all fixes --------------------------------------------------

pans = plot_grid(gplt,
          fixed_gplt,
          ori_md_plt,
          fixed_md_plt,
          ori_tag_plt,
          fixed_tag_plt,
          nrow=3)
col1 = ggdraw() + draw_label('original', fontface ='bold')
col2 = ggdraw() + draw_label('fixed', fontface ='bold')
row1 = ggdraw() + draw_label('SNPs', angle=90, fontface ='bold')
row2 = ggdraw() + draw_label('mdRAD', angle=90, fontface ='bold')
row3 = ggdraw() + draw_label('tag-seq', angle=90, fontface ='bold')
cols = plot_grid(col1, col2, nrow=1)
rows = plot_grid(row1, row2, row3, nrow=3)
right=plot_grid(cols, pans, nrow=2, rel_heights = c(1,20))
plot_grid(rows, right, rel_widths = c(1,20))


# finally, build commands to correct the mdRAD file names -----------------

rs_dat
ori = sub('.', '-', rs_dat$sampleID, fixed=TRUE)
rot = sub('.', '-', rs_dat$rot_sampleID, fixed=TRUE)

copy_commands = c('mkdir fixed')
for (i in 1:length(ori)){
  original = ori[i]
  new = rot[i]
  clane1 = paste('cp ', original, '_L001_R1_001.fastq fixed/', new, '_L001_R1_001.fastq', sep='')
  clane2 = paste('cp ', original, '_L002_R1_001.fastq fixed/', new, '_L001_R1_002.fastq', sep='')
  copy_commands = append(copy_commands, clane1)
  copy_commands = append(copy_commands, clane2)
}

#write out
fileConn<-file('troubleshoot_ids/rotate_mdRAD_fastq_commands.txt')
writeLines(copy_commands, fileConn)
close(fileConn)

#save name switches for temp fixes
save(rs_dat, file='troubleshoot_ids/name_rotations.Rdata')

