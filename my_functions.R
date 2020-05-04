#my_functions.R

library('tidyverse')
library('cowplot')
theme_set(theme_cowplot())

#function to reduce the sample names to just their integer number assignments
sample_numbers_from_names = function(name_vector){
  sample_nums = sub('mdRAD-', '', name_vector, fixed=TRUE)
  sample_nums = sub('tagSeq-', '', sample_nums, fixed=TRUE)
  sample_nums = sapply(sample_nums, function(x) strsplit(x, '_')[[1]][1])
  as.numeric(sample_nums)
}

#read in bedtools results
add_tag_for_pos = function(x){
  x['pos'] %>% 
    unite(tag=name,chr,start,stop)
}

read_bedtools = function(bedtoolsFile, useName=FALSE){
  print(paste('Reading in file', bedtoolsFile))
  counts = read.table(bedtoolsFile, header = TRUE)
  positions = counts[,1:4]
  #make a tag of the name_chr_start_end
  name_chr = paste(counts$name, counts$chr, sep='_')
  start_end = paste(counts$start, counts$end, sep='_')
  tag = paste(name_chr, start_end, sep='_')
  rownames(counts)=tag
  if (useName==TRUE){
    rownames(counts)=counts$name
  }
  counts=counts[,5:ncol(counts)]
    return(list('counts'=counts,
                'pos'=positions))
}

#function specific to the numbering in this experiment
#change order of counts file to fit sample numbering
sort_counts = function(counts_df, subout='mdRAD.'){
  print('Sorting names for counts DF:')
  print(head(counts_df))
  sample = sapply(colnames(counts_df), function(x) strsplit(x, '_')[[1]][1])
  snum = sub(subout, '', sample, fixed = TRUE)
  ordered_cols = data.frame(cname = colnames(counts_df),
                        snum = as.numeric(snum)) %>% 
    arrange(snum) %>% 
    pull(cname) %>% 
    as.character()
  ordered_counts = counts_df[,ordered_cols]
  print('new sorted names:')
  print(colnames(ordered_counts))
  return(ordered_counts)
}



#get fpkm for mdRAD counts
get_fpkm = function(n){
  print(n)
  lengths = lengthList[[n]]
  counts = countsList[[n]]
  m = apply(counts, 2, function(x) sum(x)/1e6)
  k=lengths/1e3
  fpkm=sweep(counts, 1, k, `/`) %>% 
    sweep(2, m, `/`)
  return(fpkm)
}

#funciton to get means for fpkm
get_means = function(df){
  mns=apply(df, 1, mean)
  res=data.frame(tag=names(mns),
                 mn=log(mns, 2))
  return(res)
}

#function to plot pca from normalized expression counts
#requires one df of data, and a coldata with rows matching to df's columns
addx=3
addy=2
build_pca = function(df,
                     coldata,
                     ntop = 25000,
                     pcs = 6){
  #get row varainces and select top
  df = as.matrix(df)
  rv <- rowVars(df)
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  #build pca
  pca <- prcomp(t(df[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  d = cbind(data.frame(pca$x[,1:pcs]),
            coldata)
  attr(d, "percentVar") <- percentVar[1:2]
  return(d)
}


#function to plot original and rotated comparison
plot_pca_comparison = function(rld_df, coldata_df, COLOR='timepointDescription', SHAPE='treatment', ROT_COLOR='rot_time', ROT_SHAPE='rot_treat'){
  NTOP = round(nrow(rld_df) / 10, 0)
  pca_df = build_pca(rld_df, coldata_df,
                     ntop = NTOP,
                     pcs = 10)
  
  cdat = pca_df %>% 
    filter(changed)
  
  ori_plt = pca_df %>% 
    ggplot() +
    geom_point(aes_string(x='PC1', y='PC2', color=COLOR, shape=SHAPE), size=5) +
    geom_point(data=cdat, aes_string(x='PC1', y='PC2', pch=8), color='black', pch=8) +
    labs(subtitle='original')
  rot_plt = pca_df %>% 
    ggplot() +
    geom_point(aes_string(x='PC1', y='PC2', color=ROT_COLOR, shape=ROT_SHAPE), size=5) +
    geom_point(data=cdat, aes_string(x='PC1', y='PC2'), color='black', pch=8) +
    labs(subtitle='rotation')
  plot_grid(ori_plt, rot_plt)
  
}


