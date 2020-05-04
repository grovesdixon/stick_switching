#!/usr/bin/env Rscript
#windows_around_genes_from_gff.R
#This builds a bed file indicating windows upstream and downstream from genes
#output is a bed file with first three columns as required,
#and a final 'name' column with information about each window.
#In the names, windows are numbered relative to the GENE'S start

#Example names for a forward gene:
#w1-Amillepora16770-574-1574-upstream
#w79-Amillepora16770-78574-79574-inside
#w86-Amillepora16770-86389-87389-downstream

#Example names for a forward gene:
# w1-Amillepora16774-268672-269672-upstream
# w101-Amillepora16774-168672-169672-inside
# w141-Amillepora16774-128293-129293-downstream

#PARSE ARUGMENTS
suppressMessages(library(optparse))
suppressMessages(library(tidyverse))

option_list = list(
  
  make_option(c("--i"), type="character", default=NULL, 
              help="infile"),
  
  make_option(c("--gene_id"), type="character", default='ID', 
              help="gene_id"),
  
  make_option(c("--window_size"), type="integer", default=1000, 
              help="size of the windows in bp"),
  
  make_option(c("--buffer"), type="integer", default=50000, 
              help="how far out upstream and downstream to make windows from gene starts"),
  
  make_option(c("--o"), type="character", default=NULL, 
              help="output name")
)

print("Parsing arugments...")
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
infile = opt$i
gene_id = opt$gene_id
window_size = opt$window_size
buffer = opt$buffer
outname = opt$o

print('Reading in GFF...')
gff = read_tsv(infile,
               col_names = c('chr', 'source', 'feature', 'start', 'end', '.', 'strand', '..', 'description'))

print('Getting maximum lengths for chromosomes...')
#get maximum lengths
max_length_df = gff %>% 
  group_by(chr) %>% 
  summarize(max_length = max(end))
max_length = max_length_df$max_length
names(max_length) = max_length_df$chr

#filter to only genes
gff = gff %>% 
  filter(feature=='gene')


# HANDLE THE FORWARD GENES ------------------------------------------------

#if the set of boundaries is too short, just return an empty vector
window_from_bounds_for = function(chr, gene, bounds, region){
  if (length(bounds) > 1){
    l = bounds[1:(length(bounds)-1)]
    r = bounds[2:length(bounds)]
    w = data.frame(chr,l,r)
    w$name = paste(gene, region, sep='_')
  }
  else {
    w = c()
  }
  return(w)
}


#for forward genes
forward_gene_windows = function(gff_df){
  res = data.frame()
  for (i in 1:nrow(gff_df)){
    row = gff_df[i,]
    start = as.numeric(row$start)
    end = as.numeric(row$end)
    chr = as.character(row$chr)
    description = as.character(row$description)
    gene0 = strsplit(description, gene_id)[[1]][2]
    gene1 = strsplit(gene0, ';')[[1]][1]
    gene = sub('=', '', gene1, fixed=TRUE)
    max_len = max_length[chr]
    buffer_left = start - buffer
    if (buffer_left < 0){
      buffer_left = 0
    }
    buffer_right = end + buffer
    if (buffer_right < 0){
      buffer_right = 0
    }
    left_bounds = rev(seq(start, buffer_left, by = -window_size))
    inside_bounds = seq(start, end, by = window_size)
    inside_bounds[length(inside_bounds)] = end #add leftover into final window
    right_bounds = seq(end, buffer_right, by = window_size)
    left_windows = window_from_bounds_for(chr, gene, left_bounds, 'upstream')
    inside_windows = window_from_bounds_for(chr, gene, inside_bounds, 'inside')
    right_windows = window_from_bounds_for(chr, gene, right_bounds, 'downstream')
    res0 = rbind(left_windows, inside_windows, right_windows)
    dist = res0$l - start
    tss=paste('tss',start, sep='')
    tts = paste('tts', end, sep='')
    window_nums = paste('w', 1:nrow(res0), sep='')
    res0$name = paste(res0$name, tss, tts, sep='_')
    res = rbind(res, res0)
  }
  return(res)
}

print('Buiding windows for forward genes...')
forward_bed = gff %>% 
  filter(strand == '+') %>% 
  forward_gene_windows()


# HANDLE REVERSE READS ----------------------------------------------------

#function to turn boundaries into windows
window_from_bounds_rev = function(chr, gene, bounds, region){
  if (length(bounds) > 1){
    r = bounds[1:(length(bounds)-1)]
    l = bounds[2:length(bounds)]
    w = data.frame(chr,l,r)
    w$name = paste(gene, region, sep='_')
  } else{
    w=c()
  }
  return(w)
}

#for forward genes
reverse_gene_windows = function(gff_df){
  res = data.frame()
  for (i in 1:nrow(gff_df)){
    row = gff_df[i,]
    start = as.numeric(row$start)
    end = as.numeric(row$end)
    chr = as.character(row$chr)
    description = as.character(row$description)
    gene0 = strsplit(description, gene_id)[[1]][2]
    gene1 = strsplit(gene0, ';')[[1]][1]
    gene = sub('=', '', gene1, fixed=TRUE)
    max_len = max_length[chr]
    buffer_left = start - buffer
    if (buffer_left < 0){
      buffer_left = 0
    }
    buffer_right = end + buffer
    if (buffer_right < 0){
      buffer_right = 0
    }
    left_bounds = seq(start, buffer_left, by = -window_size) #NOTE DIFFERENCE FROM forward_gene_windows
    inside_bounds = seq(end, start, by = -window_size) #NOTE DIFFERENCE FROM forward_gene_windows
    inside_bounds[length(inside_bounds)] = start #include any leftover
    right_bounds = rev(seq(end, buffer_right, by = window_size)) #NOTE DIFFERENCE FROM forward_gene_windows
    left_windows = window_from_bounds_rev(chr, gene, left_bounds, 'downstream')
    inside_windows = window_from_bounds_rev(chr, gene, inside_bounds, 'inside')
    right_windows = window_from_bounds_rev(chr, gene, right_bounds, 'upstream')
    res0 = rbind(right_windows, inside_windows, left_windows)
    dist = end - res0$r
    tss=paste('tss',end, sep='')
    tts = paste('tts', start, sep='')
    window_nums = paste('w', 1:nrow(res0), sep='')
    res0$name = paste(res0$name, tss, tts, sep='_')
    res = rbind(res, res0)
  }
  return(res)
}

print('Buiding windows for reverse genes...')
reverse_bed = gff %>% 
  filter(strand == '-') %>% 
  reverse_gene_windows()




# output ------------------------------------------------------------------

print(paste('Saving results as ', outname, '...', sep=''))
bed = rbind(forward_bed, reverse_bed)
bed %>% 
  write_tsv(outname, col_names = FALSE)
