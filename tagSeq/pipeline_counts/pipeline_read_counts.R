
library(tidyverse)
library(cowplot)
sdat = read_tsv('tagSeq/pipeline_counts/pipelineCounts.tsv', col_names=c('run', 'value', 'stat')) %>% 
  mutate(run = sub('_dupsRemoved.bam', '', run),
         run = sub('.', '-', run, fixed = TRUE))
unique(sdat$stat)
sdat$stat = factor(sdat$stat, levels=unique(sdat$stat))


#plot scatter abs
lp<-sdat %>% 
  mutate(value=as.numeric(value)) %>% 
  ggplot(aes(x=stat, y=value, color=run)) +
  geom_point() +
  geom_line(aes(group=run)) +
  theme(legend.position='none') +
  labs(y='Read count', x='Pipeline step', subtitle='read counts')
  

#plot scatter prop raw
pp<-sdat %>% 
  mutate(value=as.numeric(value)) %>% 
  group_by(run) %>% 
  mutate(prop = value/max(value) ) %>% 
  ggplot(aes(x=stat, y=prop, color=run)) +
    geom_point() +
    geom_line(aes(group=run)) +
    labs(y='Proportion raw reads', x='Pipeline step', subtitle='read proportions') +
    theme(legend.position='none')

plot_grid(lp, pp)
