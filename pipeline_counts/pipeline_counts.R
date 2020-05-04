#pipeline_counts.R

library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
rc = read_tsv('pipeline_counts/raw_read_counts.txt',
              col_names=c('Run', 'value', 'stat')) %>% 
  mutate(dat = substr(Run, start=1, stop=3))


rc %>% 
  ggplot(aes(x=value)) +
  geom_histogram() +
  labs(x='raw read count') +
  facet_grid(~dat)
