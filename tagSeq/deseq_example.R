#deseq_example.R


# load inputs for deseq ---------------------------------------------------
#(these were created in initialize_feature_counts.R)
ll = load('tagSeq/data/deseqInput.Rdata')
ll #look at the two object names that were loaded
head(counts) #the processesed counts matrix
head(coldata) #the table of sample traits

#to make manipulations easier, lets assign 



#note that we have lots of features that could influence the gene expression in these samples
#based on previous studies, and from the PCA, we expect colony to be an important source of gene expression variation
unique(coldata$colony) #we have 4 colonies in this experiment

#we also expect treatment to be important
unique(coldata$treatment)

#the age of the branch (whether it was from middle or outside of the colony) probably doesn't influence expression a lot
unique(coldata$age) #y=young; o=old

#the timepoint Temp (time during the data samples were taken, morning or afternoon hence cool or hot)
unique(coldata$timepointTemp) #note some of these have missing values, let's just remove these for now
coldata = coldata[!is.na(coldata$timepointTemp),]
counts = counts[,rownames(coldata)]
dim(coldata)
dim(counts)

#finally we had timepoints for the duration of the experiment
unique(coldata$timepointDescription)
#let's remove space from that just to be safe
for (i in seq(1,5,1)){
  coldata$timepointDescription = sub(' ', '_', coldata$timepointDescription)
}

# run deseq with very basic model ----------------------------------------------------
#first let's run deseq with a very simple model
#here we will ignore all the features in coldata except for timepointTemp
#thus we are telling DESEq to measure differential expression between samples that had "hot" or "cool" in the coldata table,
#without worrying about any of the other features.
#we code this by only including timepointTemp in the formula:
dds<-DESeqDataSetFromMatrix(counts,
                            colData = coldata,
                            design =  formula(~ timepointTemp))

#the dds object now has our counts, coldata, and a forumula associated with it
#now we can run DEseq on it
dds <- DESeq(dds)

#now we can extract the results from it
resultsNames(dds) #look at the results associated with DESeq object. "timepointTemp_hot_vs_cool" is the one we want
res = results(dds,
              contrast = c('timepointTemp', 'hot', 'cool'),
              independentFiltering = FALSE) #by designating the contrast like this we mean that positive log2 changes indicate upregulation in the "hot" samples


#now we can plot a volcano plot
#let's write a function to do it
plot_volcano = function(deseq_results, main_title=''){
  data.frame(deseq_results) %>% 
    mutate(significant = padj < 0.05,
           neg_log_pvalue = -log(pvalue, 10)) %>% 
    ggplot(aes(x=log2FoldChange, y=neg_log_pvalue, color = significant)) +
    geom_point() +
    scale_color_manual(values = c('black', 'red')) +
    labs(y=bquote('-log'[10]~'p-value'),
         x = bquote('log'[2]~'fold change'),
         title = main_title)
}

plot_volcano(res, main_title='time-of-day')
#we see here that we have a log of highly significant genes for this contrast

#we can also rescale to the Y axis to get a better look:
plot_volcano(res, main_title='time-of-day') + lims(y=c(0,15))


# run deseq with a better model -------------------------------------------
#now let's run deseq again using a more sophisticated model that includes multiple features
#this makes sense, because we know that biologically, more things than just timepointTemp should
#influence gene expression accross these samples

#set up a new deseq object with a full set of predictor features in the formula
dds<-DESeqDataSetFromMatrix(counts,
                            colData = coldata,
                            design =  formula(~ timepointTemp +
                                                treatment +
                                                timepointDescription +
                                                colony +
                                                age))


#run deseq for new formula
dds <- DESeq(dds)

#note now we have lots of results names to pick from
resultsNames(dds)

#we can only consider the log2 fold differences for one feature at a time,
#but we can create multiple results tables, one for each feature
indFilt = FALSE #set up a variable for independentFiltering in case we want to use it later (use ?results to find out what independentFiltering is about)

#get a new updated timepointTemp results table
#this one is better than the previous one, because it accounts for variation from the other features
time_hot_cool = results(dds,
                              contrast = c('timepointTemp', 'hot', 'cool'),
                              independentFiltering = indFilt)

#get results for the treatment (the heat or control bucket)
treatment_heat_control = results(dds,
                        contrast = c('treatment', 'heat', 'control'),
                        independentFiltering = indFilt)

#get results for the treatment (the heatSwitched or control bucket)
#note how this one is looking at the same feature 'treatment' but a different contrast (heatSwitched vs controls)
treatment_switched_control = results(dds,
                        contrast = c('treatment', 'heatSwitched', 'control'),
                        independentFiltering = indFilt)

#get the age
age_young_old = results(dds,
                      contrast = c('age', 'y', 'o'),
                      independentFiltering = indFilt)



#PLOT THE VOLCANO PLOT FOR NEW TIME-OF-DAY EFFECTS
plot_volcano(hotTime_vs_coolTime, main_title='new time-of-day') + lims(y=c(0,15))


#now we can get fancy and put all our results into a list and plot them together
#a list can hold lots of different types of things. Here we are using it to store 
#the deseq results tables.
results_list = list('time_hot_cool' = hotTime_vs_coolTime,
                    'treatment_heat_control' = treatment_heat_control,
                    'treatment_switched_control' = treatment_switched_control,
                    'age_young_old' = age_young_old)

#we named each thing in the list with a string matching the variable name
names(results_list)


#now we can loop through the list names and plot each one
for (results_name in names(results_list)){
  print(paste('plotting volcano plot for results table named',results_name))
  results_table = results_list[[results_name]]
  plt = plot_volcano(results_table, main_title=results_name) + lims(y=c(0,15))
  plot(plt)
}


#we can also put plots into a list
#run through the same loop, this time adding the plots to an empty list called my_volcano_plots
my_volcano_plots = list() #create the list
for (results_name in names(results_list)){
  results_table = results_list[[results_name]]
  plt = plot_volcano(results_table, main_title=results_name) + lims(y=c(0,15))
  my_volcano_plots[[results_name]] = plt #add the plot to my_volcano_plots
}


#now we can plot them in a multipanel plot
library(cowplot) 
plot_grid(plotlist = my_volcano_plots, nrow=2)