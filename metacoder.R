setwd("E:/Dropbox/Dropbox/##Soil Microbes Internship/Data/séquences/algae_23S")
getwd()


## Packages 
install.packages("metacoder")

library(tidyverse)
library(ape)
library(metacoder)


## Load Fasta file 

seqs <- read.FASTA()


## Parsing 

data <- extract_taxonomy(seqs,
                         regex = "",
                         key = c(id= "obs_info", "class"),
                         class_sep = ";")

## Manipulation and digital PCR 

pcr <- filter_taxa(data, n_obs >1) %>% 
  filter_obs(nchar(sequence) > 1000) %>%
  filter_taxa(name == "Bacteria", subtaxa = T) %>% 
  sample_n_obs(5000, taxon_weight = 1 / n_obs) %>% 
  sample_n_taxa(1000, supertaxa = T) %>% 
  primersearch(forward = "",
               reverse = "",
               mismatch = 10) %>% 
  filter_taxa(prop_amplified < 0.9, supertaxa = T)


## Visualisation

heat_tree(pcr,
          node_size = n_obs, 
          node_label = name,
          node_color = prop_amplified,
          node_color_range = c("red", "yellow", "cyan"))