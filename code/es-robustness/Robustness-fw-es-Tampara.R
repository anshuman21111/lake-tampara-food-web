##########################################
######## Food web and ES robustness calculations 
### Written by: Aislyn Keyes using functions written by Allison Barner
### From Nature Communications 2021, 10.1038/s41467-021-21824-x

### Modified for use in Tampara Lake project by Lucy Van Kleunen
### July 2023 

# load packages
library(tidyverse)
library(igraph)
library(dplyr)
library(stringi)
library(tibble)

## AB source functions for analysis----
source("NEWrobustness_functions.R")
## Other helper functions for calling robustness functions and plotting
source("robustness_helper.R")
## disable warnings 
options(warn=-1)

# for consistent results from this script with randomization
set.seed(20)

seasons <- c("Premonsoon","Postmonsoon","Monsoon")
net_props <- data.frame('season'=character(),'nodes'=numeric(),'spec_nodes'=numeric(),'es_nodes'=numeric(),'edges'=numeric(),'fw_edges'=numeric(),'es_edges'=numeric(),'basal_spec'=numeric(),'susceptible_spec'=numeric(),'R_fw_most_least'=numeric(),'R_es_most_least'=numeric(),'R_fw_rarity'=numeric(),'R_es_rarity'=numeric(),'R_fw_random_min'=numeric(),'R_fw_random_mean'=numeric(),'R_fw_random_max'=numeric(),'R_es_random_min'=numeric(),'R_es_random_mean'=numeric(),'R_es_random_max'=numeric(),stringsAsFactors=FALSE)

for (season in seasons){
  
  fw_name <- sprintf("Tampara_Lake_%s", season)
  edges <- read.csv(sprintf("../../data/%s_edgelist_NCP.csv", season)) # edgelist
  edges[1] <- NULL # remove first column
  fw_edges <- subset(edges, edges$link_type == 0) # food web edges only
  species_nodes <- unique(c(unique(fw_edges$from),unique(fw_edges$to)))
  es_edges <- subset(edges, edges$link_type == 1) # ES edges only 
  # remove any es edges that include species not in the species list 
  es_edges <- subset(es_edges, es_edges$from %in% species_nodes) 
  # full set of edges to make the es web
  es_nodes <- unique(es_edges$to)
  full_edges <- rbind(fw_edges,es_edges)
  
  fw <- graph_from_data_frame(fw_edges, directed = TRUE) # network with only food web edges
  fw_es <- graph_from_data_frame(full_edges, directed = TRUE) # network with both food web and es edges 
  
  # Remove cannibalism - TO DO - confirm we want to do this. 
  fw <- simplify(fw, remove.loops = TRUE)
  fw_es <- simplify(fw_es, remove.loops=TRUE)
  
  # convert graph objects to adjacency matrices to be used in source functions by AB
  mat.fw <- as_adjacency_matrix(fw, sparse= FALSE, attr = NULL)
  mat.fw_es <- as_adjacency_matrix(fw_es, sparse = FALSE, attr = NULL)
  
  # number of resources for each species 
  resources <- data.frame(species = species_nodes, InDegree = igraph::degree(fw, v= species_nodes, mode="in"))
  susc.spp <- resources[resources$InDegree>0,] # subset to the non-basal, susceptible species
  basal.spp <- resources[resources$InDegree==0,] # subset to the basal, NOT susceptible species
  tot.susc <- nrow(susc.spp) # denominator for all sequences on food web y axis
  
  # Most to least connected extinction sequence
  auc_fw_most_least <- fw_robustness(mat.fw, species_nodes, 1, "high", fw_name, "most_least", tot.susc, FALSE)
  auc_es_most_least <- es_robustness(mat.fw_es, species_nodes, 1, "high", fw_name, "most_least", tot.susc, FALSE, es_nodes)
  
  # Rarity sequence 
  rarity_sequence <- read.csv(sprintf("%s_rarity_sequence.csv",season), header=T)
  rarity_sequence <- rarity_sequence$ORDER
  # double check same number of species
  if (length(rarity_sequence) != length(species_nodes)){
    print("rarity sequence not the correct length")
  }
  auc_fw_rarity <- fw_robustness(mat.fw, rarity_sequence, NULL, NULL, fw_name, "rarity", tot.susc, TRUE)
  auc_es_rarity <- es_robustness(mat.fw_es, rarity_sequence, NULL, NULL, fw_name, "rarity", tot.susc, TRUE, es_nodes)
  
  # Random sequences
  auc_res <- fw_robustness_random(mat.fw, species_nodes, 1000, tot.susc, fw_name)
  auc_fw_random_min <- auc_res[1]
  auc_fw_random_mean <- auc_res[2]
  auc_fw_random_max <- auc_res[3]
  auc_res <- es_robustness_random(mat.fw_es, species_nodes, 1000, tot.susc, fw_name, es_nodes)
  auc_es_random_min <- auc_res[1]
  auc_es_random_mean <- auc_res[2]
  auc_es_random_max <- auc_res[3]
  
  net_props[nrow(net_props)+1,] <- c(season,length(c(species_nodes,es_nodes)),length(species_nodes),length(es_nodes),length(c(full_edges$from)),length(fw_edges$from),length(es_edges$from),length(basal.spp$species),length(susc.spp$species),auc_fw_most_least,auc_es_most_least,auc_fw_rarity,auc_es_rarity,auc_fw_random_min,auc_fw_random_mean,auc_fw_random_max,auc_es_random_min,auc_es_random_mean,auc_es_random_max)

}






