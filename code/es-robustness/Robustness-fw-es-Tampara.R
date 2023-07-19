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
## disable warnings 
options(warn=-1)

# Generate and save sequence information

# Rarity sequence 
#rare.seq <- read.csv("Tampara_rare_seq.csv", header=T) # rarity sequence
# TO DO -- look at data to see how this was formatted 

# # Function for FW Robustness
# fw_robustness <- function (MATRIX, BASAL, N_RAND, OUTPUT, SEQ_NAM) {
#   # TO DO -- need to see these functions 
#   # TO DO -- double check this is the same for all the sequences used
#   mat_output <- collapse_wrap(N = MATRIX, basal = BASAL, n_rand = N_RAND)
#   mat_results <- subset(mat_output, mat_output$id == OUTPUT)
#   # Y AXIS
#   mat_results$nodes_susc <- tot.susc
#   mat_results$nontarget_lost_each <- stri_count_fixed(mat_results$name_lost,";") 
#   mat_results$cum_nontarget_lost <- cumsum(mat_results$nontarget_lost_each)
#   mat_results$Y <- (mat_results$nodes_susc - mat_results$cum_nontarget_lost) / mat_results$nodes_susc
#   # X AXIS
#   target_remove <- nrow(species_nodes)
#   mat_results$prop_target_removed <- mat_results$num_removed_tot / target_remove
#   # plot to visualize
#   root <- ggplot(mat_results, aes(x=prop_target_removed, y=Y))
#   (root + geom_point(shape=19) 
#     + geom_line() + theme_bw(base_size = 14) 
#     + xlim(0,1) + ylim(0,1)
#     + xlab("Proportion of target species removed")
#     + ylab("Proportion of susceptible species remaining")
#   )
#   # TO DO -- save plot with sequence name as title 
#   # calc robustness (AUC)
#   robust_auc(x = mat_results$prop_target_removed, y = mat_results$Y) 
#   write.csv(mat_results, sprintf("Tampara_fw_%s.csv",SEQ_NAME))
# }


# # Function for ES Robustness
# es_robustness <- function(MATRIX, BASAL, N_RAND, OUTPUT, NUM_SERVICES, SEQ_NAM) {
#   # Most to least connected extinction sequence
#   # TO DO -- need to see these functions to understand output 
#   # TO DO -- double check this is the same for all the sequences used 
#   mat_output_ES <- collapse_wrap(N = MATRIX, basal = BASAL, n_rand = 1)
#   mat_results_ES <- subset(mat_output_ES, mat_output_ES$degree_type == OUTPUT)
#   # make a column to track number of services lost at each step
#   mat_results_ES$service_lost <- ifelse(str_count(string = mat_results_ES$name_lost, pattern = "TO DO FILL THIS IN FOR TAMPARA"), str_count(mat_results_ES$name_lost,".50") , NA)
#   # track number services remain to use in propES_remain column
#   mat_results_ES[["service_lost"]][is.na(mat_results_ES[["service_lost"]])] <- 0
#   mat_results_ES$cum_nontarget_lost <- cumsum(mat_results_ES$service_lost)
#   mat_results_ES$ES_remain <- c(NUM_SERVICES-cumsum(mat_results_ES$service_lost))
#   mat_results_ES$propES_remain <- c(mat_results_ES$ES_remain/NUM_SERVICES)
#   # ADJUST X AXIS
#   target_remove <- nrow(species_nodes)
#   mat_results_ES$prop_removed <- mat_results_ES$num_removed_tot / target_remove
#   # revise output to get ES AUC - secondary loss of ES instead of services
#   mat_results_ES %>%
#     # add new column that tracks whenever a service node is lost
#     # input the service node names sep by a vertical bar: "|" ("or")
#     ## NOTE: cannot have any spaces between the species names and the |
#     mutate(service_lost = str_count(string = name_lost, pattern = "TO DO FILL IN")) %>%
#     # add column with the cumulative sum of the services lost
#     mutate(service_lost_c = cumsum(service_lost)) %>%
#     # overwrite original column of proportion species remaining, for now
#     mutate(prop_remain = 1 - (service_lost_c / length(service_nodes))) %>%
#     auc_wrapper()
#   root <- ggplot(mat_results_ES, aes(x=prop_removed, y=propES_remain))
#   (root + geom_point(shape=19) 
#     + geom_line() + theme_bw(base_size = 14) 
#     + xlim(0,1) + ylim(0,1)
#     + xlab("Proportion of target species removed")
#     + ylab("Proportion of ecosystem services remaining")
#   )
#   
#   write.csv(mat_results_ES, "Tampara_es_%s.csv")
# }

# TO DO -- clean this up.. 
# fw_robustness_random_1000 <- function(){
#   rand <- data.frame(matrix(ncol=3,nrow=1000)) # data frame to hold random results
#   x <- c("Randomization","R_web","R_ES")
#   colnames(rand) <- x
#   rand$Randomization <- seq(1:1000)
#   n <- c(1:1000) # save 100 random sequences
#   sequences <- vector(mode="list",length=1000)
#   for(i in n){
#     sequences[[i]] <- sample(colnames(mat.fw))
#   }
#   # create empty list to store robustness output for 1000 randomizations
#   mat_output_list <- vector(mode="list",length=1000)
#   for(i in n) {
#     mat_output_list[[i]] <- collapse_wrap_given(N = mat.spp, basal = sequences[[i]])
#   }
#   
#   mat_results_list <- mat_output_list
#   ## need to add columns to calc new REMAIN proportions
#   for( i in seq_along(mat_results_list)){
#     mat_results_list[[i]]$nodes_susc <- tot.susc
#   }
#   for( i in seq_along(mat_results_list)){
#     mat_results_list[[i]]$nontarget_lost_each <- stri_count_fixed(mat_results_list[[i]]$name_lost,";") #count num species LOST
#   }
#   for( i in seq_along(mat_results_list)){
#     mat_results_list[[i]]$cum_nontarget_lost <- cumsum(mat_results_list[[i]]$nontarget_lost_each) # calculate the cumulative prop nontarget LOST
#   }
#   
#   for( i in seq_along(mat_results_list)){
#     mat_results_list[[i]]$Y <- (mat_results_list[[i]]$nodes_susc - mat_results_list[[i]]$cum_nontarget_lost)/mat_results_list[[i]]$nodes_susc
#   }
#   
#   ##### calculate and store robustness for each of the randomizations
#   for (i in n){
#     rand$R_web[i] <- robust_auc(x = mat_results_list[[i]]$prop_removed, y = mat_results_list[[i]]$Y)
#   }
#   
#   min(rand$R_web)
#   mean(rand$R_web)
#   max(rand$R_web)
#   # TO DO -- SAVE THIS RESULT 
#   
# }

# # TO DO -- clean this up 
# es_robustness_random_1000 <- function(){
#   
#   # create empty list to store robustness output for 1000 randomizations
#   mat_output_list_ES <- vector(mode="list",length=1000)
#   for(i in n) {
#     mat_output_list_ES[[i]] <- collapse_wrap_given(N = mat.ES, basal = sequences[[i]])
#   }
#   
#   # calculate and store AUC for each of the randomizations
#   for(i in n){
#     rand$ES_AUC[i] <- mat_output_list_ES[[i]] %>%
#       # add new column that tracks whenever a service node is lost
#       # input the service node names sep by a vertical bar: "|" ("or")
#       ## NOTE: cannot have any spaces between the species names and the |
#       mutate(service_lost = str_count(string = name_lost, pattern = "350|450|550|650|750|850|950")) %>%
#       # add column with the cumulative sum of the services lost
#       mutate(service_lost_c = cumsum(service_lost)) %>%
#       # overwrite original column of proportion species remaining, for now
#       mutate(prop_remain = 1 - (service_lost_c / length(service_nodes))) %>%
#       auc_wrapper_given()
#   }
#   
#   
#   # unlist AUC then calculate robustness and store
#   rand$R_ES <- unlist(rand$ES_AUC)
#   rand <- rand[,-4]
#   
#   min(rand$R_ES)
#   mean(rand$R_ES)
#   max(rand$R_ES)
#   # TO DO -- SAVE THIS RESULT 
# 
# }
#   

seasons <- c("Premonsoon","Postmonsoon","Monsoon")
# Save as a quick check for loading data
quick_check <- data.frame('season'=character(),'nodes'=numeric(),'spec_nodes'=numeric(),'es_nodes'=numeric(),'edges'=numeric(),'fw_edges'=numeric(),'es_edges'=numeric(),'basal_spec'=numeric(),'susceptible_spec'=numeric(),stringsAsFactors=FALSE)

for (season in seasons){
  
  edges <- read.csv(sprintf("../../data/%s_edgelist_NCP.csv", season)) # edgelist
  edges[1] <- NULL # remove first column
  fw_edges <- subset(edges, edges$link_type == 0) # food web edges only
  es_edges <- subset(edges, edges$link_type == 1) # ES edges 
  species_nodes <- unique(c(unique(fw_edges$from),unique(fw_edges$to)))
  es_nodes <- unique(es_edges$to)
  
  NUM_SERVICES <- length(es_nodes)
  
  # TO DO - check from - to 
  fw <- graph.data.frame(fw_edges, directed = T)
  fw_es <- graph.data.frame(edges, directed = T)
  
  # Remove cannibalism 
  fw <- simplify(fw, remove.loops = T)
  fw_es <- simplify(fw_es, remove.loops=T)
  
  # convert graph objects to adjacency matrices to be used in source functions by AB
  mat.fw <- get.adjacency(fw, sparse=FALSE, attr = NULL)
  mat.fw_es <- get.adjacency(fw_es, sparse = FALSE, attr = NULL)
  
  # number of resources for each species 
  resources <- data.frame(species = species_nodes, InDegree = igraph::degree(fw, mode="in"))
  susc.spp <- resources[resources$InDegree>0,] # subset to the nonbasal, susceptible species
  basal.spp <- resources[resources$InDegree==0,] # subset to the basal, NOT susceptible species
  tot.susc <- nrow(susc.spp) # denominator for all sequences on food web y axis
  
  quick_check[nrow(quick_check)+1,] <- c(season,length(c(species_nodes,es_nodes)),length(species_nodes),length(es_nodes),length(c(fw_edges$from,es_edges$from)),length(fw_edges$from),length(es_edges$from),length(basal.spp$species),length(susc.spp$species))
  
  # FOOD WEB ROBUSTNESS 
  # TO DO -- need to see these functions to understand output 
  # Most to least connected extinction sequence
  #fw_robustness(mat.fw, colnames(mat.fw), 1, "high", "most_least")
  # Rarity sequence
  #fw_robustness(mat.fw, rarity_sequence, NULL, NULL, "rarity")
  # Random sequences
  #fw_robustness_random_1000()
  
  # ES ROBUSTNESS
  # TO DO -- need to see these functions to understand output
  # Most to least connected extinction sequence
  #es_robustness(mat.fw_es, colnames(mat.fw), 1, "high", NUM_SERVICES, "most_least")
  # Rarity sequence
  #es_robustness(mat.fw_es, rarity_sequence, NULL, NULL, NUM_SERVICES,  "rarity")
  # Random sequences
  #es_robustness_random_1000()

}



















