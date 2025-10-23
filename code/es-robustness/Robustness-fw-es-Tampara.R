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
#library(patchwork)
#library(cowplot)
library(ggpubr)
library(svglite)

## AB source functions for analysis----
source("NEWrobustness_functions.R")
## Other helper functions for calling robustness functions and plotting
source("robustness_helper.R")
## disable warnings 
options(warn=-1)

NUM_RANDOM = 1000

# for consistent results from this script with randomization
set.seed(20)

seasons <- c("PREMONSOON","POSTMONSOON","MONSOON")
seasons_label <- c("Pre-monsoon","Post-monsoon","Monsoon")
net_props <- data.frame('season'=character(),
                        'spec_nodes'=numeric(),
                        'es_nodes'=numeric(),
                        'fw_edges'=numeric(),
                        'es_edges'=numeric(),
                        'basal_spec'=numeric(),
                        'susceptible_spec'=numeric(),
                        'R_fw_most_least'=numeric(),
                        'R_fw_random_mean'=numeric(),
                        'R_fw_random_sd'=numeric(),
                        'R_es_most_least'=numeric(),
                        'R_es_random_mean'=numeric(),
                        'R_es_random_sd'=numeric(),
                        stringsAsFactors=FALSE)
fw_most_least_all <- list()
es_most_least_all <- list()
fw_rand_all <- list()
es_rand_all <- list()

# Ecosystem service links
es_data_file <- read.csv("../../data/Ecosystem\ Service\ List.csv")
es_data_file <- es_data_file %>% filter(!is.na(SL)) # Filter out extra all NA rows

# Food web links
fw_data_file <- read.csv("../../data/Trophic\ interaction\ list.csv")
fw_data_file <- fw_data_file %>% filter(!is.na(SL))
# Look up table species ID to species name
species_id_map = list()
for (i in 1:length(fw_data_file$SL)){
  species_id_map[[as.character(fw_data_file$SL[i])]] <- as.character(fw_data_file$SPECIES[i])
}

# Comparing the columns that should be the same between these two files.
stopifnot(identical(es_data_file$SL, fw_data_file$SL))
stopifnot(identical(es_data_file$SPECIES, fw_data_file$SPECIES))
stopifnot(identical(es_data_file$PREMONSOON, fw_data_file$PREMONSOON))
stopifnot(identical(es_data_file$MONSOON, fw_data_file$MONSOON))
stopifnot(identical(es_data_file$POSTMONSOON, fw_data_file$POSTMONSOON))

for (season in seasons){
  
  fw_name <- sprintf("Tampara_Lake_%s", season)
  
  # ES edges from base file
  es_data_filtered <- filter(es_data_file, !!sym(season)==1)
  es_species_nodes <- es_data_filtered$SPECIES
  es_col <- es_data_filtered$DIRECT_SERVICE_PROVIDER
  
  es_edges_from <- list()
  es_edges_to <- list()
  for (i in 1:length(es_species_nodes)){
    if (!(gsub(" ","",es_col[i])=="")){ # if there are es
      es_dir <- unlist(strsplit(es_col[i], ";"))
      for (j in 1:length(es_dir)){
        es_edges_from[[length(es_edges_from) + 1]] <- es_species_nodes[i]
        es_edges_to[[length(es_edges_to) + 1]] <- as.integer(es_dir[j])
      }
    }
  }
  es_edges <- data.frame(unlist(es_edges_from), unlist(es_edges_to), rep.int(1,length(es_edges_from)))
  names(es_edges) <- c("from", "to", "link_type")
  # Remove any duplicates
  es_edges_unique <- es_edges %>% distinct()
  es_nodes <- unique(es_edges$to)
  
  # FW edges from base file
  fw_data_filtered <- filter(fw_data_file, !!sym(season)==1)
  fw_species_nodes <- fw_data_filtered$SPECIES
  stopifnot(identical(es_species_nodes, fw_species_nodes))
  prey_col <- fw_data_filtered$PREY
  
  fw_edges_from <- list()
  fw_edges_to <- list()
  for (i in 1:length(fw_species_nodes)){
    if (!(gsub(" ","", prey_col[i])=="")){ # if there is prey
      prey_list <- unlist(strsplit(prey_col[i], ";"))
      for (j in 1:length(prey_list)){
        prey_name <- as.character(species_id_map[as.character(as.integer(prey_list[j]))])
        if (prey_name %in% fw_species_nodes){
          # Point the edge from resource to consumer
          fw_edges_from[[length(fw_edges_from) + 1]] <- prey_name
          fw_edges_to[[length(fw_edges_to) + 1]] <- fw_species_nodes[[i]]
        } else{
          if (prey_name=="NULL"){
            print(prey_list[j])
          }
        }
      }
    }
  }
  fw_edges <- data.frame(unlist(fw_edges_from), unlist(fw_edges_to), rep.int(0,length(fw_edges_from)))
  names(fw_edges) <- c("from", "to", "link_type")
  fw_edges <- fw_edges %>% distinct()
  species_nodes <- unique(c(unique(fw_edges$from),unique(fw_edges$to)))
  
  full_edges <- rbind(fw_edges,es_edges)
  
  fw <- graph_from_data_frame(fw_edges, directed = TRUE) # network with only food web edges
  fw_es <- graph_from_data_frame(full_edges, directed = TRUE) # network with both food web and es edges 
  
  # Remove cannibalism 
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
  fw_most_least <- fw_robustness_high(mat.fw, species_nodes, fw_name, "most_least", tot.susc)
  auc_fw_most_least <- fw_most_least["auc"]
  fw_most_least_all[[length(fw_most_least_all) + 1]] <- fw_most_least["res"]
  es_most_least <- es_robustness(mat.fw_es, species_nodes, 1, "high", fw_name, "most_least", tot.susc, FALSE, es_nodes)
  auc_es_most_least <- es_most_least["auc"]
  es_most_least_all[[length(es_most_least_all) + 1]] <- es_most_least["res"]
  
  # Random sequences
  auc_res <- fw_robustness_random(mat.fw, species_nodes, NUM_RANDOM, tot.susc, fw_name)
  auc_fw_random_mean <- auc_res[1]
  auc_fw_random_sd <- auc_res[2]
  fw_rand_all[[length(fw_rand_all) + 1]] <- auc_res["res"]
  auc_res <- es_robustness_random(mat.fw_es, species_nodes, NUM_RANDOM, tot.susc, fw_name, es_nodes)
  auc_es_random_mean <- auc_res[1]
  auc_es_random_sd <- auc_res[2]
  es_rand_all[[length(es_rand_all) + 1]] <- auc_res["res"]
  
  net_props[nrow(net_props)+1,] <- c(season,
                                     length(species_nodes),
                                     length(es_nodes),
                                     length(fw_edges$from),
                                     length(es_edges$from),
                                     length(basal.spp$species),
                                     length(susc.spp$species),
                                     auc_fw_most_least,
                                     auc_fw_random_mean,
                                     auc_fw_random_sd,
                                     auc_es_most_least,
                                     auc_es_random_mean,
                                     auc_es_random_sd)

}

write.csv(net_props, "ES_robustness_results.csv")

# Output the results objects for easy access later to change plots
saveRDS(fw_most_least_all, file = "fw_most_least_all.RDS")
saveRDS(es_most_least_all, file = "es_most_least_all.RDS")
saveRDS(fw_rand_all, file = "fw_rand_all.RDS") 
saveRDS(es_rand_all, file = "es_rand_all.RDS") 

fw_rand_all <- readRDS("fw_rand_all.RDS")
fw_most_least_all <- readRDS("fw_most_least_all.RDS")
es_rand_all <- readRDS("es_rand_all.RDS")
es_most_least_all <- readRDS("es_most_least_all.RDS")

season_cols <- c("#05339C", "#40A2E3", "#41A67E")

fw_most_least_combo <- data.frame()
for (i in 1:length(seasons)){
  fw_most_least_all[[i]]["res"][[1]]$Season <- seasons_label[i]
  fw_most_least_combo <- rbind(fw_most_least_combo, fw_most_least_all[[i]]["res"][[1]])
}

# Subplot - FW most to least connected robustness
root <- ggplot(fw_most_least_combo, aes(x=prop_target_removed, y=Y, colour=Season))
p1 <- (root + geom_point(shape=19)
  + geom_line() + theme_bw(base_size = 14)
  + scale_colour_manual(values=season_cols)
  + xlim(0,1) + ylim(0,1)
  + theme(legend.position = "none",
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          plot.tag=element_text(size=14),
          plot.title=element_text(size=14),
          plot.tag.position=c(0.02,0.98))
  + labs(tag="A", title="Targeted")
)

es_most_least_combo <- data.frame()
for (i in 1:length(seasons)){
  es_most_least_all[[i]]["res"][[1]]$Season <- seasons_label[i]
  es_most_least_combo <- rbind(es_most_least_combo, es_most_least_all[[i]]["res"][[1]])
}

root <- ggplot(es_most_least_combo, aes(x=prop_removed, y=propES_remain, colour=Season))
p2 <- (root + geom_point(shape=19)
  + geom_line() + theme_bw(base_size = 14)
  + scale_colour_manual(values=season_cols)
  + xlim(0,1) + ylim(0,1)
  + theme(legend.position = "none",
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          plot.tag=element_text(size=14),
          plot.title=element_text(size=14),
          plot.tag.position=c(0.02,0.98))
  + labs(tag="B", title="Targeted")
)

fw_rand_combo <- data.frame()
argmin_auc <- c(NULL,NULL,NULL)
argmax_auc <- c(NULL,NULL,NULL)
mean_targets <- net_props$R_fw_random_mean
argmean_auc <- c(NULL,NULL,NULL)
ct <- 0
for (i in 1:length(seasons)){
  max_auc <- 0
  min_auc <- 100
  mean_diff <- 100
  for (j in 1:20){#NUM_RANDOM){
    fw_rand_all[[i]]$res[[j]]$Season <- seasons_label[i]
    fw_rand_all[[i]]$res[[j]]$line_group <- ct
    fw_rand_combo <- rbind(fw_rand_combo, fw_rand_all[[i]]$res[[j]])
    curr_auc <- robust_auc(x = fw_rand_all[[i]]$res[[j]]$prop_removed, y = fw_rand_all[[i]]$res[[j]]$Y)
    if (curr_auc < min_auc){
      min_auc <- curr_auc
      argmin_auc[i] <- ct
    }
    if (curr_auc > max_auc){
      max_auc <- curr_auc
      argmax_auc[i] <- ct
    }
    if (mean_diff > abs(curr_auc-mean_targets[i])){
      mean_diff <- abs(curr_auc-mean_targets[i])
      argmean_auc[i] <- ct
    }
    ct <- ct+1
  }
}
fw_rand_combo_1 <- filter(fw_rand_combo,line_group==argmin_auc[1]&Season==seasons_label[1])
fw_rand_combo_1b <- filter(fw_rand_combo,line_group==argmean_auc[1]&Season==seasons_label[1])
fw_rand_combo_1c <- filter(fw_rand_combo,line_group==argmax_auc[1]&Season==seasons_label[1])

fw_rand_combo_2 <- filter(fw_rand_combo,line_group==argmin_auc[2]&Season==seasons_label[2])
fw_rand_combo_2b <- filter(fw_rand_combo,line_group==argmean_auc[2]&Season==seasons_label[2])
fw_rand_combo_2c <- filter(fw_rand_combo,line_group==argmax_auc[3]&Season==seasons_label[3])

fw_rand_combo_3 <- filter(fw_rand_combo,line_group==argmin_auc[3]&Season==seasons_label[3])
fw_rand_combo_3b <- filter(fw_rand_combo,line_group==argmean_auc[3]&Season==seasons_label[3])
fw_rand_combo_3c <- filter(fw_rand_combo,line_group==argmax_auc[3]&Season==seasons_label[3])

fw_rand_combo <- rbind(fw_rand_combo_1,
                       fw_rand_combo_1b,
                       fw_rand_combo_1c,
                       fw_rand_combo_2,
                       fw_rand_combo_2b,
                       fw_rand_combo_2c,
                       fw_rand_combo_3,
                       fw_rand_combo_3b,
                       fw_rand_combo_3c)
fw_rand_combo$line_group <- as.factor(fw_rand_combo$line_group)


root <- ggplot(fw_rand_combo, aes(x=prop_removed, y=Y, colour=Season, 
                                  linetype = line_group,
                                  group = interaction(Season, line_group)))
p3 <- (root
  + geom_line(linewidth=0.6, alpha=0.9) + theme_bw(base_size = 14)
  + scale_colour_manual(values=season_cols)
  + xlim(0,1) + ylim(0,1)
  + theme(legend.position = "none",
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          plot.tag=element_text(size=14),
          plot.title=element_text(size=14),
          plot.tag.position=c(0.02,0.98))
  + labs(tag="C", title="Random")
)


es_rand_combo <- data.frame()
argmin_auc <- c(NULL,NULL,NULL)
argmax_auc <- c(NULL,NULL,NULL)
mean_targets <- net_props$R_es_random_mean
argmean_auc <- c(NULL,NULL,NULL)
ct <- 0
for (i in 1:length(seasons)){
  max_auc <- 0
  min_auc <- 100
  mean_diff <- 100
  for (j in 1:30){#NUM_RANDOM){
    es_rand_all[[i]]$res[[j]]$Season <- seasons_label[i]
    es_rand_all[[i]]$res[[j]]$line_group <- ct
    es_rand_combo <- rbind(es_rand_combo, es_rand_all[[i]]$res[[j]])
    curr_auc <- robust_auc(x = es_rand_all[[i]]$res[[j]]$prop_removed, y = es_rand_all[[i]]$res[[j]]$prop_remain)
    if (curr_auc < min_auc){
      min_auc <- curr_auc
      argmin_auc[i] <- ct
    }
    if (curr_auc > max_auc){
      max_auc <- curr_auc
      argmax_auc[i] <- ct
    }
    if (mean_diff > abs(curr_auc-mean_targets[i])){
      mean_diff <- abs(curr_auc-mean_targets[i])
      argmean_auc[i] <- ct
    }
    ct <- ct+1
  }
}
es_rand_combo_1 <- filter(es_rand_combo,((line_group==argmin_auc[1])|(line_group==argmax_auc[1])|(line_group==argmean_auc[1]))&Season==seasons_label[1])
es_rand_combo_2 <- filter(es_rand_combo,((line_group==argmin_auc[2])|(line_group==argmax_auc[2])|(line_group==argmean_auc[2]))&Season==seasons_label[2])
es_rand_combo_3 <- filter(es_rand_combo,((line_group==argmin_auc[3])|(line_group==argmax_auc[3])|(line_group==argmean_auc[3]))&Season==seasons_label[3])
es_rand_combo <- rbind(es_rand_combo_1,es_rand_combo_2,es_rand_combo_3)
es_rand_combo$line_group <- as.factor(es_rand_combo$line_group)

root <- ggplot(es_rand_combo, aes(x=prop_removed, y=prop_remain, colour=Season,
                                  linetype = line_group,
                                  group = interaction(Season, line_group)))
p4 <- (root
  + geom_line(linewidth=0.6, alpha=0.9) + theme_bw(base_size = 14)
  + scale_colour_manual(values=season_cols)
  + xlim(0,1) + ylim(0,1)
  + theme(legend.position = "none",
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          plot.tag=element_text(size=14),
          plot.title=element_text(size=14),
          plot.tag.position=c(0.02,0.98))
  + labs(tag="D", title="Random")
)

pp1 <- ggarrange(p1, p3, ncol=1, nrow=2)
pp1 <- annotate_figure(pp1,left = text_grob("Proportion of susceptible species remaining", rot = 90, size = 14))
pp2 <- ggarrange(p2, p4, ncol=1, nrow=2)
pp2 <- annotate_figure(pp2,left = text_grob("Proportion of NCP remaining", rot = 90, size = 14))

combo <- ggarrange(pp1, pp2, ncol=2, nrow=1)
combo <- annotate_figure(combo,
                bottom = text_grob("Proportion of species removed", size = 14))
ggsave("Robustness_Fig1.svg", combo, width = 8, height = 6, device = "svg")

# With common legend (to move into plot after)
combo2 <- ggarrange(p1, p2, p3, p4, ncol=2, nrow=2, common.legend=TRUE, legend="bottom")
ggsave("Robustness_Fig2.svg", combo2, width = 8, height = 6, device = "svg")




