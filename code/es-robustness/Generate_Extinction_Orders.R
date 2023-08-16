library(tidyverse)

# Data frames with species details
species_df <- read.csv("../../data/species_data_frame.csv", sep=';')
es_df <- read.csv("../../data/Ecosystem\ Service\ list.csv")
species_abundances <- as.list(es_df$REL_ABUN_CAT) # dictionary
names(species_abundances) <- es_df$SPECIES

# Create sequences for each food web
seasons <- c("Premonsoon","Postmonsoon","Monsoon")
for (season in seasons){
  edges <- read.csv(sprintf("../../data/%s_edgelist_NCP.csv", season)) # edgelist
  edges[1] <- NULL # remove first column
  fw_edges <- subset(edges, edges$link_type == 0) # food web edges only
  species_nodes <- unique(c(unique(fw_edges$from),unique(fw_edges$to)))
  curr_abundances <- unlist(lapply(species_nodes, function(x) {return(species_abundances[[x]])}))
  
  # sort nodes by abundance category and create rarity sequence 
  species_w_abundances <- cbind.data.frame(SPECIES=species_nodes,REL_ABUN_CAT=curr_abundances)
  species_w_abundances <- species_w_abundances[order(species_w_abundances$REL_ABUN_CAT),]
  # TO DO - create a bunch of random versions, shuffling within each abundance category
  rarity_sequence <- cbind.data.frame(ORDER=species_w_abundances$SPECIES)
  write.csv(rarity_sequence, sprintf("%s_rarity_sequence.csv",season))
}