#####################################
###
###   Convert edge lists to matrices
###
#####################################

NCPs <- read.csv("data/Tampara/Ecosystem Service list.csv") %>% 
  separate_longer_delim(cols = "DIRECT_SERVICE_PROVIDER", delim = ";") %>% 
  rename(from = SPECIES, to = DIRECT_SERVICE_PROVIDER) %>%
  filter(to != "") %>%
  mutate(to = dplyr::recode(to, '1' = "Fishery",
                     '2' = "Bird_watching",
                     '3' = "Carbon sequestration",
                     '4' = "Nutrient recycling",
                     '5' = "Aquatic recreation",
                     '6' = "Food provisioning",
                     '7' = "Biocontrol",
                     '8' = "Bioindicators"))

write.csv(NCPs, "data/Tampara/NCP_edges.csv", row.names = FALSE, quote = FALSE)

## Add a column for link type (0 = consumptive, 1 = positive effect, -1 = negative effect)
NCPs$link_type = 1

species = unique(NCPs$from)
services = unique(NCPs$to)

### Add reverse links for fishing
#NCPs_fishing <- NCPs %>% 
#  filter(to == 1) %>%
#  mutate(link_type = -1) %>%
#  rename(from = "to", to = "from") %>%
#  rbind(NCPs)







### Add some random links to give some more motifs (to be replaced by real data later!!!)
#from = c(sample(species, 10, replace = TRUE), sample(services, 10, replace = TRUE))
#to = c(sample(services, 10, replace = TRUE), sample(species, 10, replace = TRUE))
#link_type = sample(c(-1,1), 20, replace = TRUE)
#metadat_ind = map_int(c(from[1:10], to[11:20]), function(x) detect_index(NCPs$from, function(y) y == x))
#metadat = NCPs[metadat_ind,] %>%
#  select(!c(from,to,link_type))
#  
#random_links = data.frame(from = from, to = to, link_type = link_type) %>%
#  filter(from != to) %>%
#  cbind(metadat)
#
#NCP_random = NCPs_fishing %>%
#  rbind(random_links)

for(web in c("Premonsoon","Postmonsoon","Monsoon")){
  ## Read in data
  edges <- read.csv(paste0("data/Tampara/",web,"_edgelist.csv"), row.names = 1)
  colnames(edges) <- c("from","to") # Needed for the cheddar package analyses (e.g. Community())
  edges$link_type = 0
  
  # Select services occurring in this web
  NCP_sub <- NCPs[which(NCPs[[str_to_upper(web)]] == 1),] %>%
    select("from", "to", "link_type")
  
  edges <- rbind(edges,NCP_sub)
  
  write.csv(edges, file = paste0("data/Tampara/",web,"_edgelist_NCP.csv"),quote = FALSE,row.names = TRUE)
  
  # Convert to adjacency matrix (for use in Julia)
  nodes <- sort(unique(c(edges$from, edges$to)))
  g <- graph.empty(length(nodes))
  V(g)$name <- nodes
  
  for (i in seq_len(nrow(edges))) {
    g <- add.edges(g, c(edges$from[i], edges$to[i]))#, weight = edges$link_type[i])
  }
  
  adj <- as.matrix(get.adjacency(g))#, attr = "weight"))
  
  write.csv(adj, file = paste0("data/Tampara/",web,"_adjacency_matrix_NCP.csv"),quote = FALSE,row.names = TRUE)
  
}

