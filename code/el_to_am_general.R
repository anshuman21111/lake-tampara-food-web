#####################################
###
###   Convert edge lists to matrices
###
#####################################

library(tidyverse)
library(igraph)

webs <- c("BSQ", "CHB","CPP","CSF","CSM","EMB","EPB","EWB","LRL","MPC","MRM","PCR","PRV","STM","YTH")
location <- "data/Lauras/"
prefix <- "Final_"
edge_suffix <- "_edges.csv"

#webs <- c("Premonsoon","Monsoon","Postmonsoon")
#location <- "data/Tampara/"
#prefix <- ""
#edge_suffix <- "edgelist"


NCP_suffix <- ""#"_NCP" # "" (Do we want to include NCPs or not?)

for(web in webs){
  ## Read in data
  edges <- read.csv(paste0(location, prefix, web, edge_suffix)) %>%
    select(!starts_with("X")) 
  edges[,1] <- as.character(edges[,1])
  edges[,2] <- as.character(edges[,2])
  
  
  if(NCP_suffix == ""){
    NCP_edges <- edges %>% 
      filter(Type == "ES")
    write.csv(NCP_edges, file = paste0(location,prefix,web,"_NCP_edges.csv"),quote = FALSE,row.names = TRUE)
    edges <- edges %>%
      filter(Type != "ES")
  }
  
  #colnames(edges) <- c("from","to", "link_type")
  g <- graph_from_edgelist(as.matrix(edges[,1:2]), directed = TRUE)
  
  adj <- as.matrix(get.adjacency(g))
  
  write.csv(adj, file = paste0(location,prefix,web,"_adjacency_matrix",NCP_suffix,".csv"),quote = FALSE,row.names = TRUE)
  
  
}

