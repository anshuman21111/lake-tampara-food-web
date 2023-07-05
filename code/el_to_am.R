#####################################
###
###   Convert edge lists to matrices
###
#####################################


for(web in c("Premonsoon","Postmonsoon","Monsoon")){
  ## Read in data
  edges <- read.csv(paste0("data/Tampara/",web,"_edgelist.csv"), row.names = 1)
  colnames(edges) <- c("resource","consumer") # Needed for the cheddar package analyses (e.g. Community())
  
  # Convert to adjacency matrix (for use in Julia)
  nodes <- sort(unique(c(edges$resource, edges$consumer)))
  g <- graph.empty(length(nodes))
  V(g)$name <- nodes
  
  for (i in seq_len(nrow(edges))) {
    g <- add.edges(g, c(edges$resource[i], edges$consumer[i]))
  }
  
  adj <- as.matrix(get.adjacency(g))
  
  write.csv(adj, file = paste0("data/Tampara/",web,"_adjacency_matrix.csv"),quote = FALSE,row.names = TRUE)
  
}
