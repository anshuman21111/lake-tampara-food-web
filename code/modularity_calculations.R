library(tidyverse)
library(rnetcarto)

modularity_list <- list()
for(season in c("Monsoon","Premonsoon","Postmonsoon")){
  adjmat <- read.csv(paste0("data/",season,"_adjacency_matrix.csv"), row.names = 1) %>% 
    as.matrix()
  rownames(adjmat) <- rownames(adjmat) %>% str_replace_all(" ", ".")
  modularity_list[[season]] <- netcarto(adjmat)
}
saveRDS(modularity_list, "data/modularity_list.rds")
