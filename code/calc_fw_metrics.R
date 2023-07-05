#####################################
###
###
###
#####################################

## Load packages
library(magrittr)
library(tidyr)
library(dplyr)
library(cheddar) # package that has food web analyses
library(igraph)


edges_prem <- read.csv("data/Tampara/Premonsoon_edgelist.csv", row.names = 1)

nodes_prem <- edges_prem %>%
  unlist() %>%
  unique
nodes_prem_df <- data.frame(node = nodes_prem)
properties_prem = list(title = "X")
#write.csv(nodes_prem_df, file = "data/Tampara/Premonsoon_nodes.csv", row.names = FALSE, quote = FALSE)


## Analyse for general food web metrics
community_prem <- Community(nodes = nodes_prem_df, trophic.links = edges_prem, properties = properties_prem)


## Get node-level metrics

node_properties_prem <- NPS(community_prem, c(Deg ='Degree', 
                                              InDeg = 'InDegree', 
                                              OutDeg = 'OutDegree', 
                                              NormInDeg = 'NormalisedTrophicGenerality',
                                              NormOutDeg = 'NormalisedTrophicVulnerability',
                                              Top ='IsTopLevelNode',
                                              Omn = 'IsOmnivore', 
                                              TS ='TrophicSpecies', 
                                              PATL ='PreyAveragedTrophicLevel'))

## Get web-level metrics
web_properties_prem <- CPS(community_prem, c(DegDist = "DegreeDistribution",
                                             S='NumberOfNodes', 
                                             C='DirectedConnectance', 
                                             'L/S' = 'LinkageDensity',
                                             Omn   = 'FractionOmnivorous'
                                             #'#'='NumberOfNodesByClass'#, 
                                             )
                           )


# Plot web
PlotWebByLevel(community_prem)
                                              