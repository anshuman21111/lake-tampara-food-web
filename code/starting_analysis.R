###############################################################################
###
### Starting analysis
###
###############################################################################
rm(list = ls())

websT <- c("Premonsoon","Postmonsoon","Monsoon")
locT  <- "data/Tampara/"

websL <- c("BSQ", "CHB","CPP","CSF","CSM","EMB","EPB","EWB","LRL","MPC","MRM","PCR","PRV","STM","YTH")
locL  <- "data/Lauras/"

read_NCP_edges <- function(webname){
  df <- read.csv(paste0(locL,"Final_", webname, "_NCP_edges.csv")) %>% 
    mutate(web = webname)
}

NCPs_dfL <- map(websL, read_NCP_edges) %>%
  list_rbind() %>%
  rename(from = "ResourceSpeciesID", to = "ConsumerSpeciesID") 

NCPs_dfT <- read.csv("data/Tampara/NCP_edges.csv") %>%
  pivot_longer(cols = c("PREMONSOON":"POSTMONSOON"), names_to = "web", values_to = "inweb") %>%
  filter(inweb == 1) %>%
  mutate(web = str_to_title(web)) %>%
  select(from, to, web)

NCP_providers <- rbind(NCPs_dfL[,c("from","to","web")], NCPs_dfT[,c("from","to","web")]) %>%
  mutate(role = "NCP_provider") %>%
  rename(species = "from", NCP = "to") %>%
  distinct()
NCPs_df <- rbind(NCPs_dfL[,c("to","web")], NCPs_dfT[,c("to","web")]) %>% 
  mutate(role = "NCP") %>%
  rename(species = "to") %>%
  mutate(NCP = species) %>%
  distinct()
roles <- rbind(NCP_providers, NCPs_df)
                      

count_web_motifs <- function(motif_list){
  counts <- motif_list %>%
    group_by(web, motif) %>%
    summarize(count = n())
  return(counts)
}  

count_species_motifs <- function(motif_list, s){
  subs <- motif_list[which(motif_list$sp1 == s | motif_list$sp2 == s | motif_list$sp3 == s ),]
  counts <- count_web_motifs(subs) %>%
    mutate(species = s)
  return(counts)
}


  
motifs_T <- read.csv(paste0(locT, "processed/motif_lists_NCP.csv"))
motifs_L <- read.csv(paste0(locL, "processed/motif_lists_NCP.csv"))  %>%
  mutate(web = sub("Final_", "", web))
#motifs <- rbind(motifs_T, motifs_L)
motifs <- motifs_L

N <- paste(unique(roles$NCP), collapse = "|")

NCP_motifs <- motifs %>%
  filter(grepl(N, sp1) | grepl(N, sp2) |grepl(N, sp3))
sp_motifs <-  motifs %>%
  filter(!(grepl(N, sp1) | grepl(N, sp2) |grepl(N, sp3)))


web_motif_counts <- count_web_motifs(sp_motifs) %>%
  group_by(web) %>%
  mutate(z = count/sum(count))

species <- unique(unlist(sp_motifs[,1:3]))


m_sp_counts <- map(species, function(s) count_species_motifs(sp_motifs,s)) %>%
  list_rbind() %>%
  group_by(web) %>%
  mutate(z = count/sum(count)) %>%
  mutate(species = as.character(species))
m_sp_counts2 <- left_join(m_sp_counts, roles, multiple = "all", by = c("web","species")) 
m_sp_counts2[which(is.na(m_sp_counts2$role)), "role"] <- "non-provider"
m_sp_counts2[which(is.na(m_sp_counts2$NCP)), "NCP"] <- "non-provider"

NCPs <- unique(NCPs_df$species)

m_NCP_counts <- map(NCPs, function(s) count_species_motifs(NCP_motifs,s)) %>%
  list_rbind() %>%
  group_by(web) %>%
  mutate(z = count/sum(count)) %>%
  mutate(species = as.character(species))
m_NCP_counts2 <- left_join(m_NCP_counts, roles, multiple = "all", by = c("web","species")) 
m_NCP_counts2[which(is.na(m_NCP_counts2$role)), "role"] <- "non-provider"
m_NCP_counts2[which(is.na(m_NCP_counts2$NCP)), "NCP"] <- "non-provider"


m_NCP_counts2 <- m_NCP_counts2 %>%
  mutate(species= case_when(species== 600 ~ "Wave attenuation",
                        species== 700 ~ "Shorline stabilization",
                        species== 800 ~ "Carbon sequestration",
                        species== 900 ~ "Water filtration",
                        species== 1000 ~ "Commerial fishery",
                        species== 1100 ~ "Bird watching",
                        species== 1200 ~ "Waterfowl hunting",
                        species== 1300 ~ "Recreational fishery",
                        .default = as.character(species)))

m_counts <- rbind(m_NCP_counts2, m_sp_counts2)

# How do motif frequencies compare across webs? (no NCP-containing motifs)
ggplot(data=web_motif_counts[which(web_motif_counts$motif %in% c("S1","S2","S4","S5")),], aes(x = web, y = z, fill = web)) +
  geom_col(position = "dodge") + facet_wrap(.~motif)


labs <- web_motif_counts %>% group_by(motif) %>% summarize(avg = median(count), sd = sd(count)) %>% mutate(label = paste0(motif, " (", signif(avg, digits = 2), " ± ", signif(sd, digits = 2), ")"))

motif_labeller <- function(variable,value){
  return(labs[which(labs$motif == value), "label"])
}

# How do motif counts vary across NCPs, NCP providers (not counting motifs with NCPs), and non-providers?
ggplot(data = m_counts, aes(x = role, y = count, color = role)) +
  geom_boxplot() +
  #geom_jitter(aes(color = web), alpha = 0.2) +
  facet_wrap(.~motif, scales = "free_y", labeller=motif_labeller)


# What about frequency?
ggplot(data = m_counts, aes(x = role, y = z, color = role)) +
  geom_boxplot() +
  #geom_jitter(aes(color = web), alpha = 0.2) +
  facet_wrap(.~motif, scales = "free_y", labeller=motif_labeller) +
  title(main = "Motif frequency across NCPs, NCP-providers, and non-providers", sub = "Only NCP frequencies include NCP-containing motifs")


m_sp_roles_summary <- m_sp_counts2 %>% 
  group_by(motif, NCP) %>%
  summarize(avg = median(z), low = quantile(z, 0.25), high = quantile(z, 0.75))


# What about kind of NCP?
ggplot(data = m_sp_roles_summary, aes(x = NCP, color = NCP)) +#m_sp_counts2, aes(x = z, y = NCP, color = NCP)) +
  #geom_boxplot() +
  geom_point(aes(y = avg)) +
  geom_errorbar(aes(ymin = low, ymax = high)) +
  #geom_density_ridges() +
  #geom_jitter(aes(color = web), alpha = 0.2) +
  facet_wrap(.~motif, scales = "free_x", labeller=motif_labeller) + 
  #geom_text(data = web_motif_counts, x =1, y =0, color = "black", aes(label = count)) + 
  theme_pubr()+
  coord_flip()
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



ggplot(data = m_sp_counts2, aes(color = role, x = z)) +
  #geom_boxplot() +
  #geom_jitter(aes(color = web), alpha = 0.2) +
  facet_wrap(.~motif, scales = "free")

library(ggpubr)
library(ggridges)
ggdensity(data = m_sp_counts2, aes(fill = role, y = z)) +
  facet_wrap(.~motif, scales = "free")
ggplot(data = m_sp_counts2, aes(y = role, x = count)) +
  geom_density_ridges() +
  #geom_jitter(aes(color = web), alpha = 0.2) +
  facet_wrap(.~motif, scales = "free")


spT <- read.csv("data/Tampara/processed/species_metrics.csv")
spL <- read.csv("data/Lauras/processed/species_metrics.csv")
sp_metrics <- rbind(spT, spL) %>%
  pivot_longer(cols = c("degree","indegree","outdegree","omnivory","TL"), names_to = "metric", values_to = "value")

sp_metrics_means <- sp_metrics %>%
  group_by(web, metric) %>%
  summarize(medianval = median(value), sd_low = quantile(value, 0.25), sd_high = quantile(value, 0.75))

ggplot(sp_metrics_means[which(sp_metrics_means$medianval < 100000000),], aes(x = web)) +
  geom_point(aes(y = medianval, fill = web)) +
  geom_errorbar(aes(ymin = sd_low, ymax = sd_high, color = web)) +
  facet_wrap(.~metric, scales = "free_y")# +
  coord_cartesian(ylim=c(0,100))
